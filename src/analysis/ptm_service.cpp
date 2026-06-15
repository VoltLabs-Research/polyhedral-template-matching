#include <volt/analysis/reconstructed_analysis_pipeline.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cluster_graph_builder.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/orientation_cluster_rule_provider.h>
#include <volt/analysis/reconstructed_state_canonicalizer.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/analysis/ptm_service.h>
#include <volt/analysis/ptm_structure_analysis.h>
#include <volt/analysis/ptm_crystal_info_provider.h>
#include <volt/matcher/template_matcher.h>
#include <volt/matcher/template_crystal_info_provider.h>
#include <volt/cluster-rules/orientation_based.h>

#include <volt/core/lammps_parser.h>
#include <volt/core/analysis_result.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/particle_property.h>

#include <volt/structures/crystal_structure_types.h>

#include <volt/utilities/json_utils.h>
#include <volt/utilities/parquet_atom_writer.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <future>
#include <map>
#include <utility>
#include <fstream>
#include <sstream>

namespace Volt{

namespace{

Matrix3 quaternionToMatrix(const Quaternion& orientation){
    const Quaternion normalized = orientation.normalized();
    return Matrix3(
        normalized * Vector3(1.0, 0.0, 0.0),
        normalized * Vector3(0.0, 1.0, 0.0),
        normalized * Vector3(0.0, 0.0, 1.0)
    );
}

std::shared_ptr<std::vector<OrientationClusterAtomState>> buildOrientationClusterStates(
    const std::shared_ptr<const std::vector<PtmLocalAtomState>>& atomStates
){
    if(!atomStates){
        return nullptr;
    }

    auto states = std::make_shared<std::vector<OrientationClusterAtomState>>(atomStates->size());
    for(std::size_t atomIndex = 0; atomIndex < atomStates->size(); ++atomIndex){
        const PtmLocalAtomState& source = (*atomStates)[atomIndex];
        auto& target = (*states)[atomIndex];
        target.valid = source.valid;
        if(source.valid){
            target.orientation = quaternionToMatrix(source.orientation);
        }
    }
    return states;
}

}

PolyhedralTemplateMatchingService::PolyhedralTemplateMatchingService()
    : _inputCrystalStructure(LATTICE_FCC)
    , _rmsd(0.1)
    , _dissolveSmallClusters(false)
    , _latticeDirectory()
    , _cationNeighborRadius(0.0)
    , _cationMisorientation(12.0){}

void PolyhedralTemplateMatchingService::setInputCrystalStructure(LatticeStructureType structureType){
    _inputCrystalStructure = structureType;
}

void PolyhedralTemplateMatchingService::setRMSD(double rmsd){
    _rmsd = rmsd;
}

void PolyhedralTemplateMatchingService::setDissolveSmallClusters(bool dissolveSmallClusters){
    _dissolveSmallClusters = dissolveSmallClusters;
}

void PolyhedralTemplateMatchingService::setLatticesDirectory(std::string latticesDirectory){
    _latticeDirectory = std::move(latticesDirectory);
}

void PolyhedralTemplateMatchingService::setCationNeighborRadius(double radius){
    _cationNeighborRadius = radius;
}

void PolyhedralTemplateMatchingService::setCationMisorientation(double degrees){
    _cationMisorientation = degrees;
}

json PolyhedralTemplateMatchingService::compute(
    const LammpsParser::Frame& frame,
    const std::string& outputBase,
    const std::string& inputDumpPath
){
    const std::string annotatedDumpPath = outputBase.empty()
        ? inputDumpPath
        : outputBase + "_annotated.dump";

    std::string frameError;
    auto session = AnalysisPipelineUtils::prepareAnalysisSession(
        frame,
        _inputCrystalStructure,
        &frameError
    );
    if(!session){
        return AnalysisResult::failure(frameError);
    }
    AnalysisContext& context = session->context;

    try{
        StructureAnalysis analysis(context);
        auto ptmAtomStates = std::make_shared<std::vector<PtmLocalAtomState>>();

        TemplateMatcher templates;
        if(!_latticeDirectory.empty()){
            const int loaded = templates.loadDirectory(_latticeDirectory);
            if(loaded == 0){
                spdlog::warn("PTM: no user templates loaded from '{}'", _latticeDirectory);
            }
        }

        const TemplateMatcher* templatesPtr = templates.empty() ? nullptr : &templates;
        const bool clusterTemplates = (templatesPtr != nullptr && _cationNeighborRadius > 0.0);
        const double cationRadius = clusterTemplates ? _cationNeighborRadius : 0.0;

        determineLocalStructuresWithPTM(analysis, _rmsd, ptmAtomStates, templatesPtr, cationRadius);
        analysis.setClusterRuleProvider(nullptr);
        std::fill(
            context.atomSymmetryPermutations->dataInt(),
            context.atomSymmetryPermutations->dataInt() + context.atomSymmetryPermutations->size(),
            -1
        );

        std::vector<std::pair<std::size_t, int>> matchedTypes;
        if(templatesPtr){
            const std::size_t atomCount = context.atomCount();
            for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
                const int type = context.structureTypes->getInt(atomIndex);
                if(type >= TEMPLATE_STRUCTURE_TYPE_BASE){
                    matchedTypes.emplace_back(atomIndex, type);
                }
            }
        }

        if(clusterTemplates){
            analysis.setCrystalInfoProvider(
                std::make_shared<TemplateCrystalInfoProvider>(
                    PtmStructureAnalysisDetail::ptmCrystalInfoProvider(),
                    templates,
                    MAX_NEIGHBORS
                )
            );

            auto orientationStates = std::make_shared<std::vector<OrientationBasedAtomState>>(ptmAtomStates->size());
            for(std::size_t atomIndex = 0; atomIndex < ptmAtomStates->size(); ++atomIndex){
                const PtmLocalAtomState& source = (*ptmAtomStates)[atomIndex];
                (*orientationStates)[atomIndex].valid = source.valid;
                if(source.valid){
                    (*orientationStates)[atomIndex].orientation = quaternionToMatrix(source.orientation);
                }
            }

            analysis.setClusterRuleProvider(
                std::make_shared<OrientationBasedClusterRule>(
                    orientationStates,
                    context.neighborOffsets->constDataInt(),
                    context.neighborIndices ? context.neighborIndices->constDataInt() : nullptr,
                    context.structureTypes->constDataInt(),
                    context.atomCount(),
                    _cationMisorientation
                )
            );
        }else if(context.inputCrystalType == LATTICE_SC){
            analysis.setClusterRuleProvider(
                std::make_shared<OrientationClusterRuleProvider>(buildOrientationClusterStates(ptmAtomStates))
            );
        }

        ClusterBuilder clusterBuilder(analysis, context);
        clusterBuilder.build(_dissolveSmallClusters);
        normalizeReconstructedClusterGraphForExport(analysis, context);

        for(const auto& [atomIndex, type] : matchedTypes){
            context.structureTypes->setInt(atomIndex, type);
        }

        std::map<int, int> structureCounts;
        for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
            structureCounts[context.structureTypes ? context.structureTypes->getInt(atomIndex) : 0]++;
        }

        json result = AnalysisResult::success();
        result["main_listing"] = {
            {"total_atoms", frame.natoms},
            {"structure_count", static_cast<int>(structureCounts.size())},
            {"rmsd", _rmsd}
        };
        result["sub_listings"] = json::object();
        result["per-atom-properties"] = json::array();

        std::vector<AnalysisContext::ExtraScalarColumn> extraDumpColumns;
        {
            const std::size_t atomCount = context.atomCount();
            auto rmsdProperty = std::make_shared<ParticleProperty>(atomCount, DataType::Double, 1, 0, true);
            double* rmsdData = rmsdProperty->dataDouble();
            for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
                rmsdData[atomIndex] = (atomIndex < ptmAtomStates->size() && (*ptmAtomStates)[atomIndex].valid)
                    ? (*ptmAtomStates)[atomIndex].rmsd : -1.0;
            }
            extraDumpColumns.push_back({ "rmsd", rmsdProperty });
        }

        std::future<void> atomsExportFuture;

        if(!outputBase.empty()){
            const std::string analysisPath = outputBase + "_ptm_analysis.parquet";
            const std::string atomsPath = outputBase + "_atoms.parquet";

            auto ptmColumnWriter = [&ptmAtomStates](ColumnarAtomWriter& writer, std::size_t atomIndex, int){
                if(atomIndex >= ptmAtomStates->size()){
                    return;
                }
                const PtmLocalAtomState& state = (*ptmAtomStates)[atomIndex];
                writer.field("ptm_valid", state.valid);
                if(!state.valid){
                    return;
                }
                const Quaternion orientation = state.orientation.normalized();
                writer.field("rmsd", state.rmsd);
                writer.field("orientation", std::vector<double>{orientation.x(), orientation.y(), orientation.z(), orientation.w()});
                writer.field("interatomic_distance", state.interatomicDistance);
                writer.field("scaling", state.interatomicDistance);
                std::vector<double> deformationGradient(9);
                const Matrix3& gradient = state.deformationGradient;
                for(int row = 0; row < 3; ++row){
                    for(int column = 0; column < 3; ++column){
                        deformationGradient[row * 3 + column] = gradient(row, column);
                    }
                }
                writer.field("deformation_gradient", deformationGradient);
                writer.field("ordering_type", state.orderingType);
                writer.field("correspondences", static_cast<std::int64_t>(state.correspondencesCode));
                writer.field("template_index", state.bestTemplateIndex);
            };

            StructureIdentificationExport::StructureNameResolver nameResolver =
                [&templates](std::size_t, int structureType) -> std::string {
                    std::string templateName = templates.structureName(structureType);
                    if(!templateName.empty()){
                        return templateName;
                    }
                    return structureTypeName(structureType);
                };

            atomsExportFuture = std::async(
                std::launch::async,
                [&, atomsPath, ptmColumnWriter, nameResolver]{
                    StructureIdentificationExport::streamStructureIdentificationToParquet(
                        atomsPath,
                        frame,
                        analysis,
                        nameResolver,
                        ptmColumnWriter
                    );
                }
            );

            if(!AnalysisPipelineUtils::appendClusterOutputs(
                frame,
                outputBase,
                annotatedDumpPath,
                context,
                analysis,
                result,
                &frameError,
                extraDumpColumns
            )){
                atomsExportFuture.wait();
                return AnalysisResult::failure(frameError);
            }

            if(!JsonUtils::writeJsonToParquet(result, analysisPath, false)){
                atomsExportFuture.wait();
                return AnalysisResult::failure("Failed to write " + analysisPath);
            }

            atomsExportFuture.get();
        }else if(!AnalysisPipelineUtils::appendClusterOutputs(
            frame,
            outputBase,
            annotatedDumpPath,
            context,
            analysis,
            result,
            &frameError,
            extraDumpColumns
        )){
            return AnalysisResult::failure(frameError);
        }

        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("PTM analysis failed: ") + error.what());
    }
}

}
