#include <volt/analysis/reconstructed_analysis_pipeline.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cluster_graph_builder.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/orientation_cluster_rule_provider.h>
#include <volt/analysis/reconstructed_state_canonicalizer.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/analysis/ptm_service.h>
#include <volt/analysis/ptm.h>
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

#include <cctype>
#include <string_view>
#include <utility>

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

PolyhedralTemplateMatchingService::PolyhedralTemplateMatchingService()
    : _inputCrystalStructure(LATTICE_FCC)
    , _rmsd(0.1)
    , _dissolveSmallClusters(false)
    , _latticeDirectory()
    , _cationNeighborRadius(0.0)
    , _cationMisorientation(12.0)
    , _structureCheckFlags(0)
    , _identifyOrdering(false){}

void PolyhedralTemplateMatchingService::setStructureCheckFlags(std::string families){
    _structureCheckFlags = 0;
    if(families.empty()){
        return;
    }

    static const std::pair<std::string_view, int> table[] = {
        {"SC", PTM_CHECK_SC},       {"FCC", PTM_CHECK_FCC},
        {"HCP", PTM_CHECK_HCP},     {"ICO", PTM_CHECK_ICO},
        {"BCC", PTM_CHECK_BCC},     {"DCUB", PTM_CHECK_DCUB},
        {"DHEX", PTM_CHECK_DHEX},   {"GRAPHENE", PTM_CHECK_GRAPHENE},
        {"ALL", PTM_CHECK_ALL},
    };

    int flags = 0;
    std::size_t start = 0;
    while(start <= families.size()){
        std::size_t comma = families.find(',', start);
        if(comma == std::string::npos){
            comma = families.size();
        }
        std::string token = families.substr(start, comma - start);
        token.erase(0, token.find_first_not_of(" \t"));
        const auto lastChar = token.find_last_not_of(" \t");
        token.erase(lastChar == std::string::npos ? 0 : lastChar + 1);
        for(char& c : token){
            c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }

        if(!token.empty()){
            bool matched = false;
            for(const auto& [name, bit] : table){
                if(token == name){
                    flags |= bit;
                    matched = true;
                    break;
                }
            }
            if(!matched){
                spdlog::warn("PTM: unknown structure family '{}' in --ptm_structures, ignoring", token);
            }
        }
        start = comma + 1;
    }

    _structureCheckFlags = flags;
}

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

void PolyhedralTemplateMatchingService::setIdentifyOrdering(bool identifyOrdering){
    _identifyOrdering = identifyOrdering;
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

        const bool identifyOrdering = _identifyOrdering
            && frame.types.size() >= context.atomCount();
        if(_identifyOrdering && !identifyOrdering){
            spdlog::warn("PTM: ordering identification requires per-atom types; the input dump has none");
        }

        determineLocalStructuresWithPTM(analysis, _rmsd, ptmAtomStates, templatesPtr, cationRadius,
                                        _structureCheckFlags,
                                        identifyOrdering ? frame.types.data() : nullptr);
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

        json structureCountRows = json::array();
        const double totalAtomsForFraction = context.atomCount() > 0
            ? static_cast<double>(context.atomCount()) : 1.0;
        for(const auto& [structureType, count] : structureCounts){
            structureCountRows.push_back({
                {"structure_id", structureType},
                {"structure_name", structureTypeName(structureType)},
                {"count", static_cast<int64_t>(count)},
                {"fraction", static_cast<double>(count) / totalAtomsForFraction}
            });
        }

        result["sub_listings"] = json::object();
        result["sub_listings"]["structure_counts"] = std::move(structureCountRows);

        if(identifyOrdering){
            std::map<int, std::int64_t> orderingCounts;
            for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
                orderingCounts[atomIndex < ptmAtomStates->size()
                    ? (*ptmAtomStates)[atomIndex].orderingType
                    : static_cast<int>(PTM::OrderingType::ORDERING_NONE)]++;
            }

            json orderingCountRows = json::array();
            for(const auto& [orderingType, count] : orderingCounts){
                orderingCountRows.push_back({
                    {"ordering_id", orderingType},
                    {"ordering_name", PTM::orderingTypeName(orderingType)},
                    {"count", count},
                    {"fraction", static_cast<double>(count) / totalAtomsForFraction}
                });
            }
            result["sub_listings"]["ordering_counts"] = std::move(orderingCountRows);
        }

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

            auto ptmColumnWriter = [&ptmAtomStates, identifyOrdering](ColumnarAtomWriter& writer, std::size_t atomIndex, int){
                if(atomIndex >= ptmAtomStates->size()){
                    return;
                }
                const PtmLocalAtomState& state = (*ptmAtomStates)[atomIndex];
                if(!state.valid){
                    return;
                }
                const Quaternion orientation = state.orientation.normalized();
                writer.field("rmsd", state.rmsd);
                writer.field("orientation", std::vector<double>{orientation.x(), orientation.y(), orientation.z(), orientation.w()});
                writer.field("interatomic_distance", state.interatomicDistance);
                writer.field("correspondences", static_cast<std::int64_t>(state.correspondencesCode));
                if(identifyOrdering){
                    writer.field("ordering_type", static_cast<std::int64_t>(state.orderingType));
                }
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
