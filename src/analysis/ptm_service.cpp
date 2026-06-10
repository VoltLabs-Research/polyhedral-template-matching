#include <volt/analysis/reconstructed_analysis_pipeline.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cluster_graph_builder.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/orientation_cluster_rule_provider.h>
#include <volt/analysis/reconstructed_state_canonicalizer.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/core/analysis_result.h>
#include <volt/analysis/ptm_service.h>
#include <volt/analysis/ptm_structure_analysis.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/json_utils.h>
#include <volt/utilities/parquet_atom_writer.h>

#include <algorithm>
#include <future>
#include <map>
#include <utility>

namespace Volt{

namespace{

Matrix3 quaternionToMatrix(const Quaternion& orientation){ // kept for orientation cluster rules
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
    , _dissolveSmallClusters(false){}

void PolyhedralTemplateMatchingService::setInputCrystalStructure(LatticeStructureType structureType){
    _inputCrystalStructure = structureType;
}

void PolyhedralTemplateMatchingService::setRMSD(double rmsd){
    _rmsd = rmsd;
}

void PolyhedralTemplateMatchingService::setDissolveSmallClusters(bool dissolveSmallClusters){
    _dissolveSmallClusters = dissolveSmallClusters;
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

        determineLocalStructuresWithPTM(analysis, _rmsd, ptmAtomStates);
        analysis.setClusterRuleProvider(nullptr);
        std::fill(
            context.atomSymmetryPermutations->dataInt(),
            context.atomSymmetryPermutations->dataInt() + context.atomSymmetryPermutations->size(),
            -1
        );

        const bool requiresScClusterRules = context.inputCrystalType == LATTICE_SC;
        if(requiresScClusterRules){
            analysis.setClusterRuleProvider(
                std::make_shared<OrientationClusterRuleProvider>(buildOrientationClusterStates(ptmAtomStates))
            );
        }

        ClusterBuilder clusterBuilder(analysis, context);
        clusterBuilder.build(_dissolveSmallClusters);
        normalizeReconstructedClusterGraphForExport(analysis, context);

        // Count structures for summary
        std::map<int,int> structCounts;
        for(std::size_t i = 0; i < context.atomCount(); ++i)
            structCounts[context.structureTypes ? context.structureTypes->getInt(i) : 0]++;

        json result = AnalysisResult::success();
        result["main_listing"] = {
            {"total_atoms", frame.natoms},
            {"structure_count", static_cast<int>(structCounts.size())},
            {"rmsd", _rmsd}
        };
        result["sub_listings"] = json::object();
        result["per-atom-properties"] = json::array();

        std::future<void> atomsExportFuture;

        if(!outputBase.empty()){
            const std::string analysisPath = outputBase + "_ptm_analysis.parquet";
            const std::string atomsPath = outputBase + "_atoms.parquet";
            // Streaming export with PTM per-atom columns
            auto ptmColumnWriter = [&ptmAtomStates](ColumnarAtomWriter& w, std::size_t atomIndex, int){
                if(atomIndex >= ptmAtomStates->size()) return;
                const PtmLocalAtomState& state = (*ptmAtomStates)[atomIndex];
                w.field("ptm_valid", state.valid);
                if(!state.valid) return;
                const Quaternion q = state.orientation.normalized();
                w.field("rmsd", state.rmsd);
                w.field("orientation", std::vector<double>{q.x(), q.y(), q.z(), q.w()});
                w.field("interatomic_distance", state.interatomicDistance);
                w.field("scaling", state.interatomicDistance);
                std::vector<double> grad(9);
                const Matrix3& F = state.deformationGradient;
                for(int r = 0; r < 3; ++r)
                    for(int c = 0; c < 3; ++c)
                        grad[r * 3 + c] = F(r, c);
                w.field("deformation_gradient", grad);
                w.field("ordering_type", state.orderingType);
                w.field("correspondences", static_cast<std::int64_t>(state.correspondencesCode));
                w.field("template_index", state.bestTemplateIndex);
            };

            atomsExportFuture = std::async(
                std::launch::async,
                [&, atomsPath, ptmColumnWriter]{
                    StructureIdentificationExport::streamStructureIdentificationToParquet(
                        atomsPath,
                        frame,
                        analysis,
                        nullptr,
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
                &frameError
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
            &frameError
        )){
            return AnalysisResult::failure(frameError);
        }

        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("PTM analysis failed: ") + error.what());
    }
}

}
