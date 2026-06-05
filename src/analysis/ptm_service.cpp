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
#include <volt/analysis/cluster_input_preparation.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/json_utils.h>

#include <algorithm>
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
        computeMaximumNeighborDistanceFromPTM(analysis);
        ClusterInputAdapterUtils::prepareSymmetryAwareClusterInputs(
            analysis,
            context,
            false,
            [&](std::size_t atomIndex, int structureType){
                if(structureType == LATTICE_OTHER){
                    return false;
                }
                if(analysis.numberOfNeighbors(static_cast<int>(atomIndex)) == 0){
                    return false;
                }
                return context.atomAllowedSymmetryMasks->getInt64(atomIndex) == 0;
            }
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

        if(!AnalysisPipelineUtils::appendClusterOutputs(
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

        if(!outputBase.empty()){
            const std::string analysisPath = outputBase + "_ptm_analysis.msgpack";
            if(!JsonUtils::writeJsonMsgpackToFile(result, analysisPath, false)){
                return AnalysisResult::failure("Failed to write " + analysisPath);
            }

            const std::string atomsPath = outputBase + "_atoms.msgpack";
            // Streaming export with PTM per-atom fields
            auto ptmFieldWriter = [&ptmAtomStates](MsgpackWriter& w, std::size_t atomIndex, int, int& extraCount){
                if(atomIndex >= ptmAtomStates->size()){ extraCount = 0; return; }
                const PtmLocalAtomState& state = (*ptmAtomStates)[atomIndex];
                extraCount = state.valid ? 8 : 1;
                w.write_key("ptm_valid"); w.write_bool(state.valid);
                if(!state.valid) return;
                const Quaternion q = state.orientation.normalized();
                w.write_key("rmsd"); w.write_double(state.rmsd);
                w.write_key("orientation");
                w.write_array_header(4);
                w.write_double(q.x()); w.write_double(q.y()); w.write_double(q.z()); w.write_double(q.w());
                w.write_key("interatomic_distance"); w.write_double(state.interatomicDistance);
                w.write_key("scaling"); w.write_double(state.interatomicDistance);
                w.write_key("deformation_gradient");
                w.write_array_header(9);
                const Matrix3& F = state.deformationGradient;
                for(int r = 0; r < 3; ++r)
                    for(int c = 0; c < 3; ++c)
                        w.write_double(F(r, c));
                w.write_key("ordering_type"); w.write_int(state.orderingType);
                w.write_key("correspondences"); w.write_uint(state.correspondencesCode);
                w.write_key("template_index"); w.write_int(state.bestTemplateIndex);
            };

            StructureIdentificationExport::streamStructureIdentificationToFile(
                atomsPath, frame, analysis, nullptr, ptmFieldWriter
            );
        }

        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("PTM analysis failed: ") + error.what());
    }
}

}
