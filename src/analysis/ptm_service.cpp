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
#include <volt/analysis/ptm_cluster_input_adapter.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/json_utils.h>

#include <algorithm>
#include <map>
#include <utility>

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

std::string structureTypeNameForExport(int structureType){
    return structureTypeName(structureType);
}

json buildStructureListing(const AnalysisContext& context){
    std::map<int, int> counts;
    for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
        const int structureType = context.structureTypes
            ? context.structureTypes->getInt(atomIndex)
            : static_cast<int>(StructureType::OTHER);
        counts[structureType]++;
    }

    json listing = json::array();
    for(const auto& [structureType, count] : counts){
        if(count <= 0){
            continue;
        }
        listing.push_back({
            {"structure_id", structureType},
            {"structure_name", structureTypeNameForExport(structureType)},
            {"atom_count", count}
        });
    }

    std::sort(listing.begin(), listing.end(), [](const json& lhs, const json& rhs){
        return lhs.value("structure_name", "") < rhs.value("structure_name", "");
    });

    return listing;
}

json buildPerAtomProperties(
    const LammpsParser::Frame& frame,
    const StructureAnalysis& analysis,
    const std::vector<PtmLocalAtomState>& atomStates
){
    const StructureContext& context = analysis.context();
    json perAtom = json::array();

    for(std::size_t atomIndex = 0; atomIndex < static_cast<std::size_t>(frame.natoms); ++atomIndex){
        const int structureType = context.structureTypes
            ? context.structureTypes->getInt(atomIndex)
            : static_cast<int>(StructureType::OTHER);

        json atom;
        atom["id"] = atomIndex < frame.ids.size()
            ? frame.ids[atomIndex]
            : static_cast<int>(atomIndex);
        atom["structure_type"] = structureType;
        atom["structure_name"] = structureTypeNameForExport(structureType);
        atom["cluster_id"] = context.atomClusters ? context.atomClusters->getInt(atomIndex) : 0;

        if(atomIndex < frame.positions.size()){
            const auto& pos = frame.positions[atomIndex];
            atom["pos"] = {pos.x(), pos.y(), pos.z()};
        }else{
            atom["pos"] = {0.0, 0.0, 0.0};
        }

        if(atomIndex < atomStates.size()){
            const PtmLocalAtomState& state = atomStates[atomIndex];
            atom["ptm_valid"] = state.valid;
            if(state.valid){
                const Quaternion orientation = state.orientation.normalized();
                atom["rmsd"] = state.rmsd;
                atom["orientation"] = {orientation.x(), orientation.y(), orientation.z(), orientation.w()};
                atom["interatomic_distance"] = state.interatomicDistance;
                atom["scaling"] = state.interatomicDistance;
                const Matrix3& F = state.deformationGradient;
                atom["deformation_gradient"] = {
                    F(0, 0), F(0, 1), F(0, 2),
                    F(1, 0), F(1, 1), F(1, 2),
                    F(2, 0), F(2, 1), F(2, 2)
                };
                atom["ordering_type"] = state.orderingType;
                atom["correspondences"] = state.correspondencesCode;
                atom["template_index"] = state.bestTemplateIndex;
            }
        }

        perAtom.push_back(std::move(atom));
    }

    return perAtom;
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
        // Why: always collect per-atom PTM state so the `_atoms.msgpack`
        // exporter can surface RMSD / orientation / deformation gradient
        // parity with OVITO's PolyhedralTemplateMatchingModifier. Previously
        // the container was only allocated for SC input (which needs it to
        // drive cluster building); every other crystal structure silently
        // discarded PTM's per-atom output.
        auto ptmAtomStates = std::make_shared<std::vector<PtmLocalAtomState>>();

        determineLocalStructuresWithPTM(analysis, _rmsd, ptmAtomStates);
        computeMaximumNeighborDistanceFromPTM(analysis);
        PTMClusterInputAdapter clusterInputAdapter;
        clusterInputAdapter.prepare(analysis, context);

        const bool requiresScClusterRules = context.inputCrystalType == LATTICE_SC;
        if(requiresScClusterRules){
            analysis.setClusterRuleProvider(
                std::make_shared<OrientationClusterRuleProvider>(buildOrientationClusterStates(ptmAtomStates))
            );
        }

        ClusterBuilder clusterBuilder(analysis, context);
        clusterBuilder.build(_dissolveSmallClusters);
        normalizeReconstructedClusterGraphForExport(analysis, context);

        json result = AnalysisResult::success();
        const json structuresListing = buildStructureListing(context);
        result["main_listing"] = {
            {"total_atoms", frame.natoms},
            {"structure_count", static_cast<int>(structuresListing.size())},
            {"rmsd", _rmsd}
        };
        result["sub_listings"] = {
            {"structures", structuresListing}
        };
        result["per-atom-properties"] = buildPerAtomProperties(
            frame,
            analysis,
            *ptmAtomStates
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
            return AnalysisResult::failure(frameError);
        }

        if(!outputBase.empty()){
            const std::string analysisPath = outputBase + "_ptm_analysis.msgpack";
            if(!JsonUtils::writeJsonMsgpackToFile(result, analysisPath, false)){
                return AnalysisResult::failure("Failed to write " + analysisPath);
            }

            const std::string atomsPath = outputBase + "_atoms.msgpack";
            // Per-atom PTM augmentation — matches OVITO's PTM Modifier
            // output (rmsd, orientation, scaling, interatomic distance,
            // elastic deformation gradient, correspondences code).
            auto augmentAtom = [&ptmAtomStates](nlohmann::json& atom, std::size_t atomIndex, int /*structureType*/){
                if(atomIndex >= ptmAtomStates->size()){
                    return;
                }
                const PtmLocalAtomState& state = (*ptmAtomStates)[atomIndex];
                atom["ptm_valid"] = state.valid;
                if(!state.valid){
                    return;
                }
                const Quaternion orientation = state.orientation.normalized();
                atom["rmsd"] = state.rmsd;
                atom["orientation"] = {orientation.x(), orientation.y(), orientation.z(), orientation.w()};
                atom["interatomic_distance"] = state.interatomicDistance;
                // Scaling is the ratio between the template-fit interatomic
                // distance and the canonical spacing. OVITO exposes it as
                // `Scaling` in the user properties; we expose the same here.
                atom["scaling"] = state.interatomicDistance;
                const Matrix3& F = state.deformationGradient;
                atom["deformation_gradient"] = {
                    F(0, 0), F(0, 1), F(0, 2),
                    F(1, 0), F(1, 1), F(1, 2),
                    F(2, 0), F(2, 1), F(2, 2)
                };
                atom["ordering_type"] = state.orderingType;
                atom["correspondences"] = state.correspondencesCode;
                atom["template_index"] = state.bestTemplateIndex;
            };

            if(!JsonUtils::writeJsonMsgpackToFile(
                StructureIdentificationExport::buildStructureIdentificationJson(
                    frame,
                    analysis,
                    nullptr,
                    augmentAtom
                ),
                atomsPath,
                false
            )){
                return AnalysisResult::failure("Failed to write " + atomsPath);
            }
        }

        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("PTM analysis failed: ") + error.what());
    }
}

}
