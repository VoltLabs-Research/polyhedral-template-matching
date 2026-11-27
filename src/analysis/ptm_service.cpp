#include <volt/analysis/reconstructed_analysis_pipeline.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cluster_graph_builder.h>
#include <volt/analysis/orientation_cluster_rule_provider.h>
#include <volt/analysis/reconstructed_state_canonicalizer.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_identification_export.h>
#include <volt/analysis/polyhedral_template_matching_log.h>
#include <volt/analysis/ptm_service.h>
#include <volt/analysis/ptm.h>
#include <volt/analysis/ptm_structure_analysis.h>
#include <volt/analysis/ptm_crystal_info_provider.h>
#include <volt/matcher/template_matcher.h>

#include <volt/core/lammps_parser.h>
#include <volt/core/analysis_result.h>
#include <volt/core/particle_property.h>
#include <volt/core/progress.h>
#include <volt/core/type_counts.h>
#include <volt/cli/common.h>
#include <volt/utilities/json_utils.h>

#include <volt/structures/crystal_structure_types.h>

#include <cctype>
#include <string_view>
#include <utility>

#include <volt/utilities/parquet_atom_writer.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <future>
#include <map>

namespace Volt{

Matrix3 PolyhedralTemplateMatchingService::quaternionToMatrix(const Quaternion& orientation){
    const Quaternion normalized = orientation.normalized();
    return Matrix3(
        normalized * Vector3(1.0, 0.0, 0.0),
        normalized * Vector3(0.0, 1.0, 0.0),
        normalized * Vector3(0.0, 0.0, 1.0)
    );
}

std::shared_ptr<std::vector<OrientationClusterAtomState>> PolyhedralTemplateMatchingService::buildOrientationClusterStates(
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
    for(std::string token : CLI::splitList(families, ',')){
        for(char& c : token){
            c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
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

void PolyhedralTemplateMatchingService::setIdentifyOrdering(bool identifyOrdering){
    _identifyOrdering = identifyOrdering;
}

bool PolyhedralTemplateMatchingService::identifiesOrdering(const LammpsParser::Frame& frame, std::size_t atomCount) const{
    if(!_identifyOrdering){
        return false;
    }
    if(frame.types.size() < atomCount){
        spdlog::warn("PTM: ordering identification requires per-atom types; the input dump has none");
        return false;
    }
    return true;
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
        PolyhedralTemplateMatchingLog log;
        StructureAnalysis analysis(context);
        auto ptmAtomStates = std::make_shared<std::vector<PtmLocalAtomState>>();

        // The outputs weigh more than in a serial plugin: PTM runs in parallel, the neighbor
        // lattice export does not.
        ProgressSteps steps({
            {"Running PTM", 4.0},
            {"Building clusters", 2.0},
            {"Canonicalizing cluster symmetry", 1.0},
            {"Writing outputs", 3.0},
        });
        TemplateMatcher templates;
        if(!_latticeDirectory.empty()){
            log.dump(templates.loadDirectory(_latticeDirectory));
        }

        const TemplateMatcher* templatesPtr = templates.empty() ? nullptr : &templates;
        const bool identifyOrdering = identifiesOrdering(frame, context.atomCount());

        PtmRmsdCutoffReport rmsdCutoffReport;
        determineLocalStructuresWithPTM(analysis, _rmsd, ptmAtomStates, templatesPtr,
                                        _structureCheckFlags,
                                        identifyOrdering ? frame.types.data() : nullptr,
                                        &rmsdCutoffReport);
        log.dump(analysis);
        steps.next();
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

        if(context.inputCrystalType == LATTICE_SC){
            analysis.setClusterRuleProvider(
                std::make_shared<OrientationClusterRuleProvider>(buildOrientationClusterStates(ptmAtomStates))
            );
        }

        ClusterBuilder clusterBuilder(analysis, context);
        const BuiltClusters clusters = clusterBuilder.build(_dissolveSmallClusters);
        steps.next();
        normalizeReconstructedClusterGraphForExport(analysis, context);
        log.dump(clusters);
        steps.next();

        for(const auto& [atomIndex, type] : matchedTypes){
            context.structureTypes->setInt(atomIndex, type);
        }

        const std::map<int, int> structureCounts = countTypes(context.structureTypes->constIntRange());
        log.dump(structureCounts);

        json result = AnalysisResult::success();
        result["main_listing"] = {
            {"total_atoms", frame.natoms},
            {"structure_count", structureCounts.size()},
            {"rmsd", _rmsd},
            {"atoms_discarded_by_rmsd", rmsdCutoffReport.rejectedAtoms},
            {"smallest_discarded_rmsd", rmsdCutoffReport.smallestRejectedRmsd},
            {"median_discarded_rmsd", rmsdCutoffReport.medianRejectedRmsd}
        };

        result["sub_listings"]["structure_counts"] = StructureIdentificationExport::structureCountRows(
            structureCounts, context.atomCount()
        );

        if(identifyOrdering){
            std::map<int, std::int64_t> orderingCounts;
            for(std::size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
                orderingCounts[atomIndex < ptmAtomStates->size()
                    ? (*ptmAtomStates)[atomIndex].orderingType
                    : static_cast<int>(PTM::OrderingType::ORDERING_NONE)]++;
            }

            const double totalAtomsForFraction = context.atomCount() > 0
                ? static_cast<double>(context.atomCount()) : 1.0;
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
            log.dump(orderingCounts);
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
        log.dump(*ptmAtomStates);

        if(!outputBase.empty()){
            auto ptmColumnWriter = [&ptmAtomStates, identifyOrdering](ColumnarRowWriter& writer, std::size_t atomIndex, int){
                if(atomIndex >= ptmAtomStates->size()){
                    return;
                }
                const PtmLocalAtomState& state = (*ptmAtomStates)[atomIndex];
                if(!state.valid){
                    return;
                }
                writer.field("rmsd", state.rmsd);
                writer.field("orientation", state.orientation.normalized());
                writer.field("interatomic_distance", state.interatomicDistance);
                writer.field("correspondences", static_cast<std::int64_t>(state.correspondencesCode));
                if(identifyOrdering){
                    writer.field("ordering_type", state.orderingType);
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

            // Emscripten builds have no threads, so an async launch would never run.
#if defined(__EMSCRIPTEN__)
            constexpr auto atomsExportLaunch = std::launch::deferred;
#else
            constexpr auto atomsExportLaunch = std::launch::async;
#endif
            std::future<AtomsTable> atomsTable = std::async(
                atomsExportLaunch,
                [&, ptmColumnWriter, nameResolver]{
                    return StructureIdentificationExport::writeAtomsTable(
                        outputBase,
                        frame,
                        analysis,
                        nameResolver,
                        ptmColumnWriter
                    );
                }
            );

            const auto clusterOutputs = AnalysisPipelineUtils::appendClusterOutputs(
                frame,
                outputBase,
                annotatedDumpPath,
                context,
                analysis,
                result,
                &frameError,
                extraDumpColumns
            );
            if(!clusterOutputs){
                atomsTable.wait();
                return AnalysisResult::failure(frameError);
            }
            log.dump(clusterOutputs);
            log.dump(JsonUtils::writeSummaryTable(result, outputBase, "ptm_analysis"));
            log.dump(atomsTable.get());
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
