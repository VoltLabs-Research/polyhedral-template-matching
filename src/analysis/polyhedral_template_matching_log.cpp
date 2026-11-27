#include <volt/analysis/polyhedral_template_matching_log.h>

#include <volt/analysis/ptm.h>
#include <volt/analysis/ptm_structure_analysis.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/matcher/template_matcher.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/core/type_counts.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <limits>
#include <utility>

namespace Volt{

std::string PolyhedralTemplateMatchingLog::structureFamilyNames(int checkFlags){
    static const std::pair<int, StructureType> families[] = {
        {PTM_CHECK_SC, StructureType::SC},
        {PTM_CHECK_FCC, StructureType::FCC},
        {PTM_CHECK_HCP, StructureType::HCP},
        {PTM_CHECK_ICO, StructureType::ICO},
        {PTM_CHECK_BCC, StructureType::BCC},
        {PTM_CHECK_DCUB, StructureType::CUBIC_DIAMOND},
        {PTM_CHECK_DHEX, StructureType::HEX_DIAMOND},
        {PTM_CHECK_GRAPHENE, StructureType::GRAPHENE},
    };

    std::string names;
    for(const auto& [flag, structureType] : families){
        if(checkFlags & flag){
            if(!names.empty()){
                names += ", ";
            }
            names += structureTypeName(structureType);
        }
    }
    return names.empty() ? std::string("none") : names;
}

void PolyhedralTemplateMatchingLog::dump(const TemplateDirectory& directory){
    if(!directory.found){
        spdlog::warn("TemplateMatcher: lattices directory '{}' not found; no user templates loaded", directory.path);
    }
    for(const TemplateDirectory::File& file : directory.files){
        if(!file.read){
            spdlog::warn("TemplateMatcher: skipping '{}': {}", file.path, file.error);
        }else if(!file.error.empty()){
            spdlog::warn("TemplateMatcher: skipping template '{}': {}", file.templateName, file.error);
        }else{
            spdlog::info(
                "TemplateMatcher: loaded template '{}' (coordination {}, {} runtime graphs) as structure_type {}",
                file.templateName, file.coordinationNumber, file.graphs, file.structureType
            );
            _templateNames[file.structureType] = file.templateName;
        }
    }
    const std::size_t loaded = directory.loadedTemplates();
    if(loaded == 0){
        spdlog::warn("PTM: no user templates loaded from '{}'", directory.path);
    }else{
        spdlog::info("Loaded {} user templates in {}", loaded, phaseElapsed());
    }
}

void PolyhedralTemplateMatchingLog::dump(const StructureAnalysis&){
    spdlog::info("PTM finished in {}", phaseElapsed());
}

void PolyhedralTemplateMatchingLog::dump(const std::map<int, int>& structureCounts){
    TypeCounts line{"Structure types"};
    for(const auto& [structureType, count] : structureCounts){
        line.atoms += static_cast<std::size_t>(count);
        if(structureType == StructureType::OTHER){
            continue;
        }
        const auto templateName = _templateNames.find(structureType);
        line.counts.emplace_back(
            templateName != _templateNames.end() && !templateName->second.empty()
                ? templateName->second : std::string(structureTypeName(structureType)),
            count
        );
    }
    if(const auto other = structureCounts.find(StructureType::OTHER); other != structureCounts.end()){
        line.counts.emplace_back(structureTypeName(StructureType::OTHER), other->second);
    }
    dump(line);
}

void PolyhedralTemplateMatchingLog::dump(const std::map<int, std::int64_t>& orderingCounts){
    constexpr int none = static_cast<int>(PTM::OrderingType::ORDERING_NONE);
    TypeCounts line{"Ordering types"};
    for(const auto& [orderingType, count] : orderingCounts){
        line.atoms += static_cast<std::size_t>(count);
        if(orderingType != none){
            line.counts.emplace_back(PTM::orderingTypeName(orderingType), count);
        }
    }
    if(const auto unordered = orderingCounts.find(none); unordered != orderingCounts.end()){
        line.counts.emplace_back(PTM::orderingTypeName(none), unordered->second);
    }
    dump(line);
}

void PolyhedralTemplateMatchingLog::dump(const std::vector<PtmLocalAtomState>& atomStates){
    std::size_t rmsdAtoms = 0;
    double rmsdSum = 0.0;
    double rmsdMin = std::numeric_limits<double>::infinity();
    double rmsdMax = 0.0;
    for(const PtmLocalAtomState& state : atomStates){
        if(!state.valid || !(state.rmsd >= 0.0)){
            continue;
        }
        ++rmsdAtoms;
        rmsdSum += state.rmsd;
        rmsdMin = std::min(rmsdMin, state.rmsd);
        rmsdMax = std::max(rmsdMax, state.rmsd);
    }
    if(rmsdAtoms > 0){
        spdlog::info("RMSD of matched atoms: min {:.4f}, mean {:.4f}, max {:.4f}",
                     rmsdMin, rmsdSum / static_cast<double>(rmsdAtoms), rmsdMax);
    }
}

void PolyhedralTemplateMatchingLog::dump(const PtmNeighborOrder& order){
    spdlog::info("  Neighbor ordering done in {}", phaseElapsed());
    spdlog::info("  Structure templates: {} (input lattice {})",
                 structureFamilyNames(order.ptm.ptmCheckFlags()), latticeStructureTypeName(order.inputLattice));
}

void PolyhedralTemplateMatchingLog::dump(const PtmRmsdCutoffReport& matching){
    spdlog::info("  Template matching done for {} atoms in {}", matching.atoms, phaseElapsed());
    if(matching.rejectedAtoms == 0){
        return;
    }
    spdlog::info(
        "  PTM: {} of {} atoms ({:.1f}%) matched a lattice template but were discarded by --rmsd {} "
        "(smallest discarded RMSD {:.4f}, median {:.3f}); raise --rmsd to keep them",
        matching.rejectedAtoms, matching.atoms,
        100.0 * static_cast<double>(matching.rejectedAtoms) / static_cast<double>(matching.atoms),
        matching.rmsdCutoff, matching.smallestRejectedRmsd, matching.medianRejectedRmsd
    );
}

void PolyhedralTemplateMatchingLog::dump(const PtmDiamondShells& shells){
    spdlog::info("  Diamond neighbor shells reordered by CNA in {}: {} atoms kept, {} set to Other",
                 phaseElapsed(), shells.keptAtoms, shells.demotedAtoms);
}

void PolyhedralTemplateMatchingLog::dump(const StructureContext& neighborLists){
    spdlog::info("  Neighbor lists built in {}: {} entries, max neighbor distance {:.3f}",
                 phaseElapsed(), neighborLists.neighborIndices->size(), neighborLists.maximumNeighborDistance);
}

}
