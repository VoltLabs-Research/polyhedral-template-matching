#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/ptm.h>
#include <volt/analysis/internal/ptm_crystal_info_provider.h>

#include <array>
#include <cstddef>

namespace Volt::PtmStructureAnalysisDetail{

inline constexpr std::array<int, 6> kSimpleCubicTemplateToCanonicalNeighborSlot = {5, 4, 3, 2, 1, 0};

bool setupPTM(StructureContext& context, Volt::PTM& ptm, std::size_t particleCount, bool collectDefGradient = false);

}
