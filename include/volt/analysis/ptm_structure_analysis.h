#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/math/matrix3.h>
#include <volt/math/quaternion.h>

#include <cstdint>
#include <memory>
#include <vector>

namespace Volt{

class TemplateMatcher;

struct PtmLocalAtomState{
    Quaternion orientation;
    Matrix3 deformationGradient;
    double rmsd = 0.0;
    double interatomicDistance = 0.0;
    std::uint64_t correspondencesCode = 0;
    int orderingType = 0;
    int bestTemplateIndex = -1;
    bool valid = false;

    PtmLocalAtomState()
        : orientation(Quaternion::Identity{})
        , deformationGradient(Matrix3::Identity()){}
};

void computeMaximumNeighborDistanceFromPTM(StructureAnalysis& analysis);
// structureCheckFlags: PTM_CHECK_* bitmask restricting which structure families
// ptm_index() searches. 0 (the default) derives the set from the analysis
// context's input crystal structure — see PTM::defaultCheckFlagsForLattice(),
// which explains why the derived set is both faster and slightly safer than
// searching everything. Pass an explicit mask to widen (e.g. to include
// PTM_CHECK_GRAPHENE, which no LatticeStructureType maps to).
void determineLocalStructuresWithPTM(
    StructureAnalysis& analysis,
    double rmsdCutoff,
    std::shared_ptr<std::vector<PtmLocalAtomState>> atomStates = nullptr,
    const TemplateMatcher* templates = nullptr,
    double cationNeighborRadius = 0.0,
    int structureCheckFlags = 0
);

}
