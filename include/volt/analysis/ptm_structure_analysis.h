#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/math/matrix3.h>
#include <volt/math/quaternion.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace Volt{

class PTM;
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

// The atoms that matched a lattice template but were discarded for an RMSD above the cutoff, out of the atoms PTM
// matched.
struct PtmRmsdCutoffReport{
    std::size_t rejectedAtoms = 0;
    double medianRejectedRmsd = 0.0;
    double smallestRejectedRmsd = 0.0;
    std::size_t atoms = 0;
    double rmsdCutoff = 0.0;
};

// PTM once it has cached the ordered neighbors of every atom, about to match them against the templates of its
// input lattice.
struct PtmNeighborOrder{
    const PTM& ptm;
    LatticeStructureType inputLattice;
};

// The neighbor shells the CNA ordering gave the diamond atoms, whether each atom has one, and the diamond atoms
// kept and set to Other for lacking one.
struct PtmDiamondShells{
    std::vector<std::array<int, MAX_NEIGHBORS>> neighbors;
    std::vector<unsigned char> valid;
    std::size_t keptAtoms = 0;
    std::size_t demotedAtoms = 0;
};

// The same-template neighbors of each user-template atom within the cation radius.
struct PtmCationNetwork{
    std::vector<std::array<int, MAX_NEIGHBORS>> neighbors;
    double radius = 0.0;
};

void computeMaximumNeighborDistanceFromPTM(StructureAnalysis& analysis);
void determineLocalStructuresWithPTM(
    StructureAnalysis& analysis,
    double rmsdCutoff,
    std::shared_ptr<std::vector<PtmLocalAtomState>> atomStates = nullptr,
    const TemplateMatcher* templates = nullptr,
    double cationNeighborRadius = 0.0,
    int structureCheckFlags = 0,
    const int* particleTypes = nullptr,
    PtmRmsdCutoffReport* rmsdCutoffReport = nullptr
);

namespace PtmStructureAnalysisDetail{

// Sets PTM up for the context's input crystal over its positions; false when it cannot be.
bool setupPTM(StructureContext& context, PTM& ptm, size_t particleCount, bool collectDefGradient,
              int structureCheckFlags, const int* particleTypes);

}

}
