#pragma once

#include <cstdint>
#include <volt/math/matrix3.h>
#include <volt/math/quaternion.h>

namespace Volt{

// Per-atom state emitted by PTM. `valid` is true when the structure was
// identified within the RMSD cutoff; the remaining fields are only
// meaningful in that case and should be ignored when `valid` is false.
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

}
