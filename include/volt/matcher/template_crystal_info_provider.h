#pragma once

#include <volt/analysis/structure_analysis.h>
#include <volt/matcher/template_matcher.h>
#include <volt/math/matrix3.h>
#include <volt/math/vector3.h>

#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace Volt{

class TemplateCrystalInfoProvider final : public StructureAnalysisCrystalInfo{
public:
    TemplateCrystalInfoProvider(
        std::shared_ptr<const StructureAnalysisCrystalInfo> builtin,
        const TemplateMatcher& matcher,
        int cationCoordination
    );

    int findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const override;
    int coordinationNumber(int structureType) const override;
    int commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const override;
    int symmetryPermutationCount(int structureType) const override;
    int symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const override;
    const Matrix3& symmetryTransformation(int structureType, int symmetryIndex) const override;
    int symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const override;
    const Vector3& latticeVector(int structureType, int latticeVectorIndex) const override;
    std::string_view topologyName(int structureType) const override;

private:
    static bool isDefinedType(int structureType){
        return structureType >= TEMPLATE_STRUCTURE_TYPE_BASE;
    }

    struct CationCrystalData{
        int coordinationNumber = 0;
        std::string name;
        std::array<int, MAX_NEIGHBORS> identityPermutation{};
    };

    const CationCrystalData* dataFor(int structureType) const;

    std::shared_ptr<const StructureAnalysisCrystalInfo> _builtin;
    int _cationCoordination;
    std::unordered_map<int, CationCrystalData> _data;
};

}
