#pragma once

#include <volt/analysis/structure_analysis.h>

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>

namespace Volt::PtmStructureAnalysisDetail{

struct PtmSymmetryPermutation{
    std::vector<int> permutation;
    Matrix3 transformation = Matrix3::Identity();
    std::vector<int> inverseProduct;
};

struct PtmCrystalData{
    int coordinationNumber = 0;
    std::vector<Vector3> latticeVectors;
    std::vector<std::array<int, 2>> commonNeighbors;
    std::vector<PtmSymmetryPermutation> symmetries;
    std::vector<int> templateToCanonicalNeighborSlot;
};

class PtmCrystalInfoProvider final : public StructureAnalysisCrystalInfo{
public:
    PtmCrystalInfoProvider();

    std::string_view topologyName(int structureType) const override;
    int findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const override;
    int coordinationNumber(int structureType) const override;
    int commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const override;
    int symmetryPermutationCount(int structureType) const override;
    int symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const override;
    const Matrix3& symmetryTransformation(int structureType, int symmetryIndex) const override;
    int symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const override;
    const Vector3& latticeVector(int structureType, int latticeVectorIndex) const override;
    int templateToCanonicalNeighborSlot(int structureType, int templateSlot) const;

private:
    void initialize(int structureType) const;
    const PtmCrystalData& dataFor(int structureType) const;

    mutable std::unordered_map<int, PtmCrystalData> _data;
};

std::shared_ptr<const StructureAnalysisCrystalInfo> ptmCrystalInfoProvider();
int ptmTemplateToCanonicalNeighborSlot(int structureType, int templateSlot);

}
