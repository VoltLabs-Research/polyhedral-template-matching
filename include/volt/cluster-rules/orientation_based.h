#pragma once

#include <volt/analysis/cluster_rule_provider.h>
#include <volt/math/matrix3.h>

#include <memory>
#include <vector>

namespace Volt{

struct OrientationBasedAtomState{
    Matrix3 orientation = Matrix3::Identity();
    bool valid = false;
};

class OrientationBasedClusterRule final : public ClusterRuleProvider{
public:
    OrientationBasedClusterRule(
        std::shared_ptr<const std::vector<OrientationBasedAtomState>> atomStates,
        const int* neighborOffsets,
        const int* neighborIndices,
        const int* structureTypes,
        std::size_t atomCount,
        double toleranceDegrees
    );

    void initializeClusterSeed(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        Cluster& cluster,
        int seedAtomIndex,
        int structureType
    ) const override;

    bool finalizeClusterOrientation(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        Cluster& cluster,
        int seedAtomIndex,
        int structureType
    ) const override;

    ClusterRuleDecision tryAssignNeighbor(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        const Cluster& cluster,
        int currentAtomIndex,
        int neighborAtomIndex,
        int neighborIndex,
        int structureType,
        int& outNeighborSymmetry
    ) const override;

    ClusterRuleDecision tryCalculateTransition(
        const StructureAnalysis& analysis,
        const AnalysisContext& context,
        int atomIndex,
        int neighborAtomIndex,
        int neighborIndex,
        Matrix3& outTransition
    ) const override;

private:
    const OrientationBasedAtomState* stateFor(int atomIndex) const;

    bool matchesReference(const Matrix3& a, const Matrix3& b) const;

    std::shared_ptr<const std::vector<OrientationBasedAtomState>> _atomStates;

    double _toleranceRad;

    std::vector<double> _referenceAngles;
};

}
