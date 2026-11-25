#pragma once

#include <volt/analysis/cluster_rule_provider.h>
#include <volt/math/matrix3.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <vector>

namespace Volt{

struct OrientationBasedAtomState{
    Matrix3 orientation = Matrix3::Identity();
    bool valid = false;
};

class OrientationBasedClusterRule final : public ClusterRuleProvider{
public:
    // The same-template neighbor pairs the reference misorientation angles were fitted on, and how many of them the
    // kept angles cover.
    struct AngleFit{
        long pairs = 0;
        long coveredPairs = 0;
    };

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

    // None when the rule had no neighbor lists to fit the angles on.
    const std::optional<AngleFit>& angleFit() const{ return _angleFit; }
    std::size_t referenceAngleCount() const{ return _referenceAngles.size(); }

private:
    const OrientationBasedAtomState* stateFor(int atomIndex) const;

    // The angle of the rotation from orientation a to orientation b, in radians.
    static double misorientationAngle(const Matrix3& a, const Matrix3& b);

    bool matchesReference(const Matrix3& a, const Matrix3& b) const;

    std::shared_ptr<const std::vector<OrientationBasedAtomState>> _atomStates;

    double _toleranceRad;

    std::vector<double> _referenceAngles;

    std::optional<AngleFit> _angleFit;
};

}
