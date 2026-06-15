#include <volt/cluster-rules/orientation_based.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/matcher/template_matcher.h>
#include <volt/structures/cluster.h>

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <utility>

namespace Volt{

namespace{

double traceOf(const Matrix3& matrix){
    return matrix(0, 0) + matrix(1, 1) + matrix(2, 2);
}

double misorientationAngle(const Matrix3& a, const Matrix3& b){
    const Matrix3 relative = Matrix3(a.transposed() * b);
    double cosine = (traceOf(relative) - 1.0) / 2.0;
    if(cosine > 1.0){
        cosine = 1.0;
    }
    if(cosine < -1.0){
        cosine = -1.0;
    }
    return std::acos(cosine);
}

} // namespace

OrientationBasedClusterRule::OrientationBasedClusterRule(
    std::shared_ptr<const std::vector<OrientationBasedAtomState>> atomStates,
    const int* neighborOffsets,
    const int* neighborIndices,
    const int* structureTypes,
    std::size_t atomCount,
    double toleranceDegrees
) : _atomStates(std::move(atomStates)){
    _toleranceRad = toleranceDegrees * PI / 180.0;

    if(!neighborOffsets || !neighborIndices || !structureTypes || !_atomStates){
        _referenceAngles.push_back(0.0);
        return;
    }

    constexpr int kBins = 180;
    std::vector<long> histogram(kBins, 0);
    long pairCount = 0;

    for(std::size_t atomIndex = 0; atomIndex < atomCount; ++atomIndex){
        const int type = structureTypes[atomIndex];
        if(type < TEMPLATE_STRUCTURE_TYPE_BASE){
            continue;
        }

        const auto* current = stateFor(static_cast<int>(atomIndex));
        if(!current){
            continue;
        }

        for(int slot = neighborOffsets[atomIndex]; slot < neighborOffsets[atomIndex + 1]; ++slot){
            const int neighbor = neighborIndices[slot];
            if(neighbor < 0 || static_cast<std::size_t>(neighbor) <= atomIndex || structureTypes[neighbor] != type){
                continue;
            }

            const auto* neighborState = stateFor(neighbor);
            if(!neighborState){
                continue;
            }

            const double degrees = misorientationAngle(current->orientation, neighborState->orientation) * 180.0 / PI;
            int bin = static_cast<int>(degrees);
            if(bin < 0){
                bin = 0;
            }
            if(bin >= kBins){
                bin = kBins - 1;
            }
            ++histogram[static_cast<std::size_t>(bin)];
            ++pairCount;
        }
    }

    long covered = 0;
    if(pairCount > 0){
        const long minSupport = std::max<long>(1, pairCount / 50);
        for(int bin = 0; bin < kBins; ++bin){
            if(histogram[static_cast<std::size_t>(bin)] >= minSupport){
                _referenceAngles.push_back((bin + 0.5) * PI / 180.0);
                covered += histogram[static_cast<std::size_t>(bin)];
            }
        }
    }

    _referenceAngles.push_back(0.0);

    spdlog::info("PTM: {} cation pairs, {} reference angles kept, covering {}/{} pairs ({:.1f}%)",
                 pairCount, _referenceAngles.size(), covered, std::max<long>(pairCount, 1),
                 pairCount > 0 ? 100.0 * covered / pairCount : 0.0);
}

const OrientationBasedAtomState* OrientationBasedClusterRule::stateFor(int atomIndex) const{
    if(!_atomStates || atomIndex < 0 || atomIndex >= static_cast<int>(_atomStates->size())){
        return nullptr;
    }
    const auto& state = (*_atomStates)[static_cast<std::size_t>(atomIndex)];
    return state.valid ? &state : nullptr;
}

bool OrientationBasedClusterRule::matchesReference(const Matrix3& a, const Matrix3& b) const{
    const double angle = misorientationAngle(a, b);
    for(const double reference : _referenceAngles){
        if(std::fabs(angle - reference) <= _toleranceRad){
            return true;
        }
    }
    return false;
}

void OrientationBasedClusterRule::initializeClusterSeed(
    const StructureAnalysis&,
    const AnalysisContext&,
    Cluster& cluster,
    int,
    int
) const{
    cluster.symmetryTransformation = 0;
}

bool OrientationBasedClusterRule::finalizeClusterOrientation(
    const StructureAnalysis&,
    const AnalysisContext&,
    Cluster& cluster,
    int seedAtomIndex,
    int
) const{
    const auto* seed = stateFor(seedAtomIndex);
    cluster.orientation = seed ? seed->orientation : Matrix3::Identity();
    return true;
}

ClusterRuleDecision OrientationBasedClusterRule::tryAssignNeighbor(
    const StructureAnalysis&,
    const AnalysisContext& context,
    const Cluster&,
    int currentAtomIndex,
    int neighborAtomIndex,
    int,
    int structureType,
    int& outNeighborSymmetry
) const{
    if(structureType < TEMPLATE_STRUCTURE_TYPE_BASE){
        return ClusterRuleDecision::Unhandled;
    }

    if(context.structureTypes->getInt(neighborAtomIndex) != structureType){
        return ClusterRuleDecision::Rejected;
    }

    const auto* current = stateFor(currentAtomIndex);
    const auto* neighbor = stateFor(neighborAtomIndex);
    if(!current || !neighbor){
        return ClusterRuleDecision::Rejected;
    }

    if(!matchesReference(current->orientation, neighbor->orientation)){
        return ClusterRuleDecision::Rejected;
    }

    outNeighborSymmetry = 0;
    return ClusterRuleDecision::Accepted;
}

ClusterRuleDecision OrientationBasedClusterRule::tryCalculateTransition(
    const StructureAnalysis&,
    const AnalysisContext& context,
    int atomIndex,
    int neighborAtomIndex,
    int,
    Matrix3& outTransition
) const{
    if(context.structureTypes->getInt(atomIndex) < TEMPLATE_STRUCTURE_TYPE_BASE){
        return ClusterRuleDecision::Unhandled;
    }

    const auto* current = stateFor(atomIndex);
    const auto* neighbor = stateFor(neighborAtomIndex);
    if(!current || !neighbor){
        return ClusterRuleDecision::Rejected;
    }

    outTransition = Matrix3(current->orientation.transposed() * neighbor->orientation);
    return ClusterRuleDecision::Accepted;
}

}
