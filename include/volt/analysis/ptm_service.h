#pragma once

#include <volt/core/volt.h>
#include <nlohmann/json.hpp>
#include <volt/core/lammps_parser.h>
#include <volt/structures/crystal_structure_types.h>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace Volt{

using json = nlohmann::json;

class TemplateMatcher;
struct OrientationClusterAtomState;
struct PtmLocalAtomState;

class PolyhedralTemplateMatchingService{
public:
    PolyhedralTemplateMatchingService();

    void setInputCrystalStructure(LatticeStructureType structureType);
    void setRMSD(double rmsd);
    void setDissolveSmallClusters(bool dissolveSmallClusters);
    void setLatticesDirectory(std::string latticesDirectory);

    void setCationNeighborRadius(double radius);

    void setCationMisorientation(double degrees);

    void setStructureCheckFlags(std::string families);

    void setIdentifyOrdering(bool identifyOrdering);

    json compute(
        const LammpsParser::Frame& frame,
        const std::string& outputBase,
        const std::string& inputDumpPath
    );

private:
    // The cation radius the template-atom network is grown with: none without user templates.
    double cationRadiusFor(const TemplateMatcher& templates) const;
    // Whether the run identifies chemical ordering, which needs the per-atom types.
    bool identifiesOrdering(const LammpsParser::Frame& frame, std::size_t atomCount) const;

    static Matrix3 quaternionToMatrix(const Quaternion& orientation);
    // The orientation of each atom PTM matched, for the clusters grown by orientation.
    static std::shared_ptr<std::vector<OrientationClusterAtomState>> buildOrientationClusterStates(
        const std::shared_ptr<const std::vector<PtmLocalAtomState>>& atomStates
    );

    LatticeStructureType _inputCrystalStructure;
    double _rmsd;
    bool _dissolveSmallClusters;
    std::string _latticeDirectory;
    double _cationNeighborRadius;
    double _cationMisorientation;
    int _structureCheckFlags;
    bool _identifyOrdering;
};
    
}
