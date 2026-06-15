#pragma once

#include <volt/core/volt.h>
#include <nlohmann/json.hpp>
#include <volt/core/lammps_parser.h>
#include <volt/structures/crystal_structure_types.h>
#include <string>

namespace Volt{

using json = nlohmann::json;

class PolyhedralTemplateMatchingService{
public:
    PolyhedralTemplateMatchingService();

    void setInputCrystalStructure(LatticeStructureType structureType);
    void setRMSD(double rmsd);
    void setDissolveSmallClusters(bool dissolveSmallClusters);
    void setLatticesDirectory(std::string latticesDirectory);

    // Search radius (Angstrom) for the cation-cation subnetwork used to grow
    // orientation-coherent grains over templates. 0 disables grain
    // clustering of templates (per-atom RMSD/orientation only).
    void setCationNeighborRadius(double radius);

    // Misorientation tolerance for joining two cations into one grain.
    void setCationMisorientation(double degrees);
    
    json compute(
        const LammpsParser::Frame& frame,
        const std::string& outputBase,
        const std::string& inputDumpPath
    );

private:
    LatticeStructureType _inputCrystalStructure;
    double _rmsd;
    bool _dissolveSmallClusters;
    std::string _latticeDirectory;
    double _cationNeighborRadius;
    double _cationMisorientation;
};
    
}
