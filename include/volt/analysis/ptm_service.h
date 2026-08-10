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

    // Comma-separated structure families to search: SC,FCC,HCP,ICO,BCC,DCUB,DHEX,
    // GRAPHENE, or ALL. Empty (the default) derives the set from
    // --crystal_structure, which is both faster and avoids spurious diamond
    // matches in metals. Widen this only if the sample really is mixed, or to
    // reach GRAPHENE (no --crystal_structure value maps to it).
    void setStructureCheckFlags(std::string families);
    
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
    // 0 == derive from _inputCrystalStructure.
    int _structureCheckFlags;
};
    
}
