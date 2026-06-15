#include <volt/plugin/plugin_entry.h>
#include <volt/analysis/ptm_service.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = PolyhedralTemplateMatchingService;

static const std::vector<OptionBinding<S>> bindings = {
    optLattice("--crystal_structure", "SC|FCC|HCP|BCC|CUBIC_DIAMOND|HEX_DIAMOND", "FCC", &S::setInputCrystalStructure),
    opt("--rmsd", "RMSD threshold for PTM", 0.1, &S::setRMSD),
    opt("--dissolve_small_clusters", "Mark small clusters as OTHER", false, &S::setDissolveSmallClusters),
    opt("--lattices", "Directory of user-defined lattice template YAML files (cell+basis+coordination_number)", "", &S::setLatticesDirectory),
    opt("--cation_radius", "Cation-cation search radius (A) to grow orientation grains over user templates; 0 disables", 0.0, &S::setCationNeighborRadius),
    opt("--cation_misorientation", "Misorientation tolerance (deg) for joining cations into one grain", 12.0, &S::setCationMisorientation),
};

VOLT_SERVICE_PLUGIN("volt-polyhedral-template-matching", "Polyhedral Template Matching", S, bindings)
