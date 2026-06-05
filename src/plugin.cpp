#include <volt/plugin/plugin_entry.h>
#include <volt/analysis/ptm_service.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = PolyhedralTemplateMatchingService;

static const std::vector<OptionBinding<S>> bindings = {
    optLattice("--crystal_structure", "SC|FCC|HCP|BCC|CUBIC_DIAMOND|HEX_DIAMOND", "FCC", &S::setInputCrystalStructure),
    opt("--rmsd", "RMSD threshold for PTM", 0.1, &S::setRMSD),
    opt("--dissolve_small_clusters", "Mark small clusters as OTHER", false, &S::setDissolveSmallClusters),
};

VOLT_SERVICE_PLUGIN("volt-polyhedral-template-matching", "Polyhedral Template Matching", S, bindings)
