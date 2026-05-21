#include <volt/cli/common.h>
#include <volt/analysis/ptm_service.h>
#include <volt/structures/crystal_structure_types.h>

using namespace Volt;
using namespace Volt::CLI;

void showUsage(const std::string& name){
    printUsageHeader(name, "Volt - Polyhedral Template Matching");
    std::cerr
        << "  --crystal_structure <type>     Crystal structure. (SC|FCC|HCP|BCC|CUBIC_DIAMOND|HEX_DIAMOND) [default: FCC]\n"
        << "  --rmsd <float>                RMSD threshold for PTM. [default: 0.1]\n"
        << "  --dissolve_small_clusters       Mark small clusters as OTHER after building clusters.\n";
    printHelpOption();
}

int main(int argc, char* argv[]){
    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);
    if(const int startupStatus = handleHelpOrMissingInput(argc, argv, opts, filename, showUsage);
       startupStatus >= 0){
        return startupStatus;
    }

    initLogging("volt-polyhedral-template-matching");

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    PolyhedralTemplateMatchingService analyzer;
    LatticeStructureType crystal_structure = LATTICE_FCC;
    const std::string crystal_structureOption = getString(opts, "--crystal_structure", "FCC");
    if(!parseLatticeStructureType(crystal_structureOption, crystal_structure)){
        spdlog::warn("Unknown crystal structure '{}', defaulting to FCC.", crystal_structureOption);
        crystal_structure = LATTICE_FCC;
    }
    analyzer.setInputCrystalStructure(crystal_structure);
    analyzer.setRMSD(getDouble(opts, "--rmsd", 0.1));
    analyzer.setDissolveSmallClusters(hasOption(opts, "--dissolve_small_clusters"));

    spdlog::info("Starting PTM analysis...");
    json result = analyzer.compute(frame, outputBase, filename);
    if(result.value("is_failed", false)){
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }

    spdlog::info("PTM analysis completed.");
    return 0;
}
