#pragma once

#include <volt/analysis/structure_identification_log.h>

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace Volt{

class StructureAnalysis;
class StructureContext;
struct PtmDiamondShells;
struct PtmLocalAtomState;
struct PtmNeighborOrder;
struct PtmRmsdCutoffReport;
struct TemplateDirectory;

// The execution log of one PTM run, on top of the lines every structure identification prints.
class PolyhedralTemplateMatchingLog : public StructureIdentificationLog{
public:
    using StructureIdentificationLog::dump;

    // The user templates it loads name their atoms in the structure counts.
    void dump(const TemplateDirectory& directory);
    void dump(const StructureAnalysis& identified);
    // The atoms of each structure type, and of each chemical ordering type, keyed by type.
    void dump(const std::map<int, int>& structureCounts);
    void dump(const std::map<int, std::int64_t>& orderingCounts);
    void dump(const std::vector<PtmLocalAtomState>& atomStates);

    void dump(const PtmNeighborOrder& order);
    void dump(const PtmRmsdCutoffReport& matching);
    void dump(const PtmDiamondShells& shells);
    void dump(const StructureContext& neighborLists);

private:
    // "SC, FCC, HCP": the structure families PTM checks.
    static std::string structureFamilyNames(int checkFlags);

    std::map<int, std::string> _templateNames;
};

}
