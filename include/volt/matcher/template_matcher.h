#pragma once

#include <volt/math/vector3.h>
#include <volt/math/quaternion.h>

#include <ptm_constants.h>

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace Volt{

// The IDs for dynamic structures start from this value.
inline constexpr int TEMPLATE_STRUCTURE_TYPE_BASE = 1000;

struct TemplateDefinition{
    std::string name;
    // unit-cell column vectors
    std::array<Vector3, 3> cell;
    // basis-atom positions, Cartesian
    std::vector<Vector3> basisCartesian;
    // 1:1 with basisCartesian
    std::vector<int> basisSpecies;
    // neighbours that define the template shell
    int coordinationNumber = 0;
    // which basis atom is the central atom
    int referenceBasisAtomIndex = 0;
};

// Runtime-derived triangulation of the template's ideal convex hull,
// wth its canonical labelling and automorphism permutations. 
// PTM ships per built-in structure in graphs_*[] + automorphisms[]. 
struct TemplateGraph{
    // 64-bit Weinberg hash
    std::uint64_t hash = 0;
    // full Weinberg canonical code
    std::vector<std::int8_t> code;
    int numEdges = 0;
    // point-index -> canonical label
    std::array<std::int8_t, PTM_MAX_POINTS> canonicalLabelling{};
    int numFacets = 0;
    std::array<std::array<std::int8_t, 3>, PTM_MAX_FACETS> facets{};
    // point-index permutations, [0] fixed = 0
    std::vector<std::array<std::int8_t, PTM_MAX_POINTS>> automorphisms;
};

// Ready-to-match template. The normalized ideal point set + its
// runtime graph multiplicity set.
struct LoadedTemplate{
    std::string name;
    // TEMPLATE_STRUCTURE_TYPE_BASE + index
    int structureType = -1;
    // coordinationNumber + 1
    int numPoints = 0;
    int numNeighbours = 0;
    int referenceNumFacets = -1;
    // barycentre 0, mean nbr dist 1
    std::array<std::array<double, 3>, PTM_MAX_POINTS> ideal{};
    std::vector<TemplateGraph> graphs;
};

// Result of matching one measured environment against one loaded template.
struct TemplateMatch{
    double rmsd = 0.0;
    double scale = 0.0;
    Quaternion orientation = Quaternion(0.0, 0.0, 0.0, 1.0);
    bool topologicalMatch = false;
    // ideal[i] <-> measured[mapping[i]]
    std::array<std::int8_t, PTM_MAX_POINTS> mapping{};
};

class TemplateMatcher{
public:
    TemplateMatcher() = default;

    // Load every *.yml / *.yaml as a template.
    int loadDirectory(const std::filesystem::path &directory);

    // Parse + compile a single template def.
    bool addTemplate(
        const TemplateDefinition& definition,
        std::string* error = nullptr
    );

    // Match a measured environment (point[0] = central atom at origin, point[1 .. numNeighbors] = neighbours deltas)
    // against all loaded templates and return the lowest-RMSD topological match.
    TemplateMatch matchBest(
        const double (*points)[3],
        int numEnvPoints,
        int* outStructureType
    ) const;

    std::string structureName(int structureType) const;

    int coordinationNumber(int structureType) const;

    bool empty() const{ return _templates.empty(); }

    const std::vector<LoadedTemplate>& templates() const{ return _templates; }

private:
    bool compile(
        const TemplateDefinition& definition,
        LoadedTemplate& out,
        std::string* error
    ) const;

    std::vector<LoadedTemplate> _templates;

    // Number of random rotation/noise samples used to discover the template's
    // convex-hull triangulation multiplicity set. 
    static constexpr int kGraphDiscoverySamples = 4000;
};

}