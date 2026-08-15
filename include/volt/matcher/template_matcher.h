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

inline constexpr int TEMPLATE_STRUCTURE_TYPE_BASE = 1000;

struct TemplateDefinition{
    std::string name;
    std::array<Vector3, 3> cell;
    std::vector<Vector3> basisCartesian;
    std::vector<int> basisSpecies;
    int coordinationNumber = 0;
    int referenceBasisAtomIndex = 0;
};

struct TemplateGraph{
    std::uint64_t hash = 0;
    std::vector<std::int8_t> code;
    int numEdges = 0;
    std::array<std::int8_t, PTM_MAX_POINTS> canonicalLabelling{};
    int numFacets = 0;
    std::array<std::array<std::int8_t, 3>, PTM_MAX_FACETS> facets{};
    std::vector<std::array<std::int8_t, PTM_MAX_POINTS>> automorphisms;
};

struct LoadedTemplate{
    std::string name;
    int structureType = -1;
    int numPoints = 0;
    int numNeighbours = 0;
    int referenceNumFacets = -1;
    std::array<std::array<double, 3>, PTM_MAX_POINTS> ideal{};
    std::vector<TemplateGraph> graphs;
};

struct TemplateMatch{
    double rmsd = 0.0;
    double scale = 0.0;
    Quaternion orientation = Quaternion(0.0, 0.0, 0.0, 1.0);
    bool topologicalMatch = false;
    std::array<std::int8_t, PTM_MAX_POINTS> mapping{};
};

class TemplateMatcher{
public:
    TemplateMatcher() = default;

    int loadDirectory(const std::filesystem::path &directory);

    bool addTemplate(
        const TemplateDefinition& definition,
        std::string* error = nullptr
    );

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

    static constexpr int kGraphDiscoverySamples = 4000;
};

}