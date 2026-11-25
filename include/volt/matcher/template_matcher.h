#pragma once

#include <volt/math/vector3.h>
#include <volt/math/quaternion.h>

#include <ptm_constants.h>

#include <array>
#include <cstddef>
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

// What loadDirectory() made of a directory of template files: whether it exists, and each file it read, in order,
// with the template it loaded from it or the reason it was skipped.
struct TemplateDirectory{
    struct File{
        std::string path;
        // False when the file could not be read as a template definition.
        bool read = false;
        std::string templateName;
        // Empty when the template was loaded.
        std::string error;
        int coordinationNumber = 0;
        std::size_t graphs = 0;
        int structureType = -1;
    };

    std::string path;
    bool found = false;
    std::vector<File> files;

    std::size_t loadedTemplates() const;
};

class TemplateMatcher{
public:
    TemplateMatcher() = default;

    TemplateDirectory loadDirectory(const std::filesystem::path &directory);

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

    // The template's points and their graph, and the RMSD fit of a neighborhood onto them; the template file reader.
    static double calcRmsd(int numPoints, const double (*idealPoints)[3],
                           double (*normalized)[3], std::int8_t* mapping,
                           double* quaternionOut, double* scaleOut);
    static std::array<double, 3> quatRotate(const double* quaternion, const std::array<double, 3>& vector);
    static std::vector<std::array<double, 3>> buildTemplatePoints(const TemplateDefinition& definition);
    static void enumerateAutomorphisms(int neighborCount, const bool adjacency[PTM_MAX_NBRS][PTM_MAX_NBRS],
                                       const std::int8_t* degree,
                                       std::vector<std::array<std::int8_t, PTM_MAX_POINTS>>& out);
    static bool hullCanonical(int numPoints, const double (*points)[3], int wantFacets,
                              std::int8_t facetsOut[PTM_MAX_FACETS][3], int& numFacetsOut,
                              std::int8_t* degreeOut, std::int8_t* canonicalLabellingOut,
                              std::uint64_t* hashOut, std::int8_t* codeOut, int* numEdgesOut);
    static bool readTemplateFile(const std::filesystem::path& filePath, TemplateDefinition& out, std::string* error);

    std::vector<LoadedTemplate> _templates;

    static constexpr int kGraphDiscoverySamples = 4000;
};

}