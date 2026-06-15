#include <volt/matcher/template_matcher.h>

#include <ptm_convex_hull_incremental.h>
#include <ptm_graph_tools.h>
#include <ptm_canonical_coloured.h>
#include <ptm_normalize_vertices.h>
#include <ptm_quat.h>
#include <ptm_structure_matcher.h>

#include <yaml-cpp/yaml.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <random>
#include <set>
#include <string>

namespace Volt{

namespace{

using Vec3 = std::array<double, 3>;

double calcRmsd(int numPoints, const double (*idealPoints)[3],
                double (*normalized)[3], std::int8_t* mapping,
                double* quaternionOut, double* scaleOut)
{
    double gramIdeal = 0;
    double gramMeasured = 0;
    for(int pointIndex = 0; pointIndex < numPoints; ++pointIndex){
        gramIdeal += idealPoints[pointIndex][0] * idealPoints[pointIndex][0] +
                     idealPoints[pointIndex][1] * idealPoints[pointIndex][1] +
                     idealPoints[pointIndex][2] * idealPoints[pointIndex][2];
        gramMeasured += normalized[pointIndex][0] * normalized[pointIndex][0] +
                        normalized[pointIndex][1] * normalized[pointIndex][1] +
                        normalized[pointIndex][2] * normalized[pointIndex][2];
    }
    const double e0 = (gramIdeal + gramMeasured) / 2;
    return ptm::calc_rmsd(numPoints, idealPoints, normalized, mapping, gramIdeal, gramMeasured, e0, quaternionOut, scaleOut);
}

Vec3 quatRotate(const double* quaternion, const Vec3& vector){
    double rotation[9];
    ptm::quaternion_to_rotation_matrix(const_cast<double*>(quaternion), rotation);
    return {
        rotation[0] * vector[0] + rotation[1] * vector[1] + rotation[2] * vector[2],
        rotation[3] * vector[0] + rotation[4] * vector[1] + rotation[5] * vector[2],
        rotation[6] * vector[0] + rotation[7] * vector[1] + rotation[8] * vector[2]
    };
}

std::vector<Vec3> buildTemplatePoints(const TemplateDefinition& definition){
    const Vec3 center = {
        definition.basisCartesian[static_cast<std::size_t>(definition.referenceBasisAtomIndex)].x(),
        definition.basisCartesian[static_cast<std::size_t>(definition.referenceBasisAtomIndex)].y(),
        definition.basisCartesian[static_cast<std::size_t>(definition.referenceBasisAtomIndex)].z()
    };

    std::vector<std::pair<double, Vec3>> neighbors;
    const int imageRange = 3;
    for(int imageX = -imageRange; imageX <= imageRange; ++imageX){
        for(int imageY = -imageRange; imageY <= imageRange; ++imageY){
            for(int imageZ = -imageRange; imageZ <= imageRange; ++imageZ){
                const Vector3 translation =
                    definition.cell[0] * static_cast<double>(imageX) +
                    definition.cell[1] * static_cast<double>(imageY) +
                    definition.cell[2] * static_cast<double>(imageZ);
                for(const Vector3& basis : definition.basisCartesian){
                    const Vec3 delta = {
                        basis.x() + translation.x() - center[0],
                        basis.y() + translation.y() - center[1],
                        basis.z() + translation.z() - center[2]
                    };
                    const double distance = std::sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
                    if(distance < 1e-9){
                        continue;
                    }
                    neighbors.emplace_back(distance, delta);
                }
            }
        }
    }
    std::sort(neighbors.begin(), neighbors.end(), [](const auto& left, const auto& right){
        return left.first < right.first;
    });

    std::vector<Vec3> points;
    points.push_back({0, 0, 0});
    for(int neighborIndex = 0; neighborIndex < definition.coordinationNumber && neighborIndex < static_cast<int>(neighbors.size()); ++neighborIndex){
        points.push_back(neighbors[static_cast<std::size_t>(neighborIndex)].second);
    }
    return points;
}

void enumerateAutomorphisms(int neighborCount, const bool adjacency[PTM_MAX_NBRS][PTM_MAX_NBRS],
                            const std::int8_t* degree,
                            std::vector<std::array<std::int8_t, PTM_MAX_POINTS>>& out)
{
    out.clear();
    std::vector<int> permutation(static_cast<std::size_t>(neighborCount), -1);
    std::vector<char> used(static_cast<std::size_t>(neighborCount), 0);

    std::function<void(int)> recurse = [&](int node){
        if(node == neighborCount){
            std::array<std::int8_t, PTM_MAX_POINTS> automorphism;
            automorphism.fill(0);
            automorphism[0] = 0;
            for(int index = 0; index < neighborCount; ++index){
                automorphism[static_cast<std::size_t>(index + 1)] =
                    static_cast<std::int8_t>(permutation[static_cast<std::size_t>(index)] + 1);
            }
            out.push_back(automorphism);
            return;
        }
        for(int candidate = 0; candidate < neighborCount; ++candidate){
            if(used[static_cast<std::size_t>(candidate)]){
                continue;
            }
            if(degree[candidate] != degree[node]){
                continue;
            }
            bool consistent = true;
            for(int previous = 0; previous < node; ++previous){
                if(adjacency[node][previous] != adjacency[candidate][permutation[static_cast<std::size_t>(previous)]]){
                    consistent = false;
                    break;
                }
            }
            if(!consistent){
                continue;
            }
            permutation[static_cast<std::size_t>(node)] = candidate;
            used[static_cast<std::size_t>(candidate)] = 1;
            recurse(node + 1);
            used[static_cast<std::size_t>(candidate)] = 0;
            permutation[static_cast<std::size_t>(node)] = -1;
        }
    };
    recurse(0);
}

bool hullCanonical(int numPoints, const double (*points)[3], int wantFacets,
                   std::int8_t facetsOut[PTM_MAX_FACETS][3], int& numFacetsOut,
                   std::int8_t* degreeOut, std::int8_t* canonicalLabellingOut,
                   std::uint64_t* hashOut, std::int8_t* codeOut, int* numEdgesOut)
{
    const int numNeighbors = numPoints - 1;

    double hullPoints[PTM_MAX_INPUT_POINTS][3];
    double inputPoints[PTM_MAX_INPUT_POINTS][3];
    for(int pointIndex = 0; pointIndex < numPoints; ++pointIndex){
        for(int axis = 0; axis < 3; ++axis){
            inputPoints[pointIndex][axis] = points[pointIndex][axis];
        }
    }
    ptm::normalize_vertices(numPoints, inputPoints, hullPoints);

    ptm::convexhull_t convexHull;
    convexHull.ok = false;
    std::int8_t facets[PTM_MAX_FACETS][3];
    const int returnCode = ptm::get_convex_hull(numPoints, (const double (*)[3])hullPoints, &convexHull, facets);
    if(returnCode != 0){
        return false;
    }
    if(wantFacets > 0 && convexHull.num_facets != wantFacets){
        return false;
    }

    std::int8_t degree[PTM_MAX_NBRS];
    ptm::graph_degree(convexHull.num_facets, facets, numNeighbors, degree);

    std::int8_t code[2 * PTM_MAX_EDGES];
    std::int8_t colours[PTM_MAX_POINTS] = {0};
    std::int8_t canonicalLabelling[PTM_MAX_POINTS];
    std::uint64_t hash = 0;
    if(ptm::canonical_form_coloured(convexHull.num_facets, facets, numNeighbors, degree, colours, canonicalLabelling, code, &hash) != PTM_NO_ERROR){
        return false;
    }

    std::memcpy(facetsOut, facets, sizeof(std::int8_t) * 3 * static_cast<std::size_t>(convexHull.num_facets));
    numFacetsOut = convexHull.num_facets;
    std::memcpy(degreeOut, degree, sizeof(std::int8_t) * static_cast<std::size_t>(numNeighbors));
    std::memcpy(canonicalLabellingOut, canonicalLabelling, sizeof(std::int8_t) * static_cast<std::size_t>(numPoints));
    *hashOut = hash;
    const int numEdges = 3 * convexHull.num_facets / 2;
    if(codeOut){
        std::memcpy(codeOut, code, sizeof(std::int8_t) * 2 * static_cast<std::size_t>(numEdges));
    }
    if(numEdgesOut){
        *numEdgesOut = numEdges;
    }
    return true;
}

} // namespace

bool TemplateMatcher::compile(const TemplateDefinition& definition,
                              LoadedTemplate& out, std::string* error) const
{
    const auto fail = [&](const std::string& message){
        if(error){
            *error = message;
        }
        return false;
    };

    if(definition.coordinationNumber < 1 || definition.coordinationNumber > PTM_MAX_NBRS){
        return fail("coordination_number must be in [1, " + std::to_string(PTM_MAX_NBRS) + "]");
    }
    if(definition.basisCartesian.empty()){
        return fail("template has an empty basis");
    }
    if(definition.referenceBasisAtomIndex < 0 ||
       definition.referenceBasisAtomIndex >= static_cast<int>(definition.basisCartesian.size())){
        return fail("reference basis atom index out of range");
    }

    const std::vector<Vec3> templatePoints = buildTemplatePoints(definition);
    out.numPoints = static_cast<int>(templatePoints.size());
    out.numNeighbours = out.numPoints - 1;
    if(out.numNeighbours != definition.coordinationNumber){
        return fail("could not gather coordination_number neighbours from the cell/basis "
                    "(got " + std::to_string(out.numNeighbours) + ")");
    }
    if(out.numPoints > PTM_MAX_POINTS){
        return fail("coordination_number too large for PTM (max " + std::to_string(PTM_MAX_NBRS) + ")");
    }

    out.name = definition.name;

    double inputPoints[PTM_MAX_INPUT_POINTS][3];
    for(int pointIndex = 0; pointIndex < out.numPoints; ++pointIndex){
        inputPoints[pointIndex][0] = templatePoints[static_cast<std::size_t>(pointIndex)][0];
        inputPoints[pointIndex][1] = templatePoints[static_cast<std::size_t>(pointIndex)][1];
        inputPoints[pointIndex][2] = templatePoints[static_cast<std::size_t>(pointIndex)][2];
    }
    double normalizedIdeal[PTM_MAX_INPUT_POINTS][3];
    ptm::normalize_vertices(out.numPoints, inputPoints, normalizedIdeal);
    for(int pointIndex = 0; pointIndex < out.numPoints; ++pointIndex){
        out.ideal[static_cast<std::size_t>(pointIndex)] =
            {normalizedIdeal[pointIndex][0], normalizedIdeal[pointIndex][1], normalizedIdeal[pointIndex][2]};
    }

    {
        std::int8_t facets[PTM_MAX_FACETS][3];
        int numFacets = 0;
        std::int8_t degree[PTM_MAX_NBRS];
        std::int8_t canonicalLabelling[PTM_MAX_POINTS];
        std::uint64_t hash = 0;
        double idealPoints[PTM_MAX_INPUT_POINTS][3];
        for(int pointIndex = 0; pointIndex < out.numPoints; ++pointIndex){
            idealPoints[pointIndex][0] = out.ideal[static_cast<std::size_t>(pointIndex)][0];
            idealPoints[pointIndex][1] = out.ideal[static_cast<std::size_t>(pointIndex)][1];
            idealPoints[pointIndex][2] = out.ideal[static_cast<std::size_t>(pointIndex)][2];
        }
        if(!hullCanonical(out.numPoints, (const double (*)[3])idealPoints, -1, facets, numFacets, degree, canonicalLabelling, &hash, nullptr, nullptr)){
            return fail("ideal convex hull is degenerate (template shell is not a clean convex "
                        "coordination polyhedron at this coordination_number)");
        }
        out.referenceNumFacets = numFacets;
    }

    std::mt19937 randomEngine(0x7A6B3C1Du ^ static_cast<unsigned>(out.numPoints * 2654435761u));
    std::uniform_real_distribution<double> uniform(-1.0, 1.0);
    std::set<std::string> seenTriangulations;

    for(int sample = 0; sample < kGraphDiscoverySamples; ++sample){
        double quaternion[4];
        double quaternionNorm = 0;
        for(int component = 0; component < 4; ++component){
            quaternion[component] = uniform(randomEngine);
            quaternionNorm += quaternion[component] * quaternion[component];
        }
        quaternionNorm = std::sqrt(quaternionNorm);
        for(int component = 0; component < 4; ++component){
            quaternion[component] /= quaternionNorm;
        }

        const double sigma = (sample % 3 == 0) ? 0.0
            : (0.01 + 0.04 * ((static_cast<unsigned>(sample) * 2654435761u >> 8 & 0xff) / 255.0));
        std::normal_distribution<double> noise(0.0, sigma);

        double points[PTM_MAX_INPUT_POINTS][3];
        for(int pointIndex = 0; pointIndex < out.numPoints; ++pointIndex){
            const Vec3 rotated = quatRotate(quaternion, out.ideal[static_cast<std::size_t>(pointIndex)]);
            points[pointIndex][0] = rotated[0] + noise(randomEngine);
            points[pointIndex][1] = rotated[1] + noise(randomEngine);
            points[pointIndex][2] = rotated[2] + noise(randomEngine);
        }

        std::int8_t facets[PTM_MAX_FACETS][3];
        int numFacets = 0;
        std::int8_t degree[PTM_MAX_NBRS];
        std::int8_t canonicalLabelling[PTM_MAX_POINTS];
        std::uint64_t hash = 0;
        std::int8_t code[2 * PTM_MAX_EDGES];
        int numEdges = 0;
        if(!hullCanonical(out.numPoints, (const double (*)[3])points, out.referenceNumFacets,
                          facets, numFacets, degree, canonicalLabelling, &hash, code, &numEdges)){
            continue;
        }

        std::vector<std::array<int, 3>> sortedFacets;
        sortedFacets.reserve(static_cast<std::size_t>(numFacets));
        for(int facetIndex = 0; facetIndex < numFacets; ++facetIndex){
            std::array<int, 3> triangle = {facets[facetIndex][0], facets[facetIndex][1], facets[facetIndex][2]};
            std::sort(triangle.begin(), triangle.end());
            sortedFacets.push_back(triangle);
        }
        std::sort(sortedFacets.begin(), sortedFacets.end());
        std::string key;
        key.reserve(static_cast<std::size_t>(3 * numFacets));
        for(const auto& triangle : sortedFacets){
            key += static_cast<char>(triangle[0]);
            key += static_cast<char>(triangle[1]);
            key += static_cast<char>(triangle[2]);
        }
        if(!seenTriangulations.insert(key).second){
            continue;
        }

        TemplateGraph graph;
        graph.hash = hash;
        graph.numEdges = numEdges;
        graph.code.assign(code, code + 2 * numEdges);
        graph.numFacets = numFacets;
        for(int facetIndex = 0; facetIndex < numFacets; ++facetIndex){
            graph.facets[static_cast<std::size_t>(facetIndex)] =
                {facets[facetIndex][0], facets[facetIndex][1], facets[facetIndex][2]};
        }
        std::memcpy(graph.canonicalLabelling.data(), canonicalLabelling, sizeof(std::int8_t) * static_cast<std::size_t>(out.numPoints));

        bool adjacency[PTM_MAX_NBRS][PTM_MAX_NBRS];
        std::memset(adjacency, 0, sizeof(adjacency));
        for(int facetIndex = 0; facetIndex < numFacets; ++facetIndex){
            const int vertexA = facets[facetIndex][0];
            const int vertexB = facets[facetIndex][1];
            const int vertexC = facets[facetIndex][2];
            adjacency[vertexA][vertexB] = adjacency[vertexB][vertexA] = true;
            adjacency[vertexB][vertexC] = adjacency[vertexC][vertexB] = true;
            adjacency[vertexC][vertexA] = adjacency[vertexA][vertexC] = true;
        }
        enumerateAutomorphisms(out.numNeighbours, adjacency, degree, graph.automorphisms);
        out.graphs.push_back(std::move(graph));
    }

    if(out.graphs.empty()){
        return fail("no convex-hull triangulation could be generated for the template");
    }
    return true;
}

bool TemplateMatcher::addTemplate(const TemplateDefinition& definition, std::string* error){
    LoadedTemplate loaded;
    if(!compile(definition, loaded, error)){
        return false;
    }
    loaded.structureType = TEMPLATE_STRUCTURE_TYPE_BASE + static_cast<int>(_templates.size());
    _templates.push_back(std::move(loaded));
    return true;
}

namespace{

bool readTemplateFile(const std::filesystem::path& filePath, TemplateDefinition& out, std::string* error){
    const auto fail = [&](const std::string& message){
        if(error){
            *error = message + " ('" + filePath.string() + "')";
        }
        return false;
    };

    YAML::Node document;
    try{
        document = YAML::LoadFile(filePath.string());
    }catch(const std::exception& e){
        return fail(std::string("unable to parse YAML: ") + e.what());
    }
    if(!document || !document.IsMap()){
        return fail("template YAML must contain a mapping root");
    }

    try{
        out.name = document["name"] ? document["name"].as<std::string>() : filePath.stem().string();

        if(!document["coordination_number"]){
            return fail("missing coordination_number");
        }
        out.coordinationNumber = document["coordination_number"].as<int>();

        const double scale = document["scale"] ? document["scale"].as<double>() : 1.0;

        const YAML::Node cellNode = document["cell"];
        if(!cellNode || !cellNode.IsSequence() || cellNode.size() != 3){
            return fail("cell must contain exactly three vectors");
        }
        const auto parseVector = [](const YAML::Node& node){
            if(!node || !node.IsSequence() || node.size() != 3){
                throw std::runtime_error("expected a 3-component vector");
            }
            return Vector3(node[0].as<double>(), node[1].as<double>(), node[2].as<double>());
        };
        for(int axis = 0; axis < 3; ++axis){
            out.cell[static_cast<std::size_t>(axis)] = parseVector(cellNode[axis]) * scale;
        }

        std::string coordinateMode = document["coordinate_mode"] ? document["coordinate_mode"].as<std::string>() : "fractional";
        std::transform(coordinateMode.begin(), coordinateMode.end(), coordinateMode.begin(),
                       [](unsigned char character){ return static_cast<char>(std::tolower(character)); });
        const bool cartesian = (coordinateMode == "cartesian" || coordinateMode == "cart" || coordinateMode == "c");
        if(!cartesian && coordinateMode != "fractional" && coordinateMode != "reduced" && coordinateMode != "direct"){
            return fail("coordinate_mode must be 'fractional' or 'cartesian'");
        }

        const YAML::Node basisNode = document["basis"];
        if(!basisNode || !basisNode.IsSequence() || basisNode.size() == 0){
            return fail("basis must contain at least one atom");
        }
        for(const auto& atom : basisNode){
            if(!atom || !atom.IsMap()){
                return fail("basis entries must be mappings");
            }
            const int species = atom["species"] ? atom["species"].as<int>() : 1;
            Vector3 coordinate = parseVector(atom["position"]);
            if(cartesian){
                coordinate *= scale;
            }else{
                coordinate = out.cell[0] * coordinate.x() + out.cell[1] * coordinate.y() + out.cell[2] * coordinate.z();
            }
            out.basisCartesian.push_back(coordinate);
            out.basisSpecies.push_back(species);
        }

        if(document["reference_basis_atom_index"]){
            out.referenceBasisAtomIndex = document["reference_basis_atom_index"].as<int>();
        }
    }catch(const std::exception& e){
        return fail(std::string("invalid template: ") + e.what());
    }
    return true;
}

} // namespace

int TemplateMatcher::loadDirectory(const std::filesystem::path& directory){
    std::error_code errorCode;
    if(!std::filesystem::exists(directory, errorCode) || !std::filesystem::is_directory(directory, errorCode)){
        spdlog::warn("TemplateMatcher: lattices directory '{}' not found; no user templates loaded",
                     directory.string());
        return 0;
    }

    std::vector<std::filesystem::path> files;
    for(const auto& entry : std::filesystem::directory_iterator(directory, errorCode)){
        if(!entry.is_regular_file()){
            continue;
        }
        std::string extension = entry.path().extension().string();
        std::transform(extension.begin(), extension.end(), extension.begin(),
                       [](unsigned char character){ return static_cast<char>(std::tolower(character)); });
        if(extension == ".yml" || extension == ".yaml"){
            files.push_back(entry.path());
        }
    }
    std::sort(files.begin(), files.end());

    int loadedCount = 0;
    for(const auto& file : files){
        TemplateDefinition definition;
        std::string error;
        if(!readTemplateFile(file, definition, &error)){
            spdlog::warn("TemplateMatcher: skipping '{}': {}", file.string(), error);
            continue;
        }
        if(!addTemplate(definition, &error)){
            spdlog::warn("TemplateMatcher: skipping template '{}': {}", definition.name, error);
            continue;
        }
        spdlog::info("TemplateMatcher: loaded template '{}' (coordination {}, {} runtime graphs) as structure_type {}",
                     definition.name, definition.coordinationNumber, _templates.back().graphs.size(),
                     _templates.back().structureType);
        ++loadedCount;
    }
    return loadedCount;
}

TemplateMatch TemplateMatcher::matchBest(const double (*points)[3], int numEnvPoints,
                                         int* outStructureType) const
{
    TemplateMatch best;
    best.rmsd = std::numeric_limits<double>::infinity();
    int bestType = -1;

    for(const LoadedTemplate& loaded : _templates){
        if(numEnvPoints < loaded.numPoints){
            continue;
        }

        const int numPoints = loaded.numPoints;

        std::int8_t facets[PTM_MAX_FACETS][3];
        int numFacets = 0;
        std::int8_t degree[PTM_MAX_NBRS];
        std::int8_t environmentLabelling[PTM_MAX_POINTS];
        std::uint64_t hash = 0;
        std::int8_t code[2 * PTM_MAX_EDGES];
        int numEdges = 0;

        double environmentCopy[PTM_MAX_INPUT_POINTS][3];
        for(int pointIndex = 0; pointIndex < numPoints; ++pointIndex){
            for(int axis = 0; axis < 3; ++axis){
                environmentCopy[pointIndex][axis] = points[pointIndex][axis];
            }
        }
        if(!hullCanonical(numPoints, (const double (*)[3])environmentCopy, loaded.referenceNumFacets,
                          facets, numFacets, degree, environmentLabelling, &hash, code, &numEdges)){
            continue;
        }

        std::int8_t environmentLabellingInverse[PTM_MAX_POINTS];
        for(int pointIndex = 0; pointIndex < numPoints; ++pointIndex){
            environmentLabellingInverse[environmentLabelling[pointIndex]] = static_cast<std::int8_t>(pointIndex);
        }

        double normalized[PTM_MAX_POINTS][3];
        {
            double temporary[PTM_MAX_INPUT_POINTS][3];
            for(int pointIndex = 0; pointIndex < numPoints; ++pointIndex){
                for(int axis = 0; axis < 3; ++axis){
                    temporary[pointIndex][axis] = points[pointIndex][axis];
                }
            }
            ptm::subtract_barycentre(numPoints, temporary, normalized);
        }

        double idealPoints[PTM_MAX_INPUT_POINTS][3];
        for(int pointIndex = 0; pointIndex < numPoints; ++pointIndex){
            idealPoints[pointIndex][0] = loaded.ideal[static_cast<std::size_t>(pointIndex)][0];
            idealPoints[pointIndex][1] = loaded.ideal[static_cast<std::size_t>(pointIndex)][1];
            idealPoints[pointIndex][2] = loaded.ideal[static_cast<std::size_t>(pointIndex)][2];
        }

        std::int8_t mapping[PTM_MAX_POINTS];
        for(const TemplateGraph& graph : loaded.graphs){
            if(graph.hash != hash || graph.numEdges != numEdges){
                continue;
            }
            if(std::memcmp(graph.code.data(), code, sizeof(std::int8_t) * 2 * static_cast<std::size_t>(numEdges)) != 0){
                continue;
            }
            for(const auto& automorphism : graph.automorphisms){
                for(int pointIndex = 0; pointIndex < numPoints; ++pointIndex){
                    mapping[automorphism[static_cast<std::size_t>(pointIndex)]] =
                        environmentLabellingInverse[graph.canonicalLabelling[static_cast<std::size_t>(pointIndex)]];
                }
                double quaternion[4];
                double scale = 0;
                const double rmsd = calcRmsd(numPoints, (const double (*)[3])idealPoints, normalized, mapping, quaternion, &scale);
                if(rmsd < best.rmsd){
                    best.rmsd = rmsd;
                    best.scale = scale;
                    best.topologicalMatch = true;
                    best.orientation = Quaternion(quaternion[1], quaternion[2], quaternion[3], quaternion[0]);
                    std::memcpy(best.mapping.data(), mapping, sizeof(std::int8_t) * static_cast<std::size_t>(numPoints));
                    bestType = loaded.structureType;
                }
            }
        }
    }

    if(outStructureType){
        *outStructureType = bestType;
    }
    if(bestType < 0){
        best.rmsd = 0.0;
        best.topologicalMatch = false;
    }
    return best;
}

std::string TemplateMatcher::structureName(int structureType) const{
    const int index = structureType - TEMPLATE_STRUCTURE_TYPE_BASE;
    if(index < 0 || index >= static_cast<int>(_templates.size())){
        return {};
    }
    return _templates[static_cast<std::size_t>(index)].name;
}

int TemplateMatcher::coordinationNumber(int structureType) const{
    const int index = structureType - TEMPLATE_STRUCTURE_TYPE_BASE;
    if(index < 0 || index >= static_cast<int>(_templates.size())){
        return 0;
    }
    return _templates[static_cast<std::size_t>(index)].numNeighbours;
}

}
