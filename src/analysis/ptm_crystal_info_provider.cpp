#include <volt/analysis/ptm_crystal_info_provider.h>

#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/analysis/ptm.h>
#include <volt/topology/crystal_coordination_topology.h>
#include <volt/topology/crystal_coordination_topology_init.h>

#include <ptm_constants.h>
#include <ptm_initialize_data.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Volt::PtmStructureAnalysisDetail{
namespace{

inline constexpr std::array<int, 6> kSimpleCubicCanonicalToTemplateSlot = {5, 4, 3, 2, 1, 0};
inline constexpr std::array<int, 6> kSimpleCubicTemplateToCanonicalNeighborSlot = {5, 4, 3, 2, 1, 0};

const Vector3& zeroVector(){
    static const Vector3 vector = Vector3::Zero();
    return vector;
}

int normalizedStructureType(int structureType){
    switch(static_cast<StructureType>(structureType)){
        case StructureType::CUBIC_DIAMOND_FIRST_NEIGH:
        case StructureType::CUBIC_DIAMOND_SECOND_NEIGH:
            return static_cast<int>(StructureType::CUBIC_DIAMOND);
        case StructureType::HEX_DIAMOND_FIRST_NEIGH:
        case StructureType::HEX_DIAMOND_SECOND_NEIGH:
            return static_cast<int>(StructureType::HEX_DIAMOND);
        default:
            return structureType;
    }
}

std::vector<int> identityMapping(int count){
    std::vector<int> mapping(static_cast<std::size_t>(count), 0);
    for(int index = 0; index < count; ++index){
        mapping[static_cast<std::size_t>(index)] = index;
    }
    return mapping;
}

PtmCrystalData buildCanonicalDiamondCrystalData(int structureType){
    ensureCoordinationStructuresInitialized();

    PtmCrystalData data;
    const CoordinationStructure& coord = CoordinationStructures::getCoordStruct(structureType);
    const LatticeStructure& lattice = CoordinationStructures::getLatticeStruct(structureType);
    data.coordinationNumber = coord.numNeighbors;
    data.latticeVectors = lattice.latticeVectors;
    data.commonNeighbors.resize(static_cast<std::size_t>(coord.numNeighbors), std::array<int, 2>{-1, -1});
    for(int neighborIndex = 0; neighborIndex < coord.numNeighbors; ++neighborIndex){
        data.commonNeighbors[static_cast<std::size_t>(neighborIndex)] = {
            coord.commonNeighbors[neighborIndex][0],
            coord.commonNeighbors[neighborIndex][1]
        };
    }
    data.symmetries.reserve(lattice.permutations.size());
    for(const auto& symmetry : lattice.permutations){
        PtmSymmetryPermutation ptmSymmetry;
        ptmSymmetry.transformation = symmetry.transformation;
        ptmSymmetry.inverseProduct = symmetry.inverseProduct;
        ptmSymmetry.permutation.assign(
            symmetry.permutation.begin(),
            symmetry.permutation.begin() + coord.numNeighbors
        );
        data.symmetries.push_back(std::move(ptmSymmetry));
    }
    data.templateToCanonicalNeighborSlot = identityMapping(coord.numNeighbors);

    return data;
}

const ptm::refdata_t* ptmReferenceForStructureType(int structureType){
    const StructureType normalized = static_cast<StructureType>(normalizedStructureType(structureType));
    if(normalized == StructureType::OTHER){
        return nullptr;
    }

    const int ptmType = PTM::toPtmStructureType(normalized);
    if(ptmType == PTM_MATCH_NONE){
        return nullptr;
    }

    return ptm::refdata[ptmType];
}

int canonicalTemplateIndex(const ptm::refdata_t& ref){
    if(ref.num_conventional_mappings > 0 && ref.template_indices != nullptr){
        return ref.template_indices[0];
    }
    return 0;
}

const int8_t (*symmetryMappings(const ptm::refdata_t& ref))[PTM_MAX_POINTS]{
    if(ref.num_conventional_mappings > 0 && ref.mapping_conventional != nullptr){
        return ref.mapping_conventional;
    }
    return ref.mapping;
}

int symmetryMappingCount(const ptm::refdata_t& ref){
    if(ref.num_conventional_mappings > 0 && ref.mapping_conventional != nullptr){
        return ref.num_conventional_mappings;
    }
    return ref.num_mappings;
}

std::vector<Vector3> buildCanonicalLatticeVectors(int structureType, const ptm::refdata_t& ref){
    const int numNeighbors = ref.num_nbrs;
    std::vector<Vector3> latticeVectors(static_cast<std::size_t>(numNeighbors), Vector3::Zero());
    const double (*points)[3] = ref.points[canonicalTemplateIndex(ref)];

    if(static_cast<StructureType>(normalizedStructureType(structureType)) == StructureType::SC){
        for(int templateSlot = 0; templateSlot < numNeighbors; ++templateSlot){
            const int canonicalSlot = kSimpleCubicTemplateToCanonicalNeighborSlot[static_cast<std::size_t>(templateSlot)];
            const double* point = points[templateSlot + 1];
            latticeVectors[static_cast<std::size_t>(canonicalSlot)] = Vector3(point[0], point[1], point[2]);
        }
        return latticeVectors;
    }

    for(int neighborSlot = 0; neighborSlot < numNeighbors; ++neighborSlot){
        const double* point = points[neighborSlot + 1];
        latticeVectors[static_cast<std::size_t>(neighborSlot)] = Vector3(point[0], point[1], point[2]);
    }

    if(static_cast<StructureType>(normalizedStructureType(structureType)) == StructureType::BCC){
        double maxAbsComponent = 0.0;
        for(const Vector3& vector : latticeVectors){
            maxAbsComponent = std::max(maxAbsComponent, std::abs(vector.x()));
            maxAbsComponent = std::max(maxAbsComponent, std::abs(vector.y()));
            maxAbsComponent = std::max(maxAbsComponent, std::abs(vector.z()));
        }

        if(maxAbsComponent > EPSILON){
            for(Vector3& vector : latticeVectors){
                vector /= maxAbsComponent;
            }
        }
    }

    return latticeVectors;
}

std::vector<int> buildTemplateToCanonicalMapping(int structureType, int count){
    if(static_cast<StructureType>(normalizedStructureType(structureType)) == StructureType::SC){
        return std::vector<int>(
            kSimpleCubicTemplateToCanonicalNeighborSlot.begin(),
            kSimpleCubicTemplateToCanonicalNeighborSlot.end()
        );
    }
    return identityMapping(count);
}

std::vector<int> buildSymmetryPermutation(int structureType, const ptm::refdata_t& ref, int mappingIndex){
    const int numNeighbors = ref.num_nbrs;
    const int8_t* mapping = symmetryMappings(ref)[mappingIndex];
    std::vector<int> permutation(static_cast<std::size_t>(numNeighbors), 0);

    if(static_cast<StructureType>(normalizedStructureType(structureType)) == StructureType::SC){
        for(int canonicalSlot = 0; canonicalSlot < numNeighbors; ++canonicalSlot){
            const int templateSlot = kSimpleCubicCanonicalToTemplateSlot[static_cast<std::size_t>(canonicalSlot)];
            const int mappedTemplateSlot = mapping[templateSlot + 1] - 1;
            permutation[static_cast<std::size_t>(canonicalSlot)] =
                kSimpleCubicTemplateToCanonicalNeighborSlot[static_cast<std::size_t>(mappedTemplateSlot)];
        }
        return permutation;
    }

    for(int neighborSlot = 0; neighborSlot < numNeighbors; ++neighborSlot){
        permutation[static_cast<std::size_t>(neighborSlot)] = mapping[neighborSlot + 1] - 1;
    }
    return permutation;
}

std::array<int, 3> findNonCoplanarIndices(const std::vector<Vector3>& latticeVectors){
    std::array<int, 3> indices{-1, -1, -1};
    Matrix3 basis = Matrix3::Zero();
    int found = 0;

    for(int vectorIndex = 0; vectorIndex < static_cast<int>(latticeVectors.size()) && found < 3; ++vectorIndex){
        basis.column(found) = latticeVectors[static_cast<std::size_t>(vectorIndex)];

        if(found == 1){
            if(basis.column(0).cross(basis.column(1)).squaredLength() <= EPSILON){
                continue;
            }
        }else if(found == 2){
            if(std::abs(basis.determinant()) <= EPSILON){
                continue;
            }
        }

        indices[static_cast<std::size_t>(found++)] = vectorIndex;
    }

    if(found != 3){
        throw std::runtime_error("Unable to determine a non-coplanar PTM basis.");
    }

    return indices;
}

std::vector<std::array<int, 2>> buildCommonNeighbors(
    int structureType,
    const std::vector<Vector3>& latticeVectors
){
    const int numNeighbors = static_cast<int>(latticeVectors.size());
    std::vector<std::array<int, 2>> commonNeighbors(
        static_cast<std::size_t>(numNeighbors),
        std::array<int, 2>{-1, -1}
    );

    for(int neighborIndex = 0; neighborIndex < numNeighbors; ++neighborIndex){
        Matrix3 basis = Matrix3::Zero();
        basis.column(0) = latticeVectors[static_cast<std::size_t>(neighborIndex)];

        if(static_cast<StructureType>(normalizedStructureType(structureType)) == StructureType::SC){
            for(int i1 = 0; i1 < numNeighbors; ++i1){
                if(i1 == neighborIndex){
                    continue;
                }
                basis.column(1) = latticeVectors[static_cast<std::size_t>(i1)];
                for(int i2 = i1 + 1; i2 < numNeighbors; ++i2){
                    if(i2 == neighborIndex){
                        continue;
                    }
                    basis.column(2) = latticeVectors[static_cast<std::size_t>(i2)];
                    if(std::abs(basis.determinant()) > EPSILON){
                        commonNeighbors[static_cast<std::size_t>(neighborIndex)] = {i1, i2};
                        goto next_neighbor;
                    }
                }
            }
        }

        {
            double minBondDistanceSq = std::numeric_limits<double>::max();
            for(int otherIndex = 0; otherIndex < numNeighbors; ++otherIndex){
                if(otherIndex == neighborIndex){
                    continue;
                }
                const double distSq = (
                    latticeVectors[static_cast<std::size_t>(neighborIndex)] -
                    latticeVectors[static_cast<std::size_t>(otherIndex)]
                ).squaredLength();
                if(distSq > EPSILON && distSq < minBondDistanceSq){
                    minBondDistanceSq = distSq;
                }
            }

            if(!std::isfinite(minBondDistanceSq)){
                goto next_neighbor;
            }

            for(int i1 = 0; i1 < numNeighbors; ++i1){
                if(i1 == neighborIndex){
                    continue;
                }
                const double d1 = (
                    latticeVectors[static_cast<std::size_t>(neighborIndex)] -
                    latticeVectors[static_cast<std::size_t>(i1)]
                ).squaredLength();
                if(std::abs(d1 - minBondDistanceSq) > 1e-6){
                    continue;
                }
                basis.column(1) = latticeVectors[static_cast<std::size_t>(i1)];

                for(int i2 = i1 + 1; i2 < numNeighbors; ++i2){
                    if(i2 == neighborIndex){
                        continue;
                    }
                    const double d2 = (
                        latticeVectors[static_cast<std::size_t>(neighborIndex)] -
                        latticeVectors[static_cast<std::size_t>(i2)]
                    ).squaredLength();
                    if(std::abs(d2 - minBondDistanceSq) > 1e-6){
                        continue;
                    }
                    basis.column(2) = latticeVectors[static_cast<std::size_t>(i2)];
                    if(std::abs(basis.determinant()) > EPSILON){
                        commonNeighbors[static_cast<std::size_t>(neighborIndex)] = {i1, i2};
                        goto next_neighbor;
                    }
                }
            }
        }

    next_neighbor:
        continue;
    }

    return commonNeighbors;
}

PtmCrystalData buildCrystalData(int structureType){
    PtmCrystalData data;
    const ptm::refdata_t* ref = ptmReferenceForStructureType(structureType);
    if(!ref){
        return data;
    }

    data.coordinationNumber = ref->num_nbrs;
    data.templateToCanonicalNeighborSlot = buildTemplateToCanonicalMapping(structureType, data.coordinationNumber);
    data.latticeVectors = buildCanonicalLatticeVectors(structureType, *ref);
    data.commonNeighbors = buildCommonNeighbors(structureType, data.latticeVectors);

    const auto basisIndices = findNonCoplanarIndices(data.latticeVectors);
    Matrix3 basis = Matrix3::Zero();
    basis.column(0) = data.latticeVectors[static_cast<std::size_t>(basisIndices[0])];
    basis.column(1) = data.latticeVectors[static_cast<std::size_t>(basisIndices[1])];
    basis.column(2) = data.latticeVectors[static_cast<std::size_t>(basisIndices[2])];
    Matrix3 basisInverse = basis.inverse();

    const int count = symmetryMappingCount(*ref);
    data.symmetries.reserve(static_cast<std::size_t>(count));
    for(int mappingIndex = 0; mappingIndex < count; ++mappingIndex){
        PtmSymmetryPermutation symmetry;
        symmetry.permutation = buildSymmetryPermutation(structureType, *ref, mappingIndex);
        symmetry.transformation = Matrix3::Zero();
        symmetry.transformation.column(0) = data.latticeVectors[
            static_cast<std::size_t>(symmetry.permutation[static_cast<std::size_t>(basisIndices[0])])
        ];
        symmetry.transformation.column(1) = data.latticeVectors[
            static_cast<std::size_t>(symmetry.permutation[static_cast<std::size_t>(basisIndices[1])])
        ];
        symmetry.transformation.column(2) = data.latticeVectors[
            static_cast<std::size_t>(symmetry.permutation[static_cast<std::size_t>(basisIndices[2])])
        ];
        symmetry.transformation = symmetry.transformation * basisInverse;

        bool duplicate = false;
        for(const auto& existing : data.symmetries){
            if(existing.transformation.equals(symmetry.transformation)){
                duplicate = true;
                break;
            }
        }

        if(!duplicate){
            data.symmetries.push_back(std::move(symmetry));
        }
    }

    AnalysisSymmetryUtils::calculateSymmetryProducts(data.symmetries);

    return data;
}

const PtmCrystalData& emptyCrystalData(){
    static const PtmCrystalData empty;
    return empty;
}

} // namespace

PtmCrystalInfoProvider::PtmCrystalInfoProvider(){
    initialize(static_cast<int>(StructureType::ICO));
    initialize(static_cast<int>(StructureType::GRAPHENE));
}

std::string_view PtmCrystalInfoProvider::topologyName(int structureType) const{
    switch(normalizedStructureType(structureType)){
        case StructureType::SC:
            return "sc";
        case StructureType::FCC:
            return "fcc";
        case StructureType::HCP:
            return "hcp";
        case StructureType::BCC:
            return "bcc";
        case StructureType::CUBIC_DIAMOND:
            return "cubic_diamond";
        case StructureType::HEX_DIAMOND:
            return "hex_diamond";
        default:
            return {};
    }
}

int PtmCrystalInfoProvider::findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const{
    return AnalysisSymmetryUtils::findClosestSymmetryPermutation(dataFor(structureType).symmetries, rotation);
}

int PtmCrystalInfoProvider::coordinationNumber(int structureType) const{
    return dataFor(structureType).coordinationNumber;
}

int PtmCrystalInfoProvider::commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const{
    const auto& commonNeighbors = dataFor(structureType).commonNeighbors;
    if(neighborIndex < 0 || neighborIndex >= static_cast<int>(commonNeighbors.size()) ||
       commonNeighborSlot < 0 || commonNeighborSlot > 1){
        return -1;
    }
    return commonNeighbors[static_cast<std::size_t>(neighborIndex)][static_cast<std::size_t>(commonNeighborSlot)];
}

int PtmCrystalInfoProvider::symmetryPermutationCount(int structureType) const{
    return static_cast<int>(dataFor(structureType).symmetries.size());
}

int PtmCrystalInfoProvider::symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const{
    const auto& data = dataFor(structureType);
    if(symmetryIndex < 0 || symmetryIndex >= static_cast<int>(data.symmetries.size())){
        return neighborIndex;
    }
    const auto& permutation = data.symmetries[static_cast<std::size_t>(symmetryIndex)].permutation;
    if(neighborIndex < 0 || neighborIndex >= static_cast<int>(permutation.size())){
        return neighborIndex;
    }
    return permutation[static_cast<std::size_t>(neighborIndex)];
}

const Matrix3& PtmCrystalInfoProvider::symmetryTransformation(int structureType, int symmetryIndex) const{
    const auto& data = dataFor(structureType);
    if(symmetryIndex < 0 || symmetryIndex >= static_cast<int>(data.symmetries.size())){
        static const Matrix3 identity = Matrix3::Identity();
        return identity;
    }
    return data.symmetries[static_cast<std::size_t>(symmetryIndex)].transformation;
}

int PtmCrystalInfoProvider::symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const{
    const auto& data = dataFor(structureType);
    if(symmetryIndex < 0 || symmetryIndex >= static_cast<int>(data.symmetries.size())){
        return 0;
    }
    const auto& inverseProduct = data.symmetries[static_cast<std::size_t>(symmetryIndex)].inverseProduct;
    if(transformationIndex < 0 || transformationIndex >= static_cast<int>(inverseProduct.size())){
        return 0;
    }
    return inverseProduct[static_cast<std::size_t>(transformationIndex)];
}

const Vector3& PtmCrystalInfoProvider::latticeVector(int structureType, int latticeVectorIndex) const{
    const auto& data = dataFor(structureType);
    if(latticeVectorIndex < 0 || latticeVectorIndex >= static_cast<int>(data.latticeVectors.size())){
        return zeroVector();
    }
    return data.latticeVectors[static_cast<std::size_t>(latticeVectorIndex)];
}

int PtmCrystalInfoProvider::templateToCanonicalNeighborSlot(int structureType, int templateSlot) const{
    const auto& mapping = dataFor(structureType).templateToCanonicalNeighborSlot;
    if(templateSlot < 0 || templateSlot >= static_cast<int>(mapping.size())){
        return templateSlot;
    }
    return mapping[static_cast<std::size_t>(templateSlot)];
}

void PtmCrystalInfoProvider::initialize(int structureType) const{
    const int normalizedType = normalizedStructureType(structureType);
    if(_data.find(normalizedType) != _data.end()){
        return;
    }

    const StructureType normalized = static_cast<StructureType>(normalizedType);
    if(normalized == StructureType::CUBIC_DIAMOND || normalized == StructureType::HEX_DIAMOND){
        _data[normalizedType] = buildCanonicalDiamondCrystalData(structureType);
        return;
    }

    try{
        _data[normalizedType] = buildCrystalData(structureType);
    }catch(const std::exception&){
        _data[normalizedType] = PtmCrystalData{};
    }
}

const PtmCrystalData& PtmCrystalInfoProvider::dataFor(int structureType) const{
    const int normalizedType = normalizedStructureType(structureType);
    if(_data.find(normalizedType) == _data.end()){
        initialize(normalizedType);
    }
    const auto it = _data.find(normalizedType);
    return it != _data.end() ? it->second : emptyCrystalData();
}

std::shared_ptr<const PtmCrystalInfoProvider> ptmCrystalInfoProviderImpl(){
    static const auto provider = std::make_shared<PtmCrystalInfoProvider>();
    return provider;
}

std::shared_ptr<const StructureAnalysisCrystalInfo> ptmCrystalInfoProvider(){
    return ptmCrystalInfoProviderImpl();
}

int ptmTemplateToCanonicalNeighborSlot(int structureType, int templateSlot){
    return ptmCrystalInfoProviderImpl()->templateToCanonicalNeighborSlot(structureType, templateSlot);
}

}
