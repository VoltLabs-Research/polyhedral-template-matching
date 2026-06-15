#include <volt/matcher/template_crystal_info_provider.h>

#include <numeric>

namespace Volt{

namespace{

const Vector3& zeroVector(){
    static const Vector3 vector = Vector3::Zero();
    return vector;
}

const Matrix3& identityMatrix(){
    static const Matrix3 matrix = Matrix3::Identity();
    return matrix;
}

} // namespace

TemplateCrystalInfoProvider::TemplateCrystalInfoProvider(
    std::shared_ptr<const StructureAnalysisCrystalInfo> builtin,
    const TemplateMatcher& matcher,
    int cationCoordination
): _builtin(std::move(builtin)), _cationCoordination(cationCoordination){
    for(const LoadedTemplate& loaded : matcher.templates()){
        CationCrystalData data;
        data.coordinationNumber = _cationCoordination;
        data.name = loaded.name;
        std::iota(data.identityPermutation.begin(), data.identityPermutation.end(), 0);
        _data.emplace(loaded.structureType, std::move(data));
    }
}

const TemplateCrystalInfoProvider::CationCrystalData* TemplateCrystalInfoProvider::dataFor(int structureType) const{
    const auto it = _data.find(structureType);
    return it != _data.end() ? &it->second : nullptr;
}

int TemplateCrystalInfoProvider::findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const{
    if(!isDefinedType(structureType)){
        return _builtin->findClosestSymmetryPermutation(structureType, rotation);
    }
    return 0;
}

int TemplateCrystalInfoProvider::coordinationNumber(int structureType) const{
    if(!isDefinedType(structureType)){
        return _builtin->coordinationNumber(structureType);
    }
    const auto* data = dataFor(structureType);
    return data ? data->coordinationNumber : 0;
}

int TemplateCrystalInfoProvider::commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const{
    if(!isDefinedType(structureType)){
        return _builtin->commonNeighborIndex(structureType, neighborIndex, commonNeighborSlot);
    }
    return -1;
}

int TemplateCrystalInfoProvider::symmetryPermutationCount(int structureType) const{
    if(!isDefinedType(structureType)){
        return _builtin->symmetryPermutationCount(structureType);
    }
    return dataFor(structureType) ? 1 : 0;
}

int TemplateCrystalInfoProvider::symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const{
    if(!isDefinedType(structureType)){
        return _builtin->symmetryPermutationEntry(structureType, symmetryIndex, neighborIndex);
    }
    return neighborIndex;
}

const Matrix3& TemplateCrystalInfoProvider::symmetryTransformation(int structureType, int symmetryIndex) const{
    if(!isDefinedType(structureType)){
        return _builtin->symmetryTransformation(structureType, symmetryIndex);
    }
    return identityMatrix();
}

int TemplateCrystalInfoProvider::symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const{
    if(!isDefinedType(structureType)){
        return _builtin->symmetryInverseProduct(structureType, symmetryIndex, transformationIndex);
    }
    return 0;
}

const Vector3& TemplateCrystalInfoProvider::latticeVector(int structureType, int latticeVectorIndex) const{
    if(!isDefinedType(structureType)){
        return _builtin->latticeVector(structureType, latticeVectorIndex);
    }
    return zeroVector();
}

std::string_view TemplateCrystalInfoProvider::topologyName(int structureType) const{
    if(!isDefinedType(structureType)){
        return _builtin->topologyName(structureType);
    }
    const auto* data = dataFor(structureType);
    return data ? std::string_view(data->name) : std::string_view{};
}

}
