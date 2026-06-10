#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/analysis/ptm.h>
#include <volt/analysis/ptm_crystal_info_provider.h>
#include <volt/analysis/ptm_structure_analysis.h>
#include <volt/topology/crystal_coordination_topology.h>
#include <volt/topology/crystal_coordination_topology_init.h>
#include <volt/analysis/nearest_neighbor_finder.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <array>
#include <limits>
#include <memory>
#include <vector>

namespace Volt{
namespace{

bool setupPTM(StructureContext& context, PTM& ptm, size_t particleCount, bool collectDefGradient){
    // Why: the deformation gradient is the atom-level analogue of OVITO's
    // `Particles::ElasticDeformationGradientProperty`. We ask PTM for it
    // whenever a PtmLocalAtomState container is supplied so downstream
    // code (atoms.parquet exporter) can surface it. When no state
    // container is requested (e.g. cluster-only passes from internal
    // helpers) we keep the cheaper path that skips the 3x3 matrix.
    ptm.setCalculateDefGradient(collectDefGradient);
    ptm.setRmsdCutoff(std::numeric_limits<double>::infinity());
    ptm.setInputCrystalStructure(context.inputCrystalType);
    return ptm.prepare(context.positions->constDataPoint3(), particleCount, context.simCell);
}

} // namespace

void computeMaximumNeighborDistanceFromPTM(StructureAnalysis& analysis){
    StructureContext& context = analysis.context();
    const size_t N = context.atomCount();
    if(N == 0){
        context.maximumNeighborDistance = 0.0;
        return;
    }
    if(context.maximumNeighborDistance > 0.0 &&
       context.neighborCounts &&
       context.neighborOffsets &&
       context.neighborIndices){
        return;
    }

    const auto* positions = context.positions->constDataPoint3();
    const auto& inverseMatrix = context.simCell.inverseMatrix();
    const auto& directMatrix = context.simCell.matrix();
    const int* counts = context.neighborCounts->constDataInt();
    const int* offsets = context.neighborOffsets->constDataInt();
    const int* indices = context.neighborIndices->constDataInt();

    double maxDistance = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, N), 0.0,
        [&](const tbb::blocked_range<size_t>& range, double maxSoFar) -> double {
            for(size_t i = range.begin(); i < range.end(); ++i){
                const int neighborCount = counts[i];
                double localMaxDistance = 0.0;
                const int start = offsets[i];

                for(int j = 0; j < neighborCount; ++j){
                    int neighbor = indices[start + j];

                    Vector3 delta = positions[neighbor] - positions[i];
                    double f[3];
                    for(int d = 0; d < 3; ++d){
                        f[d] = inverseMatrix.prodrow(delta, d);
                        f[d] -= std::round(f[d]);
                    }

                    Vector3 minimumImage;
                    minimumImage = directMatrix.column(0) * f[0];
                    minimumImage += directMatrix.column(1) * f[1];
                    minimumImage += directMatrix.column(2) * f[2];
                    double distance = minimumImage.length();

                    if(distance > localMaxDistance){
                        localMaxDistance = distance;
                    }
                }

                if(localMaxDistance > maxSoFar){
                    maxSoFar = localMaxDistance;
                }
            }

            return maxSoFar;
        },
        [](double a, double b) -> double { return std::max(a, b); }
    );

    context.maximumNeighborDistance = maxDistance;
}

void determineLocalStructuresWithPTM(
    StructureAnalysis& analysis,
    double rmsdCutoff,
    std::shared_ptr<std::vector<PtmLocalAtomState>> atomStates
) {
    StructureContext& context = analysis.context();
    const size_t N = context.atomCount();
    if(!N){
        return;
    }
    analysis.setCrystalInfoProvider(PtmStructureAnalysisDetail::ptmCrystalInfoProvider());

    PTM ptm;
    const bool collectPerAtomState = static_cast<bool>(atomStates);
    if(!setupPTM(context, ptm, N, collectPerAtomState)){
        throw std::runtime_error("Error trying to initialize PTM.");
    }

    auto* neighborCountsData = context.neighborCounts->dataInt();
    auto* structureTypesData = context.structureTypes->dataInt();
    auto* allowedSymmetryMasksData = context.atomAllowedSymmetryMasks->dataInt64();
    std::fill(neighborCountsData,
              neighborCountsData + context.neighborCounts->size(), 0);
    std::fill(structureTypesData,
              structureTypesData + context.structureTypes->size(), LATTICE_OTHER);
    std::fill(
        allowedSymmetryMasksData,
        allowedSymmetryMasksData + context.atomAllowedSymmetryMasks->size(),
        0
    );
    context.maximumNeighborDistance = 0.0;

    if(atomStates){
        atomStates->assign(N, PtmLocalAtomState{});
    }

    std::vector<uint64_t> cached(N, 0ull);
    std::vector<uint64_t> correspondenceCodes(N, 0ull);
    std::vector<int> localCounts(N, 0);
    std::vector<std::array<int, MAX_NEIGHBORS>> canonicalDiamondNeighbors;
    std::vector<unsigned char> canonicalDiamondShellValid;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& range){
        PTM::Kernel kernel(ptm);

        for(size_t i = range.begin(); i < range.end(); ++i){
            kernel.cacheNeighbors(i, &cached[i]);
            StructureType type = kernel.identifyStructure(i, cached);
            if(type == StructureType::OTHER || kernel.rmsd() > rmsdCutoff){
                continue;
            }

            structureTypesData[i] = type;
            allowedSymmetryMasksData[i] = static_cast<std::int64_t>(
                AnalysisSymmetryUtils::fullSymmetryMask(
                    PtmStructureAnalysisDetail::ptmCrystalInfoProvider()->symmetryPermutationCount(type)
                )
            );

            const int neighborCount = kernel.numTemplateNeighbors();
            localCounts[i] = neighborCount;
            neighborCountsData[i] = neighborCount;
            correspondenceCodes[i] = kernel.correspondencesCode();

            if(atomStates){
                auto& atomState = (*atomStates)[i];
                atomState.orientation = kernel.orientation().normalized();
                atomState.rmsd = kernel.rmsd();
                atomState.deformationGradient = kernel.deformationGradient();
                atomState.interatomicDistance = kernel.interatomicDistance();
                atomState.correspondencesCode = kernel.correspondencesCode();
                atomState.orderingType = static_cast<int>(kernel.orderingType());
                atomState.bestTemplateIndex = kernel.bestTemplateIndex();
                atomState.valid = true;
            }
        }
    });

    if(
        context.inputCrystalType == LATTICE_CUBIC_DIAMOND ||
        context.inputCrystalType == LATTICE_HEX_DIAMOND
    ){
        ensureCoordinationStructuresInitialized();
        canonicalDiamondNeighbors.resize(N);
        canonicalDiamondShellValid.assign(N, 0);
        for(auto& row : canonicalDiamondNeighbors){
            row.fill(-1);
        }

        auto cnaOrderedStructureTypes = std::make_shared<ParticleProperty>(N, DataType::Int, 1, 0, true);
        CoordinationStructures diamondCoordinationStructures(
            cnaOrderedStructureTypes.get(),
            context.inputCrystalType,
            true,
            context.simCell
        );
        NearestNeighborFinder diamondOrderingFinder(MAX_NEIGHBORS);
        if(!diamondOrderingFinder.prepare(context.positions, context.simCell)){
            throw std::runtime_error("Error trying to prepare PTM diamond ordering finder.");
        }

        tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& range){
            for(size_t i = range.begin(); i < range.end(); ++i){
                const int structureType = structureTypesData[i];
                if(
                    structureType != StructureType::CUBIC_DIAMOND &&
                    structureType != StructureType::HEX_DIAMOND
                ){
                    continue;
                }
                if(localCounts[i] != 16){
                    continue;
                }

                int cnaOrderedCount = 0;
                if(diamondCoordinationStructures.determineLocalStructure(
                    diamondOrderingFinder,
                    static_cast<int>(i),
                    &cnaOrderedCount,
                    canonicalDiamondNeighbors[i].data()
                ) > 0.0){
                    const int cnaOrderedType = cnaOrderedStructureTypes->getInt(i);
                    if(
                        cnaOrderedCount == 16 &&
                        (cnaOrderedType == StructureType::CUBIC_DIAMOND || cnaOrderedType == StructureType::HEX_DIAMOND)
                    ){
                        if(cnaOrderedType != structureType){
                            structureTypesData[i] = cnaOrderedType;
                        }
                        canonicalDiamondShellValid[i] = 1;
                        continue;
                    }
                }
            }
        });

        for(size_t i = 0; i < N; ++i){
            const int structureType = structureTypesData[i];
            if(
                structureType != StructureType::CUBIC_DIAMOND &&
                structureType != StructureType::HEX_DIAMOND
            ){
                continue;
            }

            if(canonicalDiamondShellValid[i]){
                continue;
            }

            structureTypesData[i] = LATTICE_OTHER;
            neighborCountsData[i] = 0;
            localCounts[i] = 0;
        }
    }

    auto* offsets = context.neighborOffsets->dataInt();
    offsets[0] = 0;
    for(size_t i = 0; i < N; ++i){
        offsets[i + 1] = offsets[i] + localCounts[i];
    }

    const size_t totalNeighbors = static_cast<size_t>(offsets[N]);
    context.neighborIndices = std::make_shared<ParticleProperty>(totalNeighbors, DataType::Int, 1, 0, false);
    auto* indices = context.neighborIndices->dataInt();
    const auto* positions = context.positions->constDataPoint3();
    const SimulationCell& simCell = context.simCell;

    context.maximumNeighborDistance = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, N),
        0.0,
        [&](const tbb::blocked_range<size_t>& range, double maxSoFar) -> double {
            std::array<int8_t, PTM_MAX_INPUT_POINTS> correspondences{};
            std::array<int, MAX_NEIGHBORS> resolvedNeighbors{};

            for(size_t i = range.begin(); i < range.end(); ++i){
                const int count = localCounts[i];
                if(count == 0){
                    continue;
                }

                const int start = offsets[i];
                const int structureType = structureTypesData[i];
                double localMaxDistance = 0.0;

                if(
                    (structureType == StructureType::CUBIC_DIAMOND || structureType == StructureType::HEX_DIAMOND) &&
                    count == 16 &&
                    canonicalDiamondShellValid.size() == N &&
                    canonicalDiamondShellValid[i]
                ){
                    for(int s = 0; s < count; ++s){
                        const int neighborIndex = canonicalDiamondNeighbors[i][static_cast<std::size_t>(s)];
                        indices[start + s] = neighborIndex;
                        if(neighborIndex >= 0){
                            const double distance = simCell.wrapVector(
                                positions[static_cast<size_t>(neighborIndex)] - positions[i]
                            ).length();
                            if(distance > localMaxDistance){
                                localMaxDistance = distance;
                            }
                        }
                    }
                }else{
                    resolvedNeighbors.fill(-1);
                    bool assignedAllCanonicalSlots = true;
                    std::array<unsigned char, MAX_NEIGHBORS> slotAssigned{};
                    slotAssigned.fill(0);

                    int decodedTemplateIndex = 0;
                    ptm_decode_correspondences(
                        PTM::toPtmStructureType(static_cast<StructureType>(structureType)),
                        correspondenceCodes[i],
                        correspondences.data(),
                        &decodedTemplateIndex
                    );

                    for(int templateSlot = 0; templateSlot < count; ++templateSlot){
                        const int canonicalSlot = PtmStructureAnalysisDetail::ptmTemplateToCanonicalNeighborSlot(
                            structureType,
                            templateSlot
                        );
                        const int mappedIndex = correspondences[static_cast<std::size_t>(templateSlot + 1)] - 1;
                        const int neighborIndex = ptm.cachedNeighborIndex(i, mappedIndex);
                        if(canonicalSlot < 0 ||
                           canonicalSlot >= count ||
                           slotAssigned[static_cast<std::size_t>(canonicalSlot)] ||
                           neighborIndex < 0){
                            assignedAllCanonicalSlots = false;
                            break;
                        }
                        resolvedNeighbors[static_cast<std::size_t>(canonicalSlot)] = neighborIndex;
                        slotAssigned[static_cast<std::size_t>(canonicalSlot)] = 1;
                    }

                    if(!assignedAllCanonicalSlots){
                        for(int s = 0; s < count; ++s){
                            const int mappedIndex = correspondences[static_cast<std::size_t>(s + 1)] - 1;
                            resolvedNeighbors[static_cast<std::size_t>(s)] = ptm.cachedNeighborIndex(i, mappedIndex);
                        }
                    }

                    for(int s = 0; s < count; ++s){
                        const int neighborIndex = resolvedNeighbors[static_cast<std::size_t>(s)];
                        indices[start + s] = neighborIndex;
                        if(neighborIndex >= 0){
                            const double distance = simCell.wrapVector(
                                positions[static_cast<size_t>(neighborIndex)] - positions[i]
                            ).length();
                            if(distance > localMaxDistance){
                                localMaxDistance = distance;
                            }
                        }
                    }
                }

                if(localMaxDistance > maxSoFar){
                    maxSoFar = localMaxDistance;
                }
            }

            return maxSoFar;
        },
        [](double left, double right) -> double {
            return std::max(left, right);
        }
    );
}

}
