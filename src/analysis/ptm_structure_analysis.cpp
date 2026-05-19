#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/topology/crystal_coordination_topology.h>
#include <volt/topology/crystal_coordination_topology_init.h>
#include <volt/analysis/nearest_neighbor_finder.h>

#include <volt/analysis/internal/ptm_structure_analysis_detail.h>
#include <volt/analysis/ptm_local_atom_state.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <array>
#include <limits>
#include <memory>
#include <vector>

namespace Volt::PtmStructureAnalysisDetail {

bool setupPTM(StructureContext& context, Volt::PTM& ptm, size_t particleCount, bool collectDefGradient){
    // Why: the deformation gradient is the atom-level analogue of OVITO's
    // `Particles::ElasticDeformationGradientProperty`. We ask PTM for it
    // whenever a PtmLocalAtomState container is supplied so downstream
    // code (atoms.msgpack exporter) can surface it. When no state
    // container is requested (e.g. cluster-only passes from internal
    // helpers) we keep the cheaper path that skips the 3x3 matrix.
    ptm.setCalculateDefGradient(collectDefGradient);
    ptm.setRmsdCutoff(std::numeric_limits<double>::infinity());
    ptm.setInputCrystalStructure(context.inputCrystalType);
    return ptm.prepare(context.positions->constDataPoint3(), particleCount, context.simCell);
}

} // namespace Volt::PtmStructureAnalysisDetail

namespace Volt {

void computeMaximumNeighborDistanceFromPTM(StructureAnalysis& analysis){
    StructureContext& context = analysis.context();
    const size_t N = context.atomCount();
    if(N == 0){
        context.maximumNeighborDistance = 0.0;
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

    Volt::PTM ptm;
    const bool collectPerAtomState = static_cast<bool>(atomStates);
    if(!PtmStructureAnalysisDetail::setupPTM(context, ptm, N, collectPerAtomState)){
        throw std::runtime_error("Error trying to initialize PTM.");
    }

    std::fill(context.neighborCounts->dataInt(),
              context.neighborCounts->dataInt() + context.neighborCounts->size(), 0);
    std::fill(context.structureTypes->dataInt(),
              context.structureTypes->dataInt() + context.structureTypes->size(), LATTICE_OTHER);
    std::fill(
        context.atomAllowedSymmetryMasks->dataInt64(),
        context.atomAllowedSymmetryMasks->dataInt64() + context.atomAllowedSymmetryMasks->size(),
        0
    );

    if(atomStates){
        atomStates->assign(N, PtmLocalAtomState{});
    }

    std::vector<uint64_t> cached(N, 0ull);
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

            context.structureTypes->setInt(i, type);
            context.atomAllowedSymmetryMasks->setInt64(
                i,
                static_cast<std::int64_t>(
                    AnalysisSymmetryUtils::fullSymmetryMask(
                        PtmStructureAnalysisDetail::ptmCrystalInfoProvider()->symmetryPermutationCount(type)
                    )
                )
            );

            const int neighborCount = kernel.numTemplateNeighbors();
            localCounts[i] = neighborCount;
            context.neighborCounts->setInt(i, neighborCount);

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
                const int structureType = context.structureTypes->getInt(i);
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
                            context.structureTypes->setInt(i, cnaOrderedType);
                        }
                        canonicalDiamondShellValid[i] = 1;
                        continue;
                    }
                }
            }
        });

        for(size_t i = 0; i < N; ++i){
            const int structureType = context.structureTypes->getInt(i);
            if(
                structureType != StructureType::CUBIC_DIAMOND &&
                structureType != StructureType::HEX_DIAMOND
            ){
                continue;
            }

            if(canonicalDiamondShellValid[i]){
                continue;
            }

            context.structureTypes->setInt(i, LATTICE_OTHER);
            context.neighborCounts->setInt(i, 0);
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

    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& range){
        PTM::Kernel kernel(ptm);

        for(size_t i = range.begin(); i < range.end(); ++i){
            const int count = localCounts[i];
            if(count == 0){
                continue;
            }

            kernel.cacheNeighbors(i, &cached[i]);
            kernel.identifyStructure(i, cached);

            const int start = offsets[i];

            const int structureType = context.structureTypes->getInt(i);
            if(
                (structureType == StructureType::CUBIC_DIAMOND || structureType == StructureType::HEX_DIAMOND) &&
                count == 16 &&
                canonicalDiamondShellValid.size() == N &&
                canonicalDiamondShellValid[i]
            ){
                for(int neighborSlot = 0; neighborSlot < count; ++neighborSlot){
                    indices[start + neighborSlot] = canonicalDiamondNeighbors[i][static_cast<std::size_t>(neighborSlot)];
                }
                continue;
            }

            bool assignedAllCanonicalSlots = true;
            std::array<unsigned char, MAX_NEIGHBORS> slotAssigned{};
            slotAssigned.fill(0);
            for(int templateSlot = 0; templateSlot < count; ++templateSlot){
                const int canonicalSlot = PtmStructureAnalysisDetail::ptmTemplateToCanonicalNeighborSlot(
                    structureType,
                    templateSlot
                );
                if(canonicalSlot < 0 || canonicalSlot >= count || slotAssigned[static_cast<std::size_t>(canonicalSlot)]){
                    assignedAllCanonicalSlots = false;
                    break;
                }
                indices[start + canonicalSlot] = kernel.getTemplateNeighbor(templateSlot).index;
                slotAssigned[static_cast<std::size_t>(canonicalSlot)] = 1;
            }

            if(!assignedAllCanonicalSlots){
                for(int j = 0; j < count; ++j){
                    indices[start + j] = kernel.getTemplateNeighbor(j).index;
                }
            }
        }
    });

    for(size_t i = 0; i < N; ++i){
        if(context.neighborCounts->getInt(i) == 0){
            context.neighborCounts->setInt(i, localCounts[i]);
        }
    }
}

}
