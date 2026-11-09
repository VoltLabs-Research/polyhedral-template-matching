#include <volt/analysis/crystal_symmetry_utils.h>
#include <volt/analysis/ptm.h>
#include <volt/analysis/ptm_crystal_info_provider.h>
#include <volt/analysis/ptm_structure_analysis.h>
#include <volt/topology/crystal_coordination_topology.h>
#include <volt/matcher/template_matcher.h>
#include <volt/topology/crystal_coordination_topology_init.h>
#include <volt/analysis/nearest_neighbor_finder.h>
#include <volt/core/progress.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/enumerable_thread_specific.h>

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Volt{

static std::string structureFamilyNames(int checkFlags){
    static const std::pair<int, StructureType> families[] = {
        {PTM_CHECK_SC, StructureType::SC},
        {PTM_CHECK_FCC, StructureType::FCC},
        {PTM_CHECK_HCP, StructureType::HCP},
        {PTM_CHECK_ICO, StructureType::ICO},
        {PTM_CHECK_BCC, StructureType::BCC},
        {PTM_CHECK_DCUB, StructureType::CUBIC_DIAMOND},
        {PTM_CHECK_DHEX, StructureType::HEX_DIAMOND},
        {PTM_CHECK_GRAPHENE, StructureType::GRAPHENE},
    };

    std::string names;
    for(const auto& [flag, structureType] : families){
        if(checkFlags & flag){
            if(!names.empty()){
                names += ", ";
            }
            names += structureTypeName(structureType);
        }
    }
    return names.empty() ? std::string("none") : names;
}

bool setupPTM(StructureContext& context, PTM& ptm, size_t particleCount, bool collectDefGradient,
              int structureCheckFlags, const int* particleTypes){
    ptm.setCalculateDefGradient(collectDefGradient);
    ptm.setRmsdCutoff(std::numeric_limits<double>::infinity());
    ptm.setInputCrystalStructure(context.inputCrystalType);
    ptm.setStructureCheckFlags(structureCheckFlags);
    if(particleTypes){
        ptm.setParticleTypes(particleTypes);
        ptm.setIdentifyOrdering(true);
    }
    return ptm.prepare(context.positions->constDataPoint3(), particleCount, context.simCell);
}

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
    std::shared_ptr<std::vector<PtmLocalAtomState>> atomStates,
    const TemplateMatcher* templates,
    double cationNeighborRadius,
    int structureCheckFlags,
    const int* particleTypes,
    PtmRmsdCutoffReport* rmsdCutoffReport
) {
    StructureContext& context = analysis.context();
    const size_t N = context.atomCount();
    if(!N){
        return;
    }
    analysis.setCrystalInfoProvider(PtmStructureAnalysisDetail::ptmCrystalInfoProvider());

    PTM ptm;
    const bool collectPerAtomState = static_cast<bool>(atomStates);
    const bool reorderDiamondShells =
        context.inputCrystalType == LATTICE_CUBIC_DIAMOND ||
        context.inputCrystalType == LATTICE_HEX_DIAMOND;
    const bool buildCationNetwork = (templates != nullptr && !templates->empty() && cationNeighborRadius > 0.0);

    // Neighbor finding, neighbor ordering, template matching, diamond shell reordering, the
    // template-atom network and the neighbor lists; a pass this run skips weighs nothing.
    ProgressSteps steps({1.0, 2.0, 5.0, reorderDiamondShells ? 2.0 : 0.0, buildCationNetwork ? 1.0 : 0.0, 0.5});
    const Stopwatch neighborTimer;
    if(!setupPTM(context, ptm, N, collectPerAtomState, structureCheckFlags, particleTypes)){
        throw std::runtime_error("Error trying to initialize PTM.");
    }
    spdlog::info("  Nearest neighbors ready ({} per atom) in {}", PTM::MAX_INPUT_NEIGHBORS, neighborTimer.elapsed());
    steps.next();

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

    std::vector<std::array<int, MAX_NEIGHBORS>> cationNeighbors;

    tbb::enumerable_thread_specific<PTM::Kernel> kernels([&ptm]{ return PTM::Kernel(ptm); });
    static constexpr int kRmsdHistogramBins = 2000;
    static constexpr double kRmsdHistogramBinWidth = 0.001;
    struct RejectedRmsdTally{
        std::size_t count = 0;
        double smallest = std::numeric_limits<double>::infinity();
        std::array<std::size_t, kRmsdHistogramBins + 1> histogram{};
    };
    tbb::enumerable_thread_specific<RejectedRmsdTally> rejectedRmsdTallies;

    ProgressReporter orderingProgress("  Ordering neighbors", N);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& range){
        PTM::Kernel& kernel = kernels.local();
        for(size_t i = range.begin(); i < range.end(); ++i){
            kernel.cacheNeighbors(i, &cached[i]);
        }
        orderingProgress.advance(range.size());
    });
    spdlog::info("  Neighbor ordering done in {}", orderingProgress.elapsed());
    steps.next();

    spdlog::info("  Structure templates: {} (input lattice {})",
                 structureFamilyNames(ptm.ptmCheckFlags()), latticeStructureTypeName(context.inputCrystalType));
    ProgressReporter matchingProgress("  Matching templates", N);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& range){
        PTM::Kernel& kernel = kernels.local();
        double env[PTM_MAX_INPUT_POINTS][3];

        RejectedRmsdTally& rejectedRmsd = rejectedRmsdTallies.local();

        for(size_t i = range.begin(); i < range.end(); ++i){
            StructureType type = kernel.identifyStructure(i, cached);

            const bool builtinMatched = (type != StructureType::OTHER && kernel.rmsd() <= rmsdCutoff);
            if(type != StructureType::OTHER && !builtinMatched){
                const double rejected = kernel.rmsd();
                ++rejectedRmsd.count;
                rejectedRmsd.smallest = std::min(rejectedRmsd.smallest, rejected);
                const int bin = static_cast<int>(rejected / kRmsdHistogramBinWidth);
                ++rejectedRmsd.histogram[
                    static_cast<std::size_t>(std::clamp(bin, 0, kRmsdHistogramBins))
                ];
            }
            double bestRmsd = builtinMatched ? kernel.rmsd() : std::numeric_limits<double>::infinity();
            
            int definedType = -1;
            TemplateMatch templateMatch;
            if(templates && !templates->empty()){
                const int envNbrs = std::min(ptm.cachedNeighborCount(i), PTM_MAX_INPUT_POINTS - 1);
                env[0][0] = env[0][1] = env[0][2] = 0.0;
                for(int n = 0; n < envNbrs; n++){
                    const Vector3 delta = ptm.cachedNeighborDelta(i, n);
                    env[n + 1][0] = delta.x();
                    env[n + 1][1] = delta.y();
                    env[n + 1][2] = delta.z();
                }

                int candidateType = -1;
                TemplateMatch match = templates->matchBest(env, envNbrs + 1, &candidateType);
                if(match.topologicalMatch && candidateType >= 0 && match.rmsd <= rmsdCutoff && match.rmsd < bestRmsd){
                    definedType = candidateType;
                    templateMatch = match;
                    bestRmsd = match.rmsd; 
                }
            }

            if(definedType > 0){
                structureTypesData[i] = definedType;
                allowedSymmetryMasksData[i] = buildCationNetwork
                    ? static_cast<std::int64_t>(AnalysisSymmetryUtils::fullSymmetryMask(1))
                    : 0;
                
                localCounts[i] = 0;
                neighborCountsData[i] = 0;
                correspondenceCodes[i] = 0;

                if(atomStates){
                    auto& atomState = (*atomStates)[i];
                    atomState.orientation = templateMatch.orientation.normalized();
                    atomState.rmsd = templateMatch.rmsd;
                    atomState.deformationGradient = Matrix3::Identity();
                    atomState.interatomicDistance = (templateMatch.scale != 0.0) ? 1.0 / templateMatch.scale : 0.0;
                    atomState.correspondencesCode = 0;
                    atomState.orderingType = static_cast<int>(PTM::OrderingType::ORDERING_NONE);
                    atomState.bestTemplateIndex = 0;
                    atomState.valid = true;
                }

                continue;
            }

            if(!builtinMatched) continue;

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
        matchingProgress.advance(range.size());
    });
    spdlog::info("  Template matching done for {} atoms in {}", N, matchingProgress.elapsed());
    steps.next();

    {
        RejectedRmsdTally rejected;
        for(const auto& local : rejectedRmsdTallies){
            rejected.count += local.count;
            rejected.smallest = std::min(rejected.smallest, local.smallest);
            for(std::size_t bin = 0; bin < rejected.histogram.size(); ++bin){
                rejected.histogram[bin] += local.histogram[bin];
            }
        }

        if(rejected.count > 0){
            const std::size_t medianRank = rejected.count / 2;
            std::size_t seen = 0;
            std::size_t medianBin = 0;
            for(std::size_t bin = 0; bin < rejected.histogram.size(); ++bin){
                seen += rejected.histogram[bin];
                if(seen > medianRank){
                    medianBin = bin;
                    break;
                }
            }
            const double median = (static_cast<double>(medianBin) + 0.5) * kRmsdHistogramBinWidth;

            if(rmsdCutoffReport){
                rmsdCutoffReport->rejectedAtoms = rejected.count;
                rmsdCutoffReport->medianRejectedRmsd = median;
                rmsdCutoffReport->smallestRejectedRmsd = rejected.smallest;
            }
            spdlog::info(
                "  PTM: {} of {} atoms ({:.1f}%) matched a lattice template but were discarded by --rmsd {} "
                "(smallest discarded RMSD {:.4f}, median {:.3f}); raise --rmsd to keep them",
                rejected.count, N, 100.0 * static_cast<double>(rejected.count) / static_cast<double>(N),
                rmsdCutoff, rejected.smallest, median
            );
        }
    }

    if(reorderDiamondShells){
        ProgressReporter diamondProgress("  Reordering diamond neighbor shells", N);
        // Indexing the atoms for the CNA ordering, then the reordering itself.
        ProgressSteps diamondSteps({1.0, 2.0});
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
        diamondSteps.next("Reordering diamond neighbor shells");

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
            diamondProgress.advance(range.size());
        });

        std::size_t orderedDiamondAtoms = 0;
        std::size_t demotedDiamondAtoms = 0;
        for(size_t i = 0; i < N; ++i){
            const int structureType = structureTypesData[i];
            if(
                structureType != StructureType::CUBIC_DIAMOND &&
                structureType != StructureType::HEX_DIAMOND
            ){
                continue;
            }

            if(canonicalDiamondShellValid[i]){
                ++orderedDiamondAtoms;
                continue;
            }

            ++demotedDiamondAtoms;
            structureTypesData[i] = LATTICE_OTHER;
            neighborCountsData[i] = 0;
            localCounts[i] = 0;
        }
        spdlog::info("  Diamond neighbor shells reordered by CNA in {}: {} atoms kept, {} set to Other",
                     diamondProgress.elapsed(), orderedDiamondAtoms, demotedDiamondAtoms);
    }
    steps.next();

    if(buildCationNetwork){
        ProgressReporter cationProgress("  Building template-atom network", N);
        // Indexing the atoms for the cation search, then the search itself.
        ProgressSteps cationSteps({1.0, 2.0});
        cationNeighbors.assign(N, std::array<int, MAX_NEIGHBORS>{});
        for(auto& row : cationNeighbors){
            row.fill(-1);
        }

        const double radiusSq = cationNeighborRadius * cationNeighborRadius;
        static constexpr int kCationCapacity = 64;
        NearestNeighborFinder cationFinder(kCationCapacity);
        if(!cationFinder.prepare(context.positions, context.simCell)){
            throw std::runtime_error("Error preparing cation-network neighbour finder.");
        }
        cationSteps.next("Building template-atom network");
        
        tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto &range){
            NearestNeighborFinder::Query<kCationCapacity> query(cationFinder);
            for(size_t i = range.begin(); i < range.end(); i++){
                const int type = structureTypesData[i];
                if(type < TEMPLATE_STRUCTURE_TYPE_BASE) continue;

                query.findNeighbors(i, false);
                const auto& results = query.results();
                int count = 0;

                for(int n = 0; n < results.size() && count < MAX_NEIGHBORS; n++){
                    const auto &nb = results[n];
                    if(nb.distanceSq > radiusSq) continue;
                    const size_t j = nb.index;

                    if(j == i) continue;    

                    if(structureTypesData[j] != type) continue;

                    cationNeighbors[i][static_cast<std::size_t>(count)] = static_cast<int>(j);
                    count++;
                }

                localCounts[i] = count;
                neighborCountsData[i] = count;
            }
            cationProgress.advance(range.size());
        });
        spdlog::info("  Template-atom network built in {} (radius {} A)", cationProgress.elapsed(), cationNeighborRadius);
    }
    steps.next();

    ProgressReporter neighborListProgress("  Building neighbor lists", N);
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

                if(structureType >= TEMPLATE_STRUCTURE_TYPE_BASE){
                    for(int s = 0; s < count; s++){
                        const int neighborIndex = cationNeighbors[i][static_cast<std::size_t>(s)];
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
                }else if(
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

            neighborListProgress.advance(range.size());
            return maxSoFar;
        },
        [](double left, double right) -> double {
            return std::max(left, right);
        }
    );
    spdlog::info("  Neighbor lists built in {}: {} entries, max neighbor distance {:.3f}",
                 neighborListProgress.elapsed(), totalNeighbors, context.maximumNeighborDistance);
}

}
