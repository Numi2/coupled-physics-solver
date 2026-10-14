#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "numi/matter/matter.hpp"
#include "numi/matter/surgical_tissue.hpp"
#include "numi/matter/suture_thread.h"
#include "metalrobo/engine_types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

constexpr std::uint32_t kLongitudinalCells = 24u;
constexpr std::uint32_t kTransverseCells = 8u;
constexpr std::uint32_t kThroughThicknessCells = 2u;
constexpr std::uint32_t kExpectedTetrahedra =
    6u * kLongitudinalCells * kTransverseCells * kThroughThicknessCells;

void require(const bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

struct LipPair {
    std::uint32_t lower = 0u;
    std::uint32_t upper = 0u;
    double restX = 0.0;
    double restZ = 0.0;
};

std::vector<LipPair> topLipPairs(
    const numi::matter::SyntheticSkinWoundCoupon& coupon,
    const std::span<const std::uint32_t> stateIndexBySourceNode
) {
    const double topZ = 0.5 * coupon.spec.thicknessM;
    std::vector<LipPair> result;
    for (const auto& pair : coupon.metadata.incisionLipNodePairs) {
        require(pair[0] < stateIndexBySourceNode.size() &&
                    pair[1] < stateIndexBySourceNode.size() &&
                    stateIndexBySourceNode[pair[0]] != NM_INVALID_INDEX &&
                    stateIndexBySourceNode[pair[1]] != NM_INVALID_INDEX,
                "incision lip node has no compiled-state mapping");
        const auto& lower = coupon.object.femNodes.at(pair[0]);
        const auto& upper = coupon.object.femNodes.at(pair[1]);
        if (std::abs(lower[2] - topZ) > 1.0e-12 ||
            std::abs(upper[2] - topZ) > 1.0e-12) {
            continue;
        }
        result.push_back({
            stateIndexBySourceNode[pair[0]],
            stateIndexBySourceNode[pair[1]],
            lower[0],
            lower[2],
        });
    }
    std::sort(result.begin(), result.end(), [](const auto& a, const auto& b) {
        return a.restX < b.restX;
    });
    require(result.size() >= 8u,
            "top incision surface has too few paired lip nodes");
    for (std::size_t index = 1u; index < result.size(); ++index) {
        require(result[index - 1u].restX < result[index].restX,
                "top incision lip samples have duplicate longitudinal nodes");
    }
    return result;
}

std::vector<std::size_t> chooseBitePatch(
    const std::vector<LipPair>& pairs,
    const double biteX
) {
    constexpr double kPatchHalfWidthM = 0.0016;
    std::vector<std::size_t> selected;
    for (std::size_t index = 0u; index < pairs.size(); ++index) {
        if (std::abs(pairs[index].restX - biteX) <= kPatchHalfWidthM) {
            selected.push_back(index);
        }
    }
    require(selected.size() >= 3u,
            "authored suture bite patch has too few surface nodes");
    return selected;
}

double pairGap(
    const LipPair& pair,
    const std::vector<NMFEMNodeStateGPU>& nodes
) {
    const auto& lower = nodes.at(pair.lower).positionAndMass;
    const auto& upper = nodes.at(pair.upper).positionAndMass;
    const double dx = static_cast<double>(upper.x - lower.x);
    const double dy = static_cast<double>(upper.y - lower.y);
    const double dz = static_cast<double>(upper.z - lower.z);
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double maximumDisplacement(
    const std::vector<NMFEMNodeStateGPU>& nodes,
    const std::vector<NMFEMNodeStateGPU>& rest
) {
    require(nodes.size() == rest.size(),
            "skin node state changed its cooked capacity");
    double maximum = 0.0;
    for (std::size_t index = 0u; index < nodes.size(); ++index) {
        const double dx = static_cast<double>(
            nodes[index].positionAndMass.x - rest[index].positionAndMass.x);
        const double dy = static_cast<double>(
            nodes[index].positionAndMass.y - rest[index].positionAndMass.y);
        const double dz = static_cast<double>(
            nodes[index].positionAndMass.z - rest[index].positionAndMass.z);
        maximum = std::max(maximum, std::sqrt(dx * dx + dy * dy + dz * dz));
    }
    return maximum;
}

void run(
    const std::filesystem::path& outputDirectory,
    const std::filesystem::path& materialPath,
    const std::filesystem::path& metallibPath,
    const double biteForcePerSiteN,
    const bool validateOnly,
    const std::uint32_t steps,
    const std::uint32_t newtonIterations,
    const bool thread, const bool threadEnabled, const double timestep
) {
    using namespace numi::matter;
    @autoreleasepool {
        constexpr double kPlannedLoadedForceN = 0.010;
        constexpr double kContactSlopM = 1.0e-5;
        constexpr double kZeroForceDriftLimitM = 1.0e-6;
        require(std::isfinite(biteForcePerSiteN) &&
                    (biteForcePerSiteN == 0.0 ||
                     biteForcePerSiteN == kPlannedLoadedForceN),
                "bite force must be the preregistered 0 N control or 0.010 N load");
        auto material = parseMatterFile(materialPath.string());
        require(material.succeeded(),
                "could not parse the synthetic-skin material");

        SyntheticSkinWoundSpec spec;
        spec.longitudinalCells = kLongitudinalCells;
        spec.transverseCells = kTransverseCells;
        spec.throughThicknessCells = kThroughThicknessCells;
        auto coupon = makeSyntheticSkinWoundCoupon(0u, spec);
        coupon.object.mixedFEM = false; // Direct compressible neo-Hookean formulation retains the authored volumetric energy.
        WorldSource source;
        source.frameTimestep = timestep;
        source.mixedSolver.newtonIterations = newtonIterations;
        source.mixedSolver.fgmresRestart = 128;
        source.mixedSolver.fgmresIterations = 128;
        source.mixedSolver.mechanicalResidualScale = 1.0e4;
        source.gravity = {0.0, 0.0, 0.0};
        source.contactSlop = kContactSlopM;
        source.deterministic = true;
        source.fullTriangleContactManifold = true;
        coupon.object.femCapacity.deformableContacts = 65536u;
        source.materials.push_back(std::move(material.material));
        source.objects.push_back(coupon.object);
        CompileOptions implicitOptions;
        implicitOptions.maximumRateExponent = 0; // Fully implicit FEM: certify the full frame rather than CFL scheduling from parameter upper bounds.
        auto compiled = compileWorld(source, implicitOptions);
        require(compiled.succeeded(),
                "could not compile the synthetic-skin traction world");
        std::cerr << "skin_solver_budget newton=" << compiled.world.mixedSolver.nonlinearIterations.x << " krylov=" << compiled.world.mixedSolver.nonlinearIterations.z << " rate=" << compiled.world.objects.front().solver.x << "\n";
        std::string layoutError;
        require(validateCompiledWorldLayout(compiled.world, &layoutError),
                "invalid synthetic-skin traction world: " + layoutError);
        require(compiled.world.objects.size() == 1u,
                "compiled traction world changed the object count");
        const NMContinuumObjectGPU& descriptor =
            compiled.world.objects.front();
        require(descriptor.stateCount ==
                        coupon.object.femNodes.size() + 64u &&
                    compiled.world.fem.tetrahedra.size() ==
                        coupon.object.tetrahedra.size() + 128u &&
                    static_cast<std::size_t>(descriptor.stateOffset) +
                            descriptor.stateCount <=
                        compiled.world.fem.nodes.size() &&
                    compiled.world.fem.topologyNodes.size() ==
                        compiled.world.fem.nodes.size(),
                "compiled traction world has invalid FEM node capacity");
        std::vector<std::uint32_t> stateIndexBySourceNode(
            coupon.object.femNodes.size(), NM_INVALID_INDEX
        );
        const std::size_t descriptorEnd =
            static_cast<std::size_t>(descriptor.stateOffset) +
            descriptor.stateCount;
        for (std::size_t stateIndex = descriptor.stateOffset;
             stateIndex < descriptorEnd;
             ++stateIndex) {
            const NMFEMTopologyNodeGPU& topology =
                compiled.world.fem.topologyNodes[stateIndex];
            if (topology.identity.y != 0u ||
                (topology.identity.w & NM_TOPOLOGY_ACTIVE) == 0u) {
                continue;
            }
            const std::uint32_t sourceNode = topology.identity.x;
            require(sourceNode < stateIndexBySourceNode.size() &&
                        stateIndexBySourceNode[sourceNode] == NM_INVALID_INDEX,
                    "compiled traction world duplicated an authored FEM node");
            stateIndexBySourceNode[sourceNode] =
                static_cast<std::uint32_t>(stateIndex);
        }
        for (std::size_t sourceNode = 0u;
             sourceNode < stateIndexBySourceNode.size();
             ++sourceNode) {
            const std::uint32_t stateIndex =
                stateIndexBySourceNode[sourceNode];
            require(stateIndex != NM_INVALID_INDEX,
                    "compiled traction world omitted an authored FEM node");
            const auto& rest =
                compiled.world.fem.nodes[stateIndex].restAndFixed;
            const auto& sourceNodePosition = coupon.object.femNodes[sourceNode];
            require(std::abs(static_cast<double>(rest.x) -
                             sourceNodePosition[0]) <= 1.0e-7 &&
                        std::abs(static_cast<double>(rest.y) -
                                 sourceNodePosition[1]) <= 1.0e-7 &&
                        std::abs(static_cast<double>(rest.z) -
                                 sourceNodePosition[2]) <= 1.0e-7,
                    "compiled FEM node mapping changed an authored rest position");
        }
        const std::vector<LipPair> lipPairs = topLipPairs(
            coupon, stateIndexBySourceNode
        );
        const std::vector<std::size_t> proximalPatch = chooseBitePatch(
            lipPairs, -0.004
        );
        const std::vector<std::size_t> distalPatch = chooseBitePatch(
            lipPairs, 0.004
        );

        if (validateOnly) {
            std::cout << "synthetic_skin_wound_traction_preflight=ok"
                      << " source_nodes=" << coupon.object.femNodes.size()
                      << " state_node_capacity="
                      << compiled.world.fem.nodes.size()
                      << " source_tetrahedra="
                      << coupon.object.tetrahedra.size()
                      << " state_tetrahedron_capacity="
                      << compiled.world.fem.tetrahedra.size()
                      << " top_lip_pairs=" << lipPairs.size()
                      << " mapped_lip_nodes=" << 2u * lipPairs.size()
                      << " expected_tetrahedra=" << kExpectedTetrahedra
                      << " proximal_bite_nodes=" << proximalPatch.size()
                      << " distal_bite_nodes=" << distalPatch.size()
                      << " force_ramp=cosine60"
                      << " newton_iterations="
                      << compiled.world.mixedSolver.nonlinearIterations.x
                      << " gpu_dispatch=false\n";
            return;
        }

        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        require(device != nil, "no Metal device is available");
        id<MTLCommandQueue> queue = [device newCommandQueue];
        require(queue != nil, "could not create the skin traction queue");
        id<MTLBuffer> worldStatuses = [device
            newBufferWithLength:sizeof(MRMetalWorldStatusGPU)
            options:MTLResourceStorageModeShared];
        const std::size_t forceValueCount =
            compiled.world.dispatch.femNodeCount * 4u;
        id<MTLBuffer> forceBuffer = [device
            newBufferWithLength:forceValueCount * sizeof(float)
            options:MTLResourceStorageModeShared];
        require(worldStatuses != nil && forceBuffer != nil,
                "could not allocate skin traction inputs");
        auto* worldStatus = static_cast<MRMetalWorldStatusGPU*>(
            worldStatuses.contents
        );
        auto* externalForces = static_cast<float*>(forceBuffer.contents);
        require(worldStatus != nullptr && externalForces != nullptr,
                "skin traction inputs are not CPU visible");

        Runtime runtime;
        const auto initialized = runtime.initialize(
            compiled.world,
            {
                .metallib = metallibPath.string(),
                .environmentCount = 1u,
                .captureEvents = true,
                .captureDiagnostics = true,
                .automaticIdentification = false,
                .adaptiveTransfer = false,
            }
        );
        require(initialized.encoded && runtime.valid(),
                initialized.message);

        std::vector<NMSutureAnchorGPU> threadAnchors;
        std::vector<NMSutureLoopGPU> threadLoops;
        // Prethreaded two-bite stitches spanning the whole incision. Eyelets are 1 mm
        // behind each wound lip, with interpolation fixed in tissue material.
        for (double stitchX : {-0.006, -0.003, 0.0, 0.003, 0.006}) {
            std::vector<std::size_t> patch;
            for(std::size_t i=0;i<lipPairs.size();++i) if(std::abs(lipPairs[i].restX-stitchX)<=0.002) patch.push_back(i);
            require(patch.size()>=3,"thread bite patch has too few nodes");
            double sum=0; for(auto i:patch) sum+=std::exp(-0.5*std::pow((lipPairs[i].restX-stitchX)/0.001,2));
            const auto loop=static_cast<std::uint32_t>(threadLoops.size());
            std::array<std::array<double,3>,4> corners{};
            for(unsigned corner=0;corner<4;++corner) {
                const bool upper=corner==1 || corner==2;
                const bool top=corner<2;
                const double targetY=(upper?1:-1)*(0.5*spec.incisionGapM+0.001);
                const double targetZ=(top?1:-1)*0.5*spec.thicknessM;
                for(auto pairIndex:patch) {
                    const double x=lipPairs[pairIndex].restX;
                    std::size_t below=SIZE_MAX,above=SIZE_MAX;
                    double lo=-1e9,hi=1e9;
                    for(std::size_t n=0;n<coupon.object.femNodes.size();++n) {
                        const auto& v=coupon.object.femNodes[n];
                        if(std::abs(v[0]-x)>1e-10 || std::abs(v[2]-targetZ)>1e-10) continue;
                        if(v[1]<=targetY && v[1]>lo) {lo=v[1];below=n;}
                        if(v[1]>=targetY && v[1]<hi) {hi=v[1];above=n;}
                    }
                    require(below!=SIZE_MAX && above!=SIZE_MAX,"eyelet interpolation missing tissue support");
                    const double wx=std::exp(-0.5*std::pow((x-stitchX)/0.001,2))/sum;
                    const double alpha=hi>lo ? (targetY-lo)/(hi-lo) : 0;
                    for(auto [n,w]: {std::pair{below,wx*(1-alpha)},std::pair{above,wx*alpha}}) {
                        if(w==0) continue;
                        threadAnchors.push_back({{stateIndexBySourceNode[n],loop,corner,0},{static_cast<float>(w),0,0,0}});
                        for(unsigned d=0;d<3;++d) corners[corner][d]+=w*coupon.object.femNodes[n][d];
                    }
                }
            }
            double length=0;
            for(unsigned c=0;c<4;++c) {double norm=0;for(unsigned d=0;d<3;++d) norm+=std::pow(corners[(c+1)%4][d]-corners[c][d],2);length+=std::sqrt(norm);}
            threadLoops.push_back({{static_cast<float>(length),static_cast<float>(length+0.006),0.5f,0.0012f},{},{}});
        }
        id<MTLBuffer> anchors=[device newBufferWithBytes:threadAnchors.data() length:threadAnchors.size()*sizeof(NMSutureAnchorGPU) options:MTLResourceStorageModeShared];
        id<MTLBuffer> loops=[device newBufferWithBytes:threadLoops.data() length:threadLoops.size()*sizeof(NMSutureLoopGPU) options:MTLResourceStorageModeShared];
        id<MTLBuffer> markers=[device newBufferWithLength:threadLoops.size()*64*sizeof(nm_float4) options:MTLResourceStorageModeShared];
        NSError* shaderError=nil;
        id<MTLLibrary> threadLibrary=[device newLibraryWithURL:[NSURL fileURLWithPath:[NSString stringWithUTF8String:metallibPath.string().c_str()]] error:&shaderError];
        const auto pipeline=[&](NSString* name) {return [device newComputePipelineStateWithFunction:[threadLibrary newFunctionWithName:name] error:&shaderError];};
        id<MTLComputePipelineState> threadForces=pipeline(@"numi_matter_metal::nm_suture_thread_forces");
        id<MTLComputePipelineState> threadMarkers=pipeline(@"numi_matter_metal::nm_suture_thread_markers");
        require(anchors && loops && markers && threadForces && threadMarkers,"could not prepare native thread pass");
        const nm_uint4 threadCounts{compiled.world.dispatch.femNodeCount,static_cast<std::uint32_t>(threadAnchors.size()),static_cast<std::uint32_t>(threadLoops.size()),threadEnabled?1u:0u};
        const auto encodeThread=[&](id<MTLCommandBuffer> command, float progress, bool draw) {
            id<MTLComputeCommandEncoder> encoder=[command computeCommandEncoder];
            [encoder setComputePipelineState:draw?threadMarkers:threadForces];
            [encoder setBuffer:(__bridge id<MTLBuffer>)runtime.acceptedFEMNodeBuffer() offset:0 atIndex:0];
            [encoder setBuffer:anchors offset:0 atIndex:1]; [encoder setBuffer:loops offset:0 atIndex:2];
            [encoder setBuffer:draw?markers:forceBuffer offset:0 atIndex:3];
            [encoder setBytes:&threadCounts length:sizeof(threadCounts) atIndex:4];
            if(!draw) [encoder setBytes:&progress length:sizeof(progress) atIndex:5];
            [encoder dispatchThreads:MTLSizeMake(draw?threadLoops.size()*64:threadCounts.x,1,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
            [encoder endEncoding];
        };
        std::filesystem::create_directories(outputDirectory);
        std::ofstream samples(outputDirectory / "traction-samples.csv");
        std::ofstream lips(outputDirectory / "lip-geometry.csv");
        require(samples.good() && lips.good(),
                "could not create the skin traction output files");
        std::ofstream threadSamples(outputDirectory / "thread.csv"), threadGeometry(outputDirectory / "thread-geometry.csv"), tissueGeometry(outputDirectory / "tissue-geometry.csv");
        threadSamples << "step,time_s,stitch,loop_rest_m,tail_rest_m,pulled_m,tension_n,loop_length_m,energy_j,stretch\n";
        threadGeometry << "step,stitch,marker,x_m,y_m,z_m,material_m\n";
        tissueGeometry << "step,node,x_m,y_m,z_m\n";
        threadSamples<<std::setprecision(12);threadGeometry<<std::setprecision(12);tissueGeometry<<std::setprecision(12);
        std::ofstream anchorFile(outputDirectory / "thread-anchors.csv"),loopFile(outputDirectory / "thread-source.csv"),reactionFile(outputDirectory / "thread-forces.csv"),surfaceFile(outputDirectory / "surface.csv");
        anchorFile<<"node,stitch,corner,weight\n"<<std::setprecision(12);
        loopFile<<"stitch,initial_loop_rest_m,total_rest_m,EA_n,pull_m\n"<<std::setprecision(12);
        for(auto a:threadAnchors) anchorFile<<a.identity.x<<','<<a.identity.y<<','<<a.identity.z<<','<<a.weights.x<<'\n';
        for(std::size_t i=0;i<threadLoops.size();++i) {auto m=threadLoops[i].material;loopFile<<i<<','<<m.x<<','<<m.y<<','<<m.z<<','<<m.w<<'\n';}
        reactionFile<<"step,node,fx_n,fy_n,fz_n\n"<<std::setprecision(12);
        surfaceFile<<"a,b,c\n";
        for(auto f:compiled.world.fem.surfaceFaces) if(f.adjacency.y==NM_INVALID_INDEX) {
            const auto t=compiled.world.fem.tetrahedra[f.adjacency.x].nodes;
            const std::array<unsigned,4> ids{t.x,t.y,t.z,t.w};
            std::vector<unsigned> v;for(unsigned c=0;c<4;++c) if(c!=f.sides.x) v.push_back(ids[c]);
            surfaceFile<<v[0]<<','<<v[1]<<','<<v[2]<<'\n';
        }
        samples << "step,time_s,bite_force_n,mean_gap_m,minimum_mean_gap_m,"
                   "center_gap_m,minimum_signed_gap_m,maximum_displacement_m,"
                   "minimum_J,maximum_residual,"
                   "equilibrium_residual,relative_correction,volume_residual,pressure_residual,active_tetrahedra,mass_kg,gpu_ms\n";
        lips << "step,time_s,pair,rest_x_m,rest_z_m,"
                "lower_x_m,lower_y_m,lower_z_m,"
                "upper_x_m,upper_y_m,upper_z_m\n";
        samples << std::setprecision(12);
        lips << std::setprecision(12);

        const std::vector<NMFEMNodeStateGPU>& restNodes =
            compiled.world.fem.nodes;
        double initialMass = 0.0;
        for (const auto& node : restNodes) {
            initialMass += node.positionAndMass.w;
        }
        auto* matterStatuses = static_cast<NMMatterStatusGPU*>(
            ((__bridge id<MTLBuffer>)runtime.statusBuffer()).contents
        );
        require(matterStatuses != nullptr,
                "skin traction Matter status is not CPU visible");

        const std::uint32_t kSteps = steps;
        const std::uint32_t kRampSteps = static_cast<std::uint32_t>(std::lround(0.060/timestep));
        constexpr float kPi = 3.14159265358979323846f;
        double totalGpuMilliseconds = 0.0;
        double finalMeanGap = 0.0;
        double initialMeanGap = 0.0;
        double minimumMeanGap = std::numeric_limits<double>::infinity();
        double maximumMeanGapDrift = 0.0;
        double maximumDisplacementM = 0.0;
        double minimumJ = std::numeric_limits<double>::infinity();
        double minimumSignedGapM = std::numeric_limits<double>::infinity();
        double maximumResidual = 0.0;
        std::uint32_t finalActiveTetrahedra = 0u;
        double finalMass = 0.0;

        const auto setPatchForces = [&](
            const std::vector<std::size_t>& patch,
            const double biteX,
            const float forceMagnitude
        ) {
            std::vector<double> weights(patch.size(), 0.0);
            double weightSum = 0.0;
            for (std::size_t local = 0u; local < patch.size(); ++local) {
                const double distance = lipPairs[patch[local]].restX - biteX;
                weights[local] = std::exp(
                    -0.5 * distance * distance / (0.001 * 0.001)
                );
                weightSum += weights[local];
            }
            require(weightSum > 0.0,
                    "suture traction patch has no positive weight");
            for (std::size_t local = 0u; local < patch.size(); ++local) {
                const LipPair& pair = lipPairs[patch[local]];
                const float force = static_cast<float>(
                    forceMagnitude * weights[local] / weightSum
                );
                externalForces[pair.lower * 4u + 1u] += force;
                externalForces[pair.upper * 4u + 1u] -= force;
            }
        };

        auto lastAccepted = runtime.snapshot();
        require(lastAccepted.available,"initial accepted state unavailable");
        std::ofstream attempts(outputDirectory / "attempts.csv");
        attempts << "step,attempt,newton,krylov,status,residual,correction,volume,pressure,rollback_verified\n" << std::setprecision(12);
        for (std::uint32_t step = 0u; step < kSteps; ++step) {
            std::fill_n(externalForces, forceValueCount, 0.0f);
            float forceScale = 1.0f;
            if (step < kRampSteps) {
                const float rampProgress = static_cast<float>(step) /
                    static_cast<float>(kRampSteps - 1u);
                forceScale = 0.5f * (1.0f - std::cos(kPi * rampProgress));
            }
            const float biteForce =
                static_cast<float>(biteForcePerSiteN) * forceScale;
            if(!thread) {setPatchForces(proximalPatch, -0.004, biteForce);
            setPatchForces(distalPatch, 0.004, biteForce);}

            double stepGpuMilliseconds = 0.0;
            for(unsigned attempt=0;attempt<2;++attempt) {
            const SolverIterationBudgets budgets{attempt==0?newtonIterations:std::max(newtonIterations,24u),attempt==0?128u:256u};
            require(runtime.setSolverIterationBudgets(budgets),"could not set native iteration ceiling");
            *worldStatus = {};
            worldStatus->code = MR_STEP_SUCCESS;
            worldStatus->environment = 0u;
            id<MTLCommandBuffer> commandBuffer = [queue commandBuffer];
            require(commandBuffer != nil,
                    "could not allocate a skin traction command buffer");
            if(thread) encodeThread(commandBuffer,forceScale,false);
            EncodeRequest request{};
            request.commandBuffer = (__bridge void*)commandBuffer;
            request.environmentStatuses = (__bridge void*)worldStatuses;
            if(thread) {request.sutureAnchors=(__bridge void*)anchors;request.sutureLoops=(__bridge void*)loops;
                request.sutureAnchorCount=threadCounts.y;request.sutureLoopCount=threadCounts.z;
                request.suturePullProgress=forceScale;request.sutureForceEnabled=threadEnabled;}
            request.femExternalForces = (__bridge void*)forceBuffer;
            request.femExternalForceCount =
                compiled.world.dispatch.femNodeCount;
            request.phase = EncodePhase::preDynamics;
            request.controlStep = step;
            request.physicsSubstep = 0u;
            request.physicsSubsteps = 1u;
            request.timestepSeconds = runtime.timestepSeconds();
            auto encoded = runtime.encode(request);
            require(encoded.encoded,
                    "skin traction pre-dynamics: " + encoded.message);
            request.phase = EncodePhase::postCommit;
            encoded = runtime.encode(request);
            require(encoded.encoded,
                    "skin traction post-commit: " + encoded.message);
            if(thread) {encodeThread(commandBuffer,forceScale,false);encodeThread(commandBuffer,forceScale,true);}
            [commandBuffer commit];
            [commandBuffer waitUntilCompleted];
            require(commandBuffer.status == MTLCommandBufferStatusCompleted,
                    "skin traction Metal command failed");
            const CFTimeInterval gpuStart = commandBuffer.GPUStartTime;
            const CFTimeInterval gpuEnd = commandBuffer.GPUEndTime;
            const double attemptGpuMilliseconds =
                std::isfinite(gpuStart) && std::isfinite(gpuEnd) &&
                    gpuEnd >= gpuStart
                ? 1000.0 * (gpuEnd - gpuStart)
                : 0.0;
            totalGpuMilliseconds += attemptGpuMilliseconds;
            stepGpuMilliseconds += attemptGpuMilliseconds;
            if (matterStatuses[0].code != NM_STATUS_SUCCESS) {
                const NMMatterStatusGPU status = matterStatuses[0];
                std::cerr << std::setprecision(12) << "failure_certificate="
                          << status.diagnostics.x << "," << status.diagnostics.y
                          << "," << status.diagnostics.z << "," << status.diagnostics.w << "\n";
                const auto rejected = runtime.snapshot();
                const bool rolledBack = rejected.available && rejected.femNodes.size()==lastAccepted.femNodes.size() &&
                    rejected.femFields.size()==lastAccepted.femFields.size() &&
                    std::memcmp(rejected.femNodes.data(),lastAccepted.femNodes.data(),rejected.femNodes.size()*sizeof(NMFEMNodeStateGPU))==0 &&
                    std::memcmp(rejected.femFields.data(),lastAccepted.femFields.data(),rejected.femFields.size()*sizeof(NMFEMFieldStateGPU))==0;
                attempts<<step<<','<<attempt<<','<<budgets.newtonIterations<<','<<budgets.fgmresIterations<<','<<status.code<<','<<status.diagnostics.x<<','<<status.diagnostics.y<<','<<status.diagnostics.z<<','<<status.diagnostics.w<<','<<rolledBack<<'\n';attempts.flush();
                require(rolledBack,"rejected native state did not roll back exactly");
                if(attempt==0 && (status.code==NM_STATUS_NONLINEAR_SOLVER_FAILURE || status.code==NM_STATUS_LINEAR_SOLVER_FAILURE)) {
                    std::cerr<<"skin_solver_retry step="<<step<<" accepted_state_unchanged=true newton=24 krylov=256 tolerances_unchanged=true\n";
                    continue;
                }
                throw std::runtime_error(
                    "skin traction Matter status failed at step " +
                    std::to_string(step) + " with code " +
                    std::to_string(status.code) + " object=" +
                    std::to_string(status.objectIndex) + " failing_index=" +
                    std::to_string(status.failingIndex) +
                    " completed_microsteps=" +
                    std::to_string(status.completedMicrosteps) +
                    " fgmres_iterations=" +
                    std::to_string(status.fgmresIterations) +
                    " contact_count=" + std::to_string(status.contactCount) +
                    " diagnostics=" +
                    std::to_string(status.diagnostics.x) + "," +
                    std::to_string(status.diagnostics.y) + "," +
                    std::to_string(status.diagnostics.z) + "," +
                    std::to_string(status.diagnostics.w)
                );
            }

            attempts<<step<<','<<attempt<<','<<budgets.newtonIterations<<','<<budgets.fgmresIterations<<",0,0,0,0,0,1\n";attempts.flush();
            break;
            }
            const RuntimeStateSnapshot snapshot = runtime.snapshot();
            lastAccepted=snapshot;
            require(snapshot.available,
                    "skin traction snapshot: " + snapshot.message);
            require(snapshot.femNodes.size() == restNodes.size(),
                    "skin traction snapshot changed the FEM node count");
            require(!snapshot.solverCertificates.empty(),
                    "skin traction omitted its solver certificate");
            for(std::size_t n=0;n<snapshot.femNodes.size();++n) {
                const auto& p=snapshot.femNodes[n].positionAndMass;
                tissueGeometry<<step<<','<<n<<','<<p.x<<','<<p.y<<','<<p.z<<'\n';
            }
            if(thread) {
                for(unsigned n=0;n<threadCounts.x;++n) if(externalForces[n*4]!=0 || externalForces[n*4+1]!=0 || externalForces[n*4+2]!=0)
                    reactionFile<<step<<','<<n<<','<<externalForces[n*4]<<','<<externalForces[n*4+1]<<','<<externalForces[n*4+2]<<'\n';
                const auto* data=static_cast<const NMSutureLoopGPU*>(loops.contents);
                const auto* points=static_cast<const nm_float4*>(markers.contents);
                for(std::size_t loop=0;loop<threadLoops.size();++loop) {
                    threadSamples<<step<<','<<timestep*(step+1)<<','<<loop<<','<<data[loop].state.x<<','<<data[loop].state.y<<','<<data[loop].state.z<<','<<data[loop].state.w<<','<<data[loop].diagnostics.x<<','<<data[loop].diagnostics.y<<','<<data[loop].diagnostics.z<<'\n';
                    require(std::abs(data[loop].state.x+data[loop].state.y-threadLoops[loop].material.y)<1e-8,"thread material length not conserved");
                    require(std::isfinite(data[loop].state.w) && data[loop].state.w>=0,"invalid thread tension");
                    for(unsigned m=0;m<64;++m) {const auto& p=points[loop*64+m];threadGeometry<<step<<','<<loop<<','<<m<<','<<p.x<<','<<p.y<<','<<p.z<<','<<p.w<<'\n';}
                }
            }
            const NMSolverCertificateGPU& certificate =
                snapshot.solverCertificates.front();
            require(certificate.validity.w > 0.5f &&
                        std::isfinite(certificate.validity.x) &&
                        certificate.validity.x > 0.20f,
                    "skin traction violated the determinant certificate");
            minimumJ = std::min(
                minimumJ,
                static_cast<double>(certificate.validity.x)
            );
            maximumResidual = std::max(
                maximumResidual,
                static_cast<double>(certificate.nonlinear.x)
            );
            maximumDisplacementM = std::max(
                maximumDisplacementM,
                maximumDisplacement(snapshot.femNodes, restNodes)
            );
            finalMass = 0.0;
            for (const auto& node : snapshot.femNodes) {
                finalMass += node.positionAndMass.w;
            }
            finalActiveTetrahedra = 0u;
            for (const auto& tetrahedron : snapshot.femTopologyTetrahedra) {
                finalActiveTetrahedra +=
                    (tetrahedron.identity.w & NM_OBJECT_ACTIVE) != 0u;
            }

            double meanGap = 0.0;
            double centerGap = 0.0;
            double stepMinimumSignedGapM =
                std::numeric_limits<double>::infinity();
            std::size_t centerIndex = 0u;
            double centerDistance = std::numeric_limits<double>::infinity();
            for (std::size_t pairIndex = 0u;
                 pairIndex < lipPairs.size();
                 ++pairIndex) {
                const LipPair& pair = lipPairs[pairIndex];
                meanGap += pairGap(pair, snapshot.femNodes);
                if (std::abs(pair.restX) < centerDistance) {
                    centerDistance = std::abs(pair.restX);
                    centerIndex = pairIndex;
                }
                const auto& lower = snapshot.femNodes[pair.lower]
                                        .positionAndMass;
                const auto& upper = snapshot.femNodes[pair.upper]
                                        .positionAndMass;
                const double signedGapM =
                    static_cast<double>(upper.y - lower.y);
                require(signedGapM >= 0.0,
                        "synthetic wound lips crossed at step " +
                            std::to_string(step));
                stepMinimumSignedGapM = std::min(
                    stepMinimumSignedGapM,
                    signedGapM
                );
                lips << step << ',' << (timestep * (step+1)) << ',' << pairIndex
                     << ',' << pair.restX << ',' << pair.restZ << ','
                     << lower.x << ',' << lower.y << ',' << lower.z << ','
                     << upper.x << ',' << upper.y << ',' << upper.z << '\n';
            }
            meanGap /= static_cast<double>(lipPairs.size());
            centerGap = pairGap(lipPairs[centerIndex], snapshot.femNodes);
            if (step == 0u) {
                initialMeanGap = meanGap;
            }
            minimumMeanGap = std::min(minimumMeanGap, meanGap);
            maximumMeanGapDrift = std::max(
                maximumMeanGapDrift,
                std::abs(meanGap - initialMeanGap)
            );
            minimumSignedGapM = std::min(
                minimumSignedGapM,
                stepMinimumSignedGapM
            );
            finalMeanGap = meanGap;
            samples << step << ',' << (timestep * (step + 1u)) << ','
                    << (thread?0.0f:biteForce) << ',' << meanGap << ',' << minimumMeanGap
                    << ',' << centerGap << ',' << stepMinimumSignedGapM << ','
                    << maximumDisplacementM << ',' << minimumJ << ','
                    << maximumResidual << ',' << certificate.nonlinear.x << ','
                    << certificate.nonlinear.y << ',' << certificate.nonlinear.z << ','
                    << certificate.nonlinear.w << ',' << finalActiveTetrahedra << ','
                    << finalMass << ',' << stepGpuMilliseconds << '\n';
            samples.flush(); lips.flush();
        }

        const double relativeMassError = initialMass > 0.0
            ? std::abs(finalMass - initialMass) / initialMass
            : 1.0;
        require(finalActiveTetrahedra == coupon.object.tetrahedra.size(),
                "skin traction removed or deactivated tetrahedra");
        require(relativeMassError <= 1.0e-5,
                "skin traction did not conserve FEM mass");
        if ((!thread && biteForcePerSiteN == 0.0) || (thread && !threadEnabled)) {
            require(maximumMeanGapDrift <= kZeroForceDriftLimitM,
                    "zero-force control drift exceeded 1 um");
        }
        std::cout << std::setprecision(9)
                  << "synthetic_skin_wound_traction=ok"
                  << " arm=" << (thread ? (threadEnabled?"prethreaded_pullthrough":"no_thread_control") : (biteForcePerSiteN == 0.0
                          ? "zero_force_control"
                          : "paired_bite_traction"))
                  << " mesh=" << kLongitudinalCells << 'x'
                  << kTransverseCells << 'x' << kThroughThicknessCells
                  << " top_lip_pairs=" << lipPairs.size()
                  << " prescribed_patch_force_per_site_n=" << (thread?0.0:biteForcePerSiteN)
                  << " force_ramp=cosine60"
                  << " newton_iterations="
                  << compiled.world.mixedSolver.nonlinearIterations.x
                  << " initial_mean_gap_m=" << initialMeanGap
                  << " minimum_mean_gap_m=" << minimumMeanGap
                  << " final_mean_gap_m=" << finalMeanGap
                  << " maximum_mean_gap_drift_m=" << maximumMeanGapDrift
                  << " minimum_signed_gap_m=" << minimumSignedGapM
                  << " peak_gap_reduction_fraction="
                  << (initialMeanGap > 0.0
                          ? (initialMeanGap - minimumMeanGap) / initialMeanGap
                          : 0.0)
                  << " maximum_tissue_displacement_m="
                  << maximumDisplacementM
                  << " minimum_J=" << minimumJ
                  << " maximum_solver_residual=" << maximumResidual
                  << " active_tetrahedra=" << finalActiveTetrahedra
                  << " mass_relative_error=" << relativeMassError
                  << " total_gpu_ms=" << totalGpuMilliseconds
                  << " physical_calibration=none"
                  << " thread_present=" << (thread?"true":"false")
                  << " thread_coupling=" << (threadEnabled?"on":"off")
                  << " timestep_s=" << timestep
                  << " robot_execution=false\n";
    }
}

} // namespace

int main(const int argc, char* argv[]) {
    try {
        if (argc < 4) {
            std::cout << "usage: coupled_skin_suture OUTPUT MATERIAL METALLIB [--force-per-site-n 0|0.01] [--thread] [--no-thread-coupling] [--steps 1..240] [--timestep 0.001|0.0005] [--newton 1..32] [--validate-only]\n";
            return argc == 2 && std::string_view{argv[1]} == "--help" ? 0 : 2;
        }
        double force = 0.010;
        bool validate = false;
        std::uint32_t steps = 60, newton = 7; bool thread=false,threadEnabled=true; double timestep=0.001;
        for (int i=4;i<argc;++i) {
            const std::string flag{argv[i]};
            if(flag=="--thread") {thread=true;continue;}
            if(flag=="--no-thread-coupling") {threadEnabled=false;continue;}
            if (flag == "--validate-only") {validate=true;continue;}
            require(i+1<argc,"missing option value");
            const std::string value{argv[++i]}; std::size_t parsed=0;
            const double number=std::stod(value,&parsed);
            require(parsed==value.size() && std::isfinite(number),"invalid number");
            if(flag=="--timestep") {require(number==0.001 || number==0.0005,"timestep must be 1 or 0.5 ms");timestep=number;}
            else if (flag=="--force-per-site-n") force=number;
            else if (flag=="--steps" || flag=="--newton") {
                require(number>=1 && number==std::floor(number) && number<=(flag=="--steps"?240:32),"invalid step/iteration count");
                (flag=="--steps"?steps:newton)=static_cast<std::uint32_t>(number);
            } else throw std::runtime_error("unknown option: "+flag);
        }
        run(argv[1],argv[2],argv[3],force,validate,steps,newton,thread,threadEnabled,timestep);
        return 0;
    } catch(const std::exception& e) { std::cerr << "skin_suture=failed reason=" << e.what() << "\n";return 1; }
}
