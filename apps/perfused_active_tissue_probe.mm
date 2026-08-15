#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "numi/matter/matter.hpp"
#include "numi/matter/surgical_tissue.hpp"
#include "metalrobo/engine_types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#ifndef NUMI_JEJUNUM_MATERIAL
#define NUMI_JEJUNUM_MATERIAL ""
#endif

#ifndef NUMI_MATTER_METALLIB
#define NUMI_MATTER_METALLIB ""
#endif

namespace {

void require(const bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

numi::matter::IdentifiedJejunalScalar identified(
    const double value,
    const double relativeBound = 0.1
) {
    const double width = std::max(std::abs(value) * relativeBound, 1.0e-12);
    return {
        .value = value,
        .lowerBound = value == 0.0 ? 0.0 : value - width,
        .upperBound = value + width,
        .basis = numi::matter::IdentifiedJejunalBasis::hierarchicalFit,
        .evidenceId = std::string(64u, 'c'),
    };
}

numi::matter::PerfusedJejunalLayerSpec layer(
    const numi::matter::JejunalLayer id,
    const double fraction,
    const double density,
    const double fungC,
    const double activeTension,
    const double conductivity,
    const std::array<double, 3>& fibre
) {
    numi::matter::PerfusedJejunalLayerSpec result;
    result.layer = id;
    result.thicknessFraction = identified(fraction);
    result.densityKgPerM3 = identified(density);
    result.fungCpa = identified(fungC);
    result.longitudinalCoefficient = identified(81.2);
    result.circumferentialCoefficient = identified(72.4);
    result.couplingCoefficient = identified(19.7);
    result.groundShearPa = identified(250.0);
    result.viscosityPaS = identified(5.0);
    // The nodal pressure variable requires one shared scale across layers.
    result.bulkModulusPa = identified(150000.0);
    result.poreStoragePerPa = identified(1.0e-6);
    result.poreMobilityM2PerPaS = identified(1.0e-10);
    result.electricalConductivitySPerM = identified(conductivity);
    result.activationDiffusivityM2PerS = identified(1.0e-6);
    result.activationOnRatePerS = identified(2.0);
    result.activationOffRatePerS = identified(1.0);
    result.maximumActiveTensionPa = identified(activeTension);
    result.activationThresholdV = identified(0.1);
    result.activationSlopePerV = identified(8.0);
    result.cohesiveStrengthPa = identified(10000.0);
    result.fractureEnergyJPerM2 = identified(100.0);
    result.staticFriction = identified(0.45);
    result.dynamicFriction = identified(0.35);
    result.fibreDirection = fibre;
    return result;
}

numi::matter::PerfusedActiveJejunumSpec tissueSpec() {
    numi::matter::PerfusedActiveJejunumSpec spec;
    spec.trialManifestSha256 = std::string(64u, 'a');
    spec.parameterBundleSha256 = std::string(64u, 'b');
    spec.lengthM = identified(0.030);
    spec.widthM = identified(0.024);
    spec.thicknessM = identified(0.0008);
    spec.incisionLengthM = identified(0.016);
    spec.incisionGapM = identified(0.0006);
    spec.perfusionTemperatureK = identified(310.15);
    spec.arterialPressurePa = identified(10000.0);
    spec.perfusateFlowM3PerS = identified(1.0e-6);
    spec.initialPorePressurePa = identified(1000.0);
    spec.layers = {{
        layer(numi::matter::JejunalLayer::mucosaSubmucosa,
              0.35, 1020.0, 760.0, 0.0, 0.08, {1.0, 0.0, 0.0}),
        layer(numi::matter::JejunalLayer::circularMuscle,
              0.25, 1050.0, 820.0, 1800.0, 0.16, {0.0, 1.0, 0.0}),
        layer(numi::matter::JejunalLayer::longitudinalMuscle,
              0.25, 1040.0, 700.0, 1500.0, 0.12, {1.0, 0.0, 0.0}),
        layer(numi::matter::JejunalLayer::serosa,
              0.15, 980.0, 520.0, 0.0, 0.05, {0.0, 1.0, 0.0}),
    }};
    spec.longitudinalCells = 4u;
    spec.circumferentialCells = 4u;
    spec.throughThicknessCells = 4u;
    return spec;
}

std::string diagnostics(const std::vector<numi::matter::Diagnostic>& values) {
    std::string result;
    for (const auto& value : values) {
        if (!result.empty()) result += "; ";
        result += value.message;
    }
    return result;
}

struct PreparedTissue {
    numi::matter::CompiledWorld world;
    numi::matter::PerfusedActiveJejunumClosureCoupon coupon;
};

PreparedTissue prepareTissue() {
    const auto spec = tissueSpec();
    std::string error;
    require(numi::matter::validatePerfusedActiveJejunumSpec(spec, &error), error);
    numi::matter::WorldSource source;
    for (const auto& layerSpec : spec.layers) {
        const auto parsed = numi::matter::parseMatterFile(NUMI_JEJUNUM_MATERIAL);
        require(parsed.succeeded(), diagnostics(parsed.diagnostics));
        auto material = parsed.material;
        require(numi::matter::configurePerfusedActiveJejunumLayerMaterial(
            material, layerSpec, &error), error);
        source.materials.push_back(std::move(material));
    }
    PreparedTissue result;
    result.coupon = numi::matter::makePerfusedActiveJejunumClosureCoupon(
        {0u, 1u, 2u, 3u}, spec
    );
    auto& object = result.coupon.object;
    double minimumX = object.femNodes.front()[0];
    double maximumX = minimumX;
    for (const auto& node : object.femNodes) {
        minimumX = std::min(minimumX, node[0]);
        maximumX = std::max(maximumX, node[0]);
    }
    for (std::uint32_t node = 0u; node < object.femNodes.size(); ++node) {
        const double x = object.femNodes[node][0];
        if (std::abs(x - minimumX) <= 1.0e-12 ||
            std::abs(x - maximumX) <= 1.0e-12) {
            object.femFixedNodes.push_back(node);
        }
    }
    numi::matter::FieldBoundarySource ground;
    ground.node = 0u;
    ground.flags = NM_FIELD_DIRICHLET_ELECTRIC_POTENTIAL;
    ground.stableIdentifier = 1u;
    ground.value = {0.0, 0.0, 0.0, 0.0};
    object.fieldBoundaries.push_back(ground);
    numi::matter::FieldBoundarySource pacing;
    pacing.node = result.coupon.metadata.nodeCount - 1u;
    pacing.flags = NM_FIELD_DIRICHLET_ELECTRIC_POTENTIAL;
    pacing.stableIdentifier = 2u;
    pacing.value = {0.0, 0.0, 1.0, 0.0};
    object.fieldBoundaries.push_back(pacing);

    source.environmentCount = 1u;
    source.frameTimestep = 1.0 / 4000.0;
    source.gravity = {0.0, 0.0, 0.0};
    source.deterministic = true;
    source.mixedSolver.newtonIterations = 12u;
    source.mixedSolver.fgmresRestart = 16u;
    source.mixedSolver.fgmresIterations = 64u;
    source.mixedSolver.lineSearchSteps = 12u;
    source.mixedSolver.relativeResidual = 5.0e-4;
    source.mixedSolver.volumeTolerance = 5.0e-4;
    source.mixedSolver.pressureTolerance = 5.0e-4;
    source.mixedSolver.transportTolerance = 5.0e-4;
    source.objects.push_back(object);
    numi::matter::CompileOptions options;
    options.maximumRateExponent = 4u;
    auto compiled = numi::matter::compileWorld(source, options);
    require(compiled.succeeded(), diagnostics(compiled.diagnostics));
    result.world = std::move(compiled.world);
    return result;
}

struct TissueRun {
    numi::matter::RuntimeStateSnapshot snapshot;
    double gpuMilliseconds = 0.0;
    std::uint32_t completedMicrosteps = 0u;
    std::uint32_t fgmresIterations = 0u;
    float minimumDeterminant = 0.0f;
};

TissueRun runTissue(const numi::matter::CompiledWorld& world) {
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        require(device != nil, "no Apple Metal device is available");
        id<MTLCommandQueue> queue = [device newCommandQueue];
        require(queue != nil, "failed to create perfused tissue queue");
        id<MTLBuffer> worldStatus = [device
            newBufferWithLength:sizeof(MRMetalWorldStatusGPU)
            options:MTLResourceStorageModeShared];
        require(worldStatus != nil, "failed to allocate world status");
        auto* rigidStatus = static_cast<MRMetalWorldStatusGPU*>(
            worldStatus.contents
        );
        *rigidStatus = {};
        rigidStatus->code = MR_STEP_SUCCESS;

        numi::matter::Runtime runtime;
        const auto initialized = runtime.initialize(world, {
            .metallib = NUMI_MATTER_METALLIB,
            .environmentCount = 1u,
            .captureEvents = true,
            .captureDiagnostics = true,
            .automaticIdentification = false,
            .adaptiveTransfer = false,
        });
        require(initialized.encoded && runtime.valid(), initialized.message);
        id<MTLBuffer> matterStatusBuffer =
            (__bridge id<MTLBuffer>)runtime.statusBuffer();
        require(matterStatusBuffer != nil &&
                matterStatusBuffer.contents != nullptr,
            "Matter status buffer is unavailable");
        id<MTLCommandBuffer> commandBuffer = [queue commandBuffer];
        require(commandBuffer != nil, "failed to create command buffer");
        numi::matter::EncodeRequest request{};
        request.commandBuffer = (__bridge void*)commandBuffer;
        request.environmentStatuses = (__bridge void*)worldStatus;
        request.phase = numi::matter::EncodePhase::preDynamics;
        request.controlStep = 0u;
        request.physicsSubstep = 0u;
        request.physicsSubsteps = 1u;
        request.timestepSeconds = runtime.timestepSeconds();
        auto encoded = runtime.encode(request);
        require(encoded.encoded, encoded.message);
        request.phase = numi::matter::EncodePhase::postCommit;
        encoded = runtime.encode(request);
        require(encoded.encoded, encoded.message);
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];
        require(commandBuffer.status == MTLCommandBufferStatusCompleted,
            "perfused tissue command buffer failed");
        const auto* matterStatus = static_cast<const NMMatterStatusGPU*>(
            matterStatusBuffer.contents
        );
        require(rigidStatus->code == MR_STEP_SUCCESS &&
                matterStatus->code == NM_STATUS_SUCCESS,
            "perfused tissue transaction failed: status=" +
                std::to_string(matterStatus->code) + " object=" +
                std::to_string(matterStatus->objectIndex) + " index=" +
                std::to_string(matterStatus->failingIndex));
        TissueRun result;
        result.snapshot = runtime.snapshot();
        require(result.snapshot.available, result.snapshot.message);
        result.gpuMilliseconds = 1000.0 * std::max(
            0.0, commandBuffer.GPUEndTime - commandBuffer.GPUStartTime
        );
        result.completedMicrosteps = matterStatus->completedMicrosteps;
        result.fgmresIterations = matterStatus->fgmresIterations;
        result.minimumDeterminant = matterStatus->diagnostics.x;
        return result;
    }
}

template <typename Value>
bool bitIdentical(
    const std::vector<Value>& first,
    const std::vector<Value>& second
) {
    return first.size() == second.size() &&
        (first.empty() || std::memcmp(
            first.data(), second.data(), first.size() * sizeof(Value)
        ) == 0);
}

std::uint64_t fingerprintFields(
    const std::vector<NMFEMFieldStateGPU>& fields
) {
    constexpr std::uint64_t kPrime = 1099511628211ull;
    std::uint64_t result = 1469598103934665603ull;
    const auto* bytes = reinterpret_cast<const unsigned char*>(fields.data());
    for (std::size_t index = 0u;
         index < fields.size() * sizeof(NMFEMFieldStateGPU); ++index) {
        result = (result ^ bytes[index]) * kPrime;
    }
    return result;
}

} // namespace

int main() {
    try {
        const auto prepared = prepareTissue();
        const TissueRun first = runTissue(prepared.world);
        const TissueRun replay = runTissue(prepared.world);
        require(bitIdentical(first.snapshot.femNodes, replay.snapshot.femNodes) &&
                bitIdentical(first.snapshot.femFields, replay.snapshot.femFields) &&
                bitIdentical(first.snapshot.femTopologyTetrahedra,
                    replay.snapshot.femTopologyTetrahedra),
            "four-layer perfused tissue replay is not bit-identical");
        double maximumActivation = 0.0;
        for (const auto& field : first.snapshot.femFields) {
            maximumActivation = std::max(
                maximumActivation, static_cast<double>(field.secondary.x)
            );
        }
        require(first.completedMicrosteps > 0u &&
                first.fgmresIterations > 0u &&
                first.minimumDeterminant > 0.0f &&
                maximumActivation > 0.0,
            "four-layer perfused tissue did not execute coupled physics");
        std::cout
            << "{\"schema\":\"numi.matter.perfused-active.v1\""
            << ",\"layers\":4"
            << ",\"nodes\":" << prepared.coupon.metadata.nodeCount
            << ",\"tetrahedra\":"
            << prepared.coupon.metadata.tetrahedronCount
            << ",\"field_fingerprint\":"
            << fingerprintFields(first.snapshot.femFields)
            << ",\"maximum_activation\":" << maximumActivation
            << ",\"minimum_J\":" << first.minimumDeterminant
            << ",\"completed_microsteps\":"
            << first.completedMicrosteps
            << ",\"fgmres_iterations\":" << first.fgmresIterations
            << ",\"gpu_ms\":" << first.gpuMilliseconds
            << ",\"mutation_enabled\":false"
            << ",\"failed_steps\":0}\n";
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "perfused-active tissue probe failed: "
                  << exception.what() << '\n';
        return 1;
    }
}
