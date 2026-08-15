#include "numi/matter/matter.hpp"
#include "numi/matter/surgical_tissue.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef NUMI_JEJUNUM_MATERIAL
#define NUMI_JEJUNUM_MATERIAL ""
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
        .lowerBound = value - width,
        .upperBound = value + width,
        .basis = numi::matter::IdentifiedJejunalBasis::hierarchicalFit,
        .evidenceId = std::string(64u, 'b'),
    };
}

numi::matter::PerfusedJejunalLayerSpec layer(
    const numi::matter::JejunalLayer id,
    const double fraction,
    const std::array<double, 3>& fibre
) {
    numi::matter::PerfusedJejunalLayerSpec result;
    result.layer = id;
    result.thicknessFraction = identified(fraction);
    result.densityKgPerM3 = identified(1000.0);
    result.fungCpa = identified(690.0);
    result.longitudinalCoefficient = identified(81.2);
    result.circumferentialCoefficient = identified(72.4);
    result.couplingCoefficient = identified(19.7);
    result.groundShearPa = identified(250.0);
    result.viscosityPaS = identified(5.0);
    result.bulkModulusPa = identified(150000.0);
    result.poreStoragePerPa = identified(1.0e-6);
    result.poreMobilityM2PerPaS = identified(1.0e-10);
    result.electricalConductivitySPerM = identified(0.1);
    result.activationDiffusivityM2PerS = identified(1.0e-6);
    result.activationOnRatePerS = identified(2.0);
    result.activationOffRatePerS = identified(1.0);
    result.maximumActiveTensionPa = identified(1000.0);
    result.activationThresholdV = identified(0.0);
    result.activationSlopePerV = identified(1.0);
    result.cohesiveStrengthPa = identified(10000.0);
    result.fractureEnergyJPerM2 = identified(100.0);
    result.staticFriction = identified(0.45);
    result.dynamicFriction = identified(0.35);
    result.fibreDirection = fibre;
    return result;
}

} // namespace

int main() {
    try {
        numi::matter::PerfusedActiveJejunumSpec missing;
        std::string error;
        require(
            !numi::matter::validatePerfusedActiveJejunumSpec(missing, &error) &&
                !error.empty(),
            "default living-tissue specification was silently accepted"
        );

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
            layer(numi::matter::JejunalLayer::mucosaSubmucosa, 0.25, {1.0, 0.0, 0.0}),
            layer(numi::matter::JejunalLayer::circularMuscle, 0.25, {0.0, 1.0, 0.0}),
            layer(numi::matter::JejunalLayer::longitudinalMuscle, 0.25, {1.0, 0.0, 0.0}),
            layer(numi::matter::JejunalLayer::serosa, 0.25, {0.0, 1.0, 0.0}),
        }};
        spec.longitudinalCells = 18u;
        spec.circumferentialCells = 16u;
        spec.throughThicknessCells = 8u;
        require(
            numi::matter::validatePerfusedActiveJejunumSpec(spec, &error),
            error
        );

        const auto parsed = numi::matter::parseMatterFile(NUMI_JEJUNUM_MATERIAL);
        require(parsed.succeeded(), "failed to parse porcine jejunum base material");
        auto material = parsed.material;
        require(
            numi::matter::configurePerfusedActiveJejunumLayerMaterial(
                material,
                spec.layers[1],
                &error
            ),
            error
        );
        require(
            material.name == "porcine_jejunum_perfused_circular-muscle" &&
                material.mixed.maximumActiveTension > 0.0 &&
                material.mixed.poreMobility > 0.0 &&
                material.mixed.electricalConductivity > 0.0,
            "identified layer did not publish coupled multiphysics"
        );

        auto duplicate = spec;
        duplicate.layers[3].layer = numi::matter::JejunalLayer::circularMuscle;
        require(
            !numi::matter::validatePerfusedActiveJejunumSpec(duplicate, &error),
            "duplicate canonical layer ownership was accepted"
        );

        std::cout << "perfused-active tissue contract passed: layers="
                  << spec.layers.size() << '\n';
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "perfused-active tissue contract failed: "
                  << exception.what() << '\n';
        return 1;
    }
}
