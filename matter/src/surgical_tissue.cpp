#include "numi/matter/surgical_tissue.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace numi::matter {
namespace {

constexpr double kMinimumDimension = 1.0e-9;

[[nodiscard]] bool finitePositive(const JejunalScalar& value) {
    return value.value > 0.0 && std::isfinite(value.value);
}

[[nodiscard]] bool validIdentified(
    const IdentifiedJejunalScalar& value,
    const bool positive
) {
    return value.basis != IdentifiedJejunalBasis::unset &&
        !value.evidenceId.empty() &&
        std::isfinite(value.value) &&
        std::isfinite(value.lowerBound) &&
        std::isfinite(value.upperBound) &&
        value.lowerBound <= value.value &&
        value.value <= value.upperBound &&
        (!positive || value.lowerBound > 0.0);
}

[[nodiscard]] bool validHash(const std::string_view value) {
    return value.size() == 64u && std::ranges::all_of(
        value,
        [](const char character) {
            return (character >= '0' && character <= '9') ||
                (character >= 'a' && character <= 'f');
        }
    );
}

[[nodiscard]] double tetrahedronVolume(
    const std::array<double, 3>& a,
    const std::array<double, 3>& b,
    const std::array<double, 3>& c,
    const std::array<double, 3>& d
) {
    const std::array<double, 3> ab{
        b[0] - a[0], b[1] - a[1], b[2] - a[2],
    };
    const std::array<double, 3> ac{
        c[0] - a[0], c[1] - a[1], c[2] - a[2],
    };
    const std::array<double, 3> ad{
        d[0] - a[0], d[1] - a[1], d[2] - a[2],
    };
    const std::array<double, 3> cross{
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    };
    return (
        cross[0] * ad[0] +
        cross[1] * ad[1] +
        cross[2] * ad[2]
    ) / 6.0;
}

[[nodiscard]] bool setParameter(
    MaterialProgram& material,
    const std::string_view name,
    const double value
) {
    const auto found = std::find_if(
        material.parameters.begin(),
        material.parameters.end(),
        [&](const Parameter& parameter) {
            return parameter.name == name;
        }
    );
    if (found == material.parameters.end() ||
        !std::isfinite(value) ||
        value < found->lower ||
        value > found->upper) {
        return false;
    }
    found->defaultValue = value;
    return true;
}

} // namespace

std::string_view identifiedJejunalBasisName(
    const IdentifiedJejunalBasis basis
) noexcept {
    switch (basis) {
    case IdentifiedJejunalBasis::unset: return "unset";
    case IdentifiedJejunalBasis::trialMeasured: return "trial-measured";
    case IdentifiedJejunalBasis::hierarchicalFit: return "hierarchical-fit";
    case IdentifiedJejunalBasis::derivedGeometry: return "derived-geometry";
    case IdentifiedJejunalBasis::boundedPrior: return "bounded-prior";
    }
    return "invalid";
}

std::string_view jejunalLayerName(const JejunalLayer layer) noexcept {
    switch (layer) {
    case JejunalLayer::mucosaSubmucosa: return "mucosa-submucosa";
    case JejunalLayer::circularMuscle: return "circular-muscle";
    case JejunalLayer::longitudinalMuscle: return "longitudinal-muscle";
    case JejunalLayer::serosa: return "serosa";
    }
    return "invalid";
}

bool validatePerfusedActiveJejunumSpec(
    const PerfusedActiveJejunumSpec& spec,
    std::string* error
) {
    const auto fail = [&](const std::string& message) {
        if (error != nullptr) *error = message;
        return false;
    };
    if (!validHash(spec.trialManifestSha256) ||
        !validHash(spec.parameterBundleSha256)) {
        return fail(
            "perfused-active jejunum requires trial-manifest and parameter-bundle SHA-256 provenance"
        );
    }
    const std::array<const IdentifiedJejunalScalar*, 9> globalScalars{{
        &spec.lengthM,
        &spec.widthM,
        &spec.thicknessM,
        &spec.incisionLengthM,
        &spec.incisionGapM,
        &spec.perfusionTemperatureK,
        &spec.arterialPressurePa,
        &spec.perfusateFlowM3PerS,
        &spec.initialPorePressurePa,
    }};
    for (const auto* scalar : globalScalars) {
        if (!validIdentified(*scalar, true)) {
            return fail(
                "every perfused-active geometry and perfusion scalar requires a positive bounded value and evidence owner"
            );
        }
    }
    if (spec.incisionLengthM.value >= spec.lengthM.value ||
        spec.incisionGapM.value >= spec.widthM.value ||
        spec.longitudinalCells < 4u ||
        spec.longitudinalCells > 256u ||
        spec.circumferentialCells < 4u ||
        spec.circumferentialCells > 256u ||
        spec.circumferentialCells % 2u != 0u ||
        spec.throughThicknessCells < 4u ||
        spec.throughThicknessCells > 64u) {
        return fail(
            "perfused-active jejunum geometry or four-layer mesh resolution is invalid"
        );
    }

    std::array<bool, 4> ownedLayers{};
    double thicknessSum = 0.0;
    for (const auto& layer : spec.layers) {
        const auto layerIndex = static_cast<std::uint32_t>(layer.layer);
        if (layerIndex >= ownedLayers.size() || ownedLayers[layerIndex]) {
            return fail("perfused-active jejunum layers must own each canonical layer exactly once");
        }
        ownedLayers[layerIndex] = true;
        const std::array<const IdentifiedJejunalScalar*, 22> layerScalars{{
            &layer.thicknessFraction,
            &layer.densityKgPerM3,
            &layer.fungCpa,
            &layer.longitudinalCoefficient,
            &layer.circumferentialCoefficient,
            &layer.couplingCoefficient,
            &layer.groundShearPa,
            &layer.viscosityPaS,
            &layer.bulkModulusPa,
            &layer.poreStoragePerPa,
            &layer.poreMobilityM2PerPaS,
            &layer.electricalConductivitySPerM,
            &layer.activationDiffusivityM2PerS,
            &layer.activationOnRatePerS,
            &layer.activationOffRatePerS,
            &layer.maximumActiveTensionPa,
            &layer.activationSlopePerV,
            &layer.cohesiveStrengthPa,
            &layer.fractureEnergyJPerM2,
            &layer.staticFriction,
            &layer.dynamicFriction,
            &layer.activationThresholdV,
        }};
        if (!std::ranges::all_of(layerScalars, [](const auto* value) {
                return validIdentified(*value, false);
            })) {
            return fail("every layer scalar requires bounded provenance");
        }
        const std::array<const IdentifiedJejunalScalar*, 12> strictlyPositive{{
            &layer.thicknessFraction,
            &layer.densityKgPerM3,
            &layer.fungCpa,
            &layer.longitudinalCoefficient,
            &layer.circumferentialCoefficient,
            &layer.groundShearPa,
            &layer.viscosityPaS,
            &layer.bulkModulusPa,
            &layer.poreStoragePerPa,
            &layer.activationSlopePerV,
            &layer.cohesiveStrengthPa,
            &layer.fractureEnergyJPerM2,
        }};
        if (!std::ranges::all_of(strictlyPositive, [](const auto* value) {
                return value->lowerBound > 0.0;
            })) {
            return fail("layer mechanical, storage, damage, and slope parameters must be strictly positive");
        }
        const std::array<const IdentifiedJejunalScalar*, 9> nonnegative{{
            &layer.couplingCoefficient,
            &layer.poreMobilityM2PerPaS,
            &layer.electricalConductivitySPerM,
            &layer.activationDiffusivityM2PerS,
            &layer.activationOnRatePerS,
            &layer.activationOffRatePerS,
            &layer.maximumActiveTensionPa,
            &layer.staticFriction,
            &layer.dynamicFriction,
        }};
        if (!std::ranges::all_of(nonnegative, [](const auto* value) {
                return value->lowerBound >= 0.0;
            })) {
            return fail("layer coupling, transport, activation, and friction parameters must be nonnegative");
        }
        const bool muscle =
            layer.layer == JejunalLayer::circularMuscle ||
            layer.layer == JejunalLayer::longitudinalMuscle;
        if (muscle && layer.maximumActiveTensionPa.lowerBound <= 0.0) {
            return fail("both canonical muscle layers require positive identified active tension");
        }
        if (layer.dynamicFriction.value > layer.staticFriction.value ||
            layer.staticFriction.value > 2.0 ||
            layer.thicknessFraction.value >= 1.0) {
            return fail("layer fractions or friction coefficients are physically inconsistent");
        }
        const double fibreNorm = std::sqrt(
            layer.fibreDirection[0] * layer.fibreDirection[0] +
            layer.fibreDirection[1] * layer.fibreDirection[1] +
            layer.fibreDirection[2] * layer.fibreDirection[2]
        );
        if (!std::isfinite(fibreNorm) || std::abs(fibreNorm - 1.0) > 1.0e-6) {
            return fail("layer fibre direction must be finite and unit length");
        }
        thicknessSum += layer.thicknessFraction.value;
    }
    if (!std::ranges::all_of(ownedLayers, [](const bool value) { return value; }) ||
        std::abs(thicknessSum - 1.0) > 1.0e-6) {
        return fail("four canonical layer thickness fractions must sum to one");
    }
    if (error != nullptr) error->clear();
    return true;
}

bool configurePerfusedActiveJejunumLayerMaterial(
    MaterialProgram& material,
    const PerfusedJejunalLayerSpec& layer,
    std::string* error
) {
    if (jejunalLayerName(layer.layer) == "invalid") {
        if (error != nullptr) *error = "perfused jejunum layer identifier is invalid";
        return false;
    }
    const std::array<const IdentifiedJejunalScalar*, 21> required{{
        &layer.densityKgPerM3,
        &layer.fungCpa,
        &layer.longitudinalCoefficient,
        &layer.circumferentialCoefficient,
        &layer.couplingCoefficient,
        &layer.groundShearPa,
        &layer.viscosityPaS,
        &layer.bulkModulusPa,
        &layer.poreStoragePerPa,
        &layer.poreMobilityM2PerPaS,
        &layer.electricalConductivitySPerM,
        &layer.activationDiffusivityM2PerS,
        &layer.activationOnRatePerS,
        &layer.activationOffRatePerS,
        &layer.maximumActiveTensionPa,
        &layer.activationThresholdV,
        &layer.activationSlopePerV,
        &layer.cohesiveStrengthPa,
        &layer.fractureEnergyJPerM2,
        &layer.staticFriction,
        &layer.dynamicFriction,
    }};
    const double fibreNorm = std::sqrt(
        layer.fibreDirection[0] * layer.fibreDirection[0] +
        layer.fibreDirection[1] * layer.fibreDirection[1] +
        layer.fibreDirection[2] * layer.fibreDirection[2]
    );
    const bool muscle =
        layer.layer == JejunalLayer::circularMuscle ||
        layer.layer == JejunalLayer::longitudinalMuscle;
    if (!std::ranges::all_of(required, [](const auto* value) {
            return validIdentified(*value, false);
        }) ||
        layer.densityKgPerM3.value <= 0.0 ||
        layer.fungCpa.value <= 0.0 ||
        layer.longitudinalCoefficient.value <= 0.0 ||
        layer.circumferentialCoefficient.value <= 0.0 ||
        layer.couplingCoefficient.value < 0.0 ||
        layer.groundShearPa.value <= 0.0 ||
        layer.viscosityPaS.value <= 0.0 ||
        layer.bulkModulusPa.value <= 0.0 ||
        layer.poreStoragePerPa.value <= 0.0 ||
        layer.poreMobilityM2PerPaS.value < 0.0 ||
        layer.electricalConductivitySPerM.value < 0.0 ||
        layer.activationDiffusivityM2PerS.value < 0.0 ||
        layer.activationOnRatePerS.value < 0.0 ||
        layer.activationOffRatePerS.value < 0.0 ||
        layer.maximumActiveTensionPa.value < 0.0 ||
        layer.activationSlopePerV.value <= 0.0 ||
        layer.cohesiveStrengthPa.value <= 0.0 ||
        layer.fractureEnergyJPerM2.value <= 0.0 ||
        layer.staticFriction.value < 0.0 ||
        layer.dynamicFriction.value < 0.0 ||
        layer.dynamicFriction.value > layer.staticFriction.value ||
        (muscle && layer.maximumActiveTensionPa.value <= 0.0) ||
        !std::isfinite(fibreNorm) || std::abs(fibreNorm - 1.0) > 1.0e-6) {
        if (error != nullptr) {
            *error = "perfused jejunum layer contains an unowned or invalid identified value";
        }
        return false;
    }

    MaterialProgram staged = material;
    const bool configured =
        staged.name == "porcine_jejunum_fung" &&
        setParameter(staged, "density", layer.densityKgPerM3.value) &&
        setParameter(staged, "fung_c", layer.fungCpa.value) &&
        setParameter(staged, "a_longitudinal", layer.longitudinalCoefficient.value) &&
        setParameter(staged, "a_circumferential", layer.circumferentialCoefficient.value) &&
        setParameter(staged, "a_coupling", layer.couplingCoefficient.value) &&
        setParameter(staged, "ground_shear", layer.groundShearPa.value) &&
        setParameter(staged, "numerical_viscosity", layer.viscosityPaS.value);
    if (!configured) {
        if (error != nullptr) {
            *error = "material is not the porcine jejunum program or a layer value is outside its compiled range";
        }
        return false;
    }
    staged.name = "porcine_jejunum_perfused_" +
        std::string(jejunalLayerName(layer.layer));
    staged.mixed.bulkModulus = layer.bulkModulusPa.value;
    staged.mixed.poreStorage = layer.poreStoragePerPa.value;
    staged.mixed.poreMobility = layer.poreMobilityM2PerPaS.value;
    staged.mixed.electricalConductivity =
        layer.electricalConductivitySPerM.value;
    staged.mixed.activationDiffusivity =
        layer.activationDiffusivityM2PerS.value;
    staged.mixed.activationOnRate = layer.activationOnRatePerS.value;
    staged.mixed.activationOffRate = layer.activationOffRatePerS.value;
    staged.mixed.maximumActiveTension = layer.maximumActiveTensionPa.value;
    staged.mixed.activationThreshold = layer.activationThresholdV.value;
    staged.mixed.activationSlope = layer.activationSlopePerV.value;
    staged.mixed.cohesiveStrength = layer.cohesiveStrengthPa.value;
    staged.mixed.fractureEnergy = layer.fractureEnergyJPerM2.value;
    staged.mixed.fibreDirection = layer.fibreDirection;
    staged.staticFriction = layer.staticFriction.value;
    staged.dynamicFriction = layer.dynamicFriction.value;
    material = std::move(staged);
    if (error != nullptr) error->clear();
    return true;
}

std::string_view jejunalValueBasisName(
    const JejunalValueBasis basis
) noexcept {
    switch (basis) {
    case JejunalValueBasis::belliniPorcineBiaxialStudy:
        return "sourced:Bellini-porcine-jejunum-biaxial";
    case JejunalValueBasis::derivedGeometry:
        return "derived-from-authored-geometry";
    case JejunalValueBasis::researchDefault:
        return "explicit-research-default";
    }
    return "invalid";
}

std::string_view jejunalValueSourceReference(
    const JejunalValueBasis basis
) noexcept {
    switch (basis) {
    case JejunalValueBasis::belliniPorcineBiaxialStudy:
        return "https://doi.org/10.1016/j.jmbbm.2011.05.030";
    case JejunalValueBasis::derivedGeometry:
        return "structured tetrahedral closure-coupon construction";
    case JejunalValueBasis::researchDefault:
        return "MetalRobo research default; specimen calibration required";
    }
    return "invalid";
}

bool configurePorcineJejunumFungMaterial(
    MaterialProgram& material,
    const PorcineJejunumFungSpec& spec,
    std::string* error
) {
    MaterialProgram staged = material;
    const bool configured =
        staged.name == "porcine_jejunum_fung" &&
        setParameter(staged, "density", spec.densityKgPerM3.value) &&
        setParameter(staged, "fung_c", spec.fungCpa.value) &&
        setParameter(
            staged,
            "a_longitudinal",
            spec.longitudinalCoefficient.value
        ) &&
        setParameter(
            staged,
            "a_circumferential",
            spec.circumferentialCoefficient.value
        ) &&
        setParameter(
            staged,
            "a_coupling",
            spec.couplingCoefficient.value
        ) &&
        setParameter(
            staged,
            "ground_shear",
            spec.groundShearPa.value
        ) &&
        setParameter(
            staged,
            "numerical_viscosity",
            spec.numericalViscosityPaS.value
        ) &&
        finitePositive(spec.bulkModulusPa);
    if (!configured) {
        if (error != nullptr) {
            *error =
                "material is not the expected porcine jejunum Fung "
                "program or a calibrated value is outside its range";
        }
        return false;
    }
    staged.mixed.bulkModulus = spec.bulkModulusPa.value;
    material = std::move(staged);
    if (error != nullptr) {
        error->clear();
    }
    return true;
}

PorcineJejunumFungResponse evaluatePorcineJejunumFungResponse(
    const PorcineJejunumFungSpec& spec,
    const double longitudinalStretch,
    const double circumferentialStretch
) {
    if (!finitePositive(spec.fungCpa) ||
        !finitePositive(spec.longitudinalCoefficient) ||
        !finitePositive(spec.circumferentialCoefficient) ||
        !finitePositive(spec.couplingCoefficient) ||
        !(longitudinalStretch > 0.0) ||
        !(circumferentialStretch > 0.0) ||
        !std::isfinite(longitudinalStretch) ||
        !std::isfinite(circumferentialStretch)) {
        throw std::invalid_argument(
            "porcine jejunum Fung response input is invalid"
        );
    }
    const double longitudinalGreen =
        0.5 * (longitudinalStretch * longitudinalStretch - 1.0);
    const double circumferentialGreen =
        0.5 * (
            circumferentialStretch * circumferentialStretch - 1.0
        );
    const double q =
        spec.longitudinalCoefficient.value *
            longitudinalGreen * longitudinalGreen +
        spec.circumferentialCoefficient.value *
            circumferentialGreen * circumferentialGreen +
        2.0 * spec.couplingCoefficient.value *
            longitudinalGreen * circumferentialGreen;
    if (q > 80.0 || !std::isfinite(q)) {
        throw std::overflow_error(
            "porcine jejunum Fung response exceeds the calibrated domain"
        );
    }
    const double exponential = std::exp(q);
    PorcineJejunumFungResponse result;
    result.energyDensityPa =
        0.5 * spec.fungCpa.value * (exponential - 1.0);
    result.longitudinalSecondPiolaPa =
        spec.fungCpa.value * exponential * (
            spec.longitudinalCoefficient.value * longitudinalGreen +
            spec.couplingCoefficient.value * circumferentialGreen
        );
    result.circumferentialSecondPiolaPa =
        spec.fungCpa.value * exponential * (
            spec.circumferentialCoefficient.value *
                circumferentialGreen +
            spec.couplingCoefficient.value * longitudinalGreen
        );
    return result;
}

namespace {

PorcineJejunumClosureCoupon makeJejunumClosureCoupon(
    const std::uint32_t materialIndex,
    const PorcineJejunumFungSpec& spec,
    const std::vector<double>& zCoordinates,
    const std::vector<std::uint32_t>& slabMaterials,
    const bool mutableTopology,
    const std::string_view objectName,
    const std::string_view fidelityBoundary
) {
    if (!finitePositive(spec.lengthM) ||
        !finitePositive(spec.widthM) ||
        !finitePositive(spec.thicknessM) ||
        !finitePositive(spec.incisionLengthM) ||
        !finitePositive(spec.incisionGapM) ||
        !finitePositive(spec.densityKgPerM3) ||
        spec.incisionLengthM.value >= spec.lengthM.value ||
        spec.incisionGapM.value >=
            spec.widthM.value /
                static_cast<double>(spec.circumferentialCells) ||
        spec.longitudinalCells < 4u ||
        spec.circumferentialCells < 4u ||
        spec.circumferentialCells % 2u != 0u ||
        spec.throughThicknessCells < 1u ||
        spec.longitudinalCells > 256u ||
        spec.circumferentialCells > 256u ||
        spec.throughThicknessCells > 64u ||
        (!zCoordinates.empty() &&
         zCoordinates.size() != spec.throughThicknessCells + 1u) ||
        (!slabMaterials.empty() &&
         slabMaterials.size() != spec.throughThicknessCells)) {
        throw std::invalid_argument(
            "porcine jejunum closure-coupon specification is invalid"
        );
    }

    const std::uint32_t nx = spec.longitudinalCells;
    const std::uint32_t ny = spec.circumferentialCells;
    const std::uint32_t nz = spec.throughThicknessCells;
    const std::uint32_t seamY = ny / 2u;
    const double dx = spec.lengthM.value / static_cast<double>(nx);
    const double dy = spec.widthM.value / static_cast<double>(ny);
    const double uniformDz =
        spec.thicknessM.value / static_cast<double>(nz);
    const auto zPosition = [&](const std::uint32_t iz) {
        return zCoordinates.empty()
            ? -0.5 * spec.thicknessM.value + uniformDz * iz
            : zCoordinates[iz];
    };
    double minimumDz = std::numeric_limits<double>::infinity();
    for (std::uint32_t iz = 0u; iz < nz; ++iz) {
        const double interval = zPosition(iz + 1u) - zPosition(iz);
        if (!(interval > kMinimumDimension) || !std::isfinite(interval)) {
            throw std::invalid_argument(
                "jejunal through-thickness coordinates are not strictly increasing"
            );
        }
        minimumDz = std::min(minimumDz, interval);
    }
    const std::uint32_t incisionCellCount = std::clamp(
        static_cast<std::uint32_t>(std::llround(
            spec.incisionLengthM.value / dx
        )),
        2u,
        nx - 2u
    );
    const std::uint32_t incisionBegin =
        (nx - incisionCellCount) / 2u;
    const std::uint32_t incisionEnd =
        incisionBegin + incisionCellCount;
    const auto longitudinalPosition = [&](const std::uint32_t ix) {
        const double halfLength = 0.5 * spec.lengthM.value;
        const double halfIncision = 0.5 * spec.incisionLengthM.value;
        if (ix <= incisionBegin) {
            return -halfLength +
                (halfLength - halfIncision) *
                    static_cast<double>(ix) /
                    static_cast<double>(incisionBegin);
        }
        if (ix <= incisionEnd) {
            return -halfIncision + spec.incisionLengthM.value *
                static_cast<double>(ix - incisionBegin) /
                static_cast<double>(incisionCellCount);
        }
        return halfIncision +
            (halfLength - halfIncision) *
                static_cast<double>(ix - incisionEnd) /
                static_cast<double>(nx - incisionEnd);
    };

    PorcineJejunumClosureCoupon result;
    result.spec = spec;
    ObjectSource& object = result.object;
    object.name = std::string(objectName);
    object.materialIndex = materialIndex;
    object.representation = Representation::fem;
    object.twoWayCoupling = true;
    object.characteristicLength = std::min({dx, dy, minimumDz});
    object.mixedFEM = true;
    object.mutationPolicy.enabled = mutableTopology;
    object.mutationPolicy.cohesiveFracture = mutableTopology;
    // Automatic puncture remains disabled until needle/tissue impulse has a
    // specimen-specific calibration; explicit puncture commands are reserved.
    object.mutationPolicy.punctureImpulseThreshold = 0.0;

    const std::size_t baseNodeCount =
        static_cast<std::size_t>(nx + 1u) *
        static_cast<std::size_t>(ny + 1u) *
        static_cast<std::size_t>(nz + 1u);
    object.femNodes.reserve(
        baseNodeCount +
        static_cast<std::size_t>(incisionCellCount - 1u) * (nz + 1u)
    );
    std::vector<std::uint32_t> baseIndices(baseNodeCount);
    const auto flat = [=](
        const std::uint32_t ix,
        const std::uint32_t iy,
        const std::uint32_t iz
    ) {
        return (
            static_cast<std::size_t>(iz) * (ny + 1u) + iy
        ) * (nx + 1u) + ix;
    };
    for (std::uint32_t iz = 0u; iz <= nz; ++iz) {
        for (std::uint32_t iy = 0u; iy <= ny; ++iy) {
            for (std::uint32_t ix = 0u; ix <= nx; ++ix) {
                std::array<double, 3> position{
                    longitudinalPosition(ix),
                    -0.5 * spec.widthM.value + dy * iy,
                    zPosition(iz),
                };
                const bool lowerIncisionLip =
                    iy == seamY &&
                    ix > incisionBegin && ix < incisionEnd;
                if (lowerIncisionLip) {
                    position[1] = -0.5 * spec.incisionGapM.value;
                }
                baseIndices[flat(ix, iy, iz)] =
                    static_cast<std::uint32_t>(object.femNodes.size());
                object.femNodes.push_back(position);
            }
        }
    }

    std::vector<std::uint32_t> upperSeam(
        static_cast<std::size_t>(nx + 1u) * (nz + 1u),
        NM_INVALID_INDEX
    );
    const auto seamFlat = [=](
        const std::uint32_t ix,
        const std::uint32_t iz
    ) {
        return static_cast<std::size_t>(iz) * (nx + 1u) + ix;
    };
    for (std::uint32_t ix = incisionBegin + 1u;
         ix < incisionEnd;
         ++ix) {
        for (std::uint32_t iz = 0u; iz <= nz; ++iz) {
            std::array<double, 3> position =
                object.femNodes[baseIndices[flat(ix, seamY, iz)]];
            position[1] = 0.5 * spec.incisionGapM.value;
            upperSeam[seamFlat(ix, iz)] =
                static_cast<std::uint32_t>(object.femNodes.size());
            object.femNodes.push_back(position);
            result.metadata.incisionLipNodePairs.push_back({
                baseIndices[flat(ix, seamY, iz)],
                upperSeam[seamFlat(ix, iz)],
            });
            ++result.metadata.duplicatedIncisionNodeCount;
        }
    }

    const auto node = [&](const std::uint32_t ix,
                          const std::uint32_t iy,
                          const std::uint32_t iz,
                          const bool upperCell) {
        if (upperCell && iy == seamY &&
            ix > incisionBegin && ix < incisionEnd) {
            return upperSeam[seamFlat(ix, iz)];
        }
        return baseIndices[flat(ix, iy, iz)];
    };
    object.tetrahedra.reserve(
        static_cast<std::size_t>(nx) * ny * nz * 6u
    );
    constexpr std::array<std::array<std::uint32_t, 4>, 6> kTets{{
        {{0u, 1u, 3u, 7u}},
        {{0u, 3u, 2u, 7u}},
        {{0u, 2u, 6u, 7u}},
        {{0u, 6u, 4u, 7u}},
        {{0u, 4u, 5u, 7u}},
        {{0u, 5u, 1u, 7u}},
    }};
    double minimumVolume = std::numeric_limits<double>::infinity();
    for (std::uint32_t iz = 0u; iz < nz; ++iz) {
        for (std::uint32_t iy = 0u; iy < ny; ++iy) {
            const bool upperCell = iy >= seamY;
            for (std::uint32_t ix = 0u; ix < nx; ++ix) {
                const std::array<std::uint32_t, 8> corners{{
                    node(ix, iy, iz, upperCell),
                    node(ix + 1u, iy, iz, upperCell),
                    node(ix, iy + 1u, iz, upperCell),
                    node(ix + 1u, iy + 1u, iz, upperCell),
                    node(ix, iy, iz + 1u, upperCell),
                    node(ix + 1u, iy, iz + 1u, upperCell),
                    node(ix, iy + 1u, iz + 1u, upperCell),
                    node(ix + 1u, iy + 1u, iz + 1u, upperCell),
                }};
                for (const auto& local : kTets) {
                    TetrahedronSource tetrahedron{{
                        corners[local[0]],
                        corners[local[1]],
                        corners[local[2]],
                        corners[local[3]],
                    }};
                    tetrahedron.materialIndex = slabMaterials.empty()
                        ? materialIndex : slabMaterials[iz];
                    const double volume = tetrahedronVolume(
                        object.femNodes[tetrahedron.nodes[0]],
                        object.femNodes[tetrahedron.nodes[1]],
                        object.femNodes[tetrahedron.nodes[2]],
                        object.femNodes[tetrahedron.nodes[3]]
                    );
                    if (!(volume > kMinimumDimension * kMinimumDimension *
                                      kMinimumDimension) ||
                        !std::isfinite(volume)) {
                        throw std::logic_error(
                            "porcine jejunum mesh contains an inverted or "
                            "degenerate tetrahedron"
                        );
                    }
                    minimumVolume = std::min(minimumVolume, volume);
                    object.tetrahedra.push_back(tetrahedron);
                }
            }
        }
    }

    // Contact is restricted to the actual tetrahedral boundary. Cooking every
    // node would expose interior quadrature support as a collision surface and
    // let a needle or tool contact material that is still buried in the wall.
    constexpr std::array<std::array<std::uint32_t, 3>, 4> kFaces{{
        {{1u, 2u, 3u}},
        {{0u, 3u, 2u}},
        {{0u, 1u, 3u}},
        {{0u, 2u, 1u}},
    }};
    std::map<std::array<std::uint32_t, 3>, std::uint32_t> faceIncidence;
    for (const TetrahedronSource& tetrahedron : object.tetrahedra) {
        for (const auto& localFace : kFaces) {
            std::array<std::uint32_t, 3> key{
                tetrahedron.nodes[localFace[0]],
                tetrahedron.nodes[localFace[1]],
                tetrahedron.nodes[localFace[2]],
            };
            std::ranges::sort(key);
            std::uint32_t& incidence = faceIncidence[key];
            ++incidence;
            if (incidence > 2u) {
                throw std::logic_error(
                    "porcine jejunum mesh contains a non-manifold face"
                );
            }
        }
    }
    std::vector<bool> boundaryNode(object.femNodes.size(), false);
    for (const auto& [face, incidence] : faceIncidence) {
        if (incidence != 1u) {
            continue;
        }
        for (const std::uint32_t nodeIndex : face) {
            boundaryNode[nodeIndex] = true;
        }
    }
    for (std::uint32_t index = 0u; index < boundaryNode.size(); ++index) {
        if (boundaryNode[index]) {
            object.femContactNodes.push_back(index);
        }
    }
    if (spec.fixLongitudinalEnds) {
        for (std::uint32_t iz = 0u; iz <= nz; ++iz) {
            for (std::uint32_t iy = 0u; iy <= ny; ++iy) {
                object.femFixedNodes.push_back(
                    baseIndices[flat(0u, iy, iz)]
                );
                object.femFixedNodes.push_back(
                    baseIndices[flat(nx, iy, iz)]
                );
            }
        }
        std::ranges::sort(object.femFixedNodes);
        const auto unique = std::ranges::unique(object.femFixedNodes);
        object.femFixedNodes.erase(unique.begin(), unique.end());
    }

    const std::uint32_t nodeCount =
        static_cast<std::uint32_t>(object.femNodes.size());
    const std::uint32_t tetrahedronCount =
        static_cast<std::uint32_t>(object.tetrahedra.size());
    if (mutableTopology) {
        object.femCapacity.nodes = nodeCount + 64u;
        object.femCapacity.tetrahedra = tetrahedronCount + 128u;
        // A conforming tetrahedral volume has fewer than 2*T internal faces.
        // Reserve that topology-derived bound rather than a blanket factor.
        object.femCapacity.cohesiveFaces = 2u * tetrahedronCount;
        object.femCapacity.punctureChannels = 64u;
        object.femCapacity.mutationCommands = 64u;
    }

    result.metadata.nodeCount = nodeCount;
    result.metadata.tetrahedronCount = tetrahedronCount;
    result.metadata.fixedNodeCount =
        static_cast<std::uint32_t>(object.femFixedNodes.size());
    result.metadata.minimumRestTetrahedronVolumeM3 = minimumVolume;
    result.metadata.fidelityBoundary = std::string(fidelityBoundary);
    return result;
}

} // namespace

PorcineJejunumClosureCoupon makePorcineJejunumClosureCoupon(
    const std::uint32_t materialIndex,
    const PorcineJejunumFungSpec& spec
) {
    return makeJejunumClosureCoupon(
        materialIndex,
        spec,
        {},
        {},
        true,
        "porcine_jejunum_enterotomy_coupon",
        "Source-parameterized porcine jejunal planar hyperelasticity; "
        "research 3-D regularization, density, incision and fixture; no "
        "patient, perfusion, layered histology, failure or puncture "
        "calibration."
    );
}

PerfusedActiveJejunumClosureCoupon
makePerfusedActiveJejunumClosureCoupon(
    const std::array<std::uint32_t, 4>& materialIndices,
    const PerfusedActiveJejunumSpec& spec,
    const bool punctureChannels
) {
    std::string validationError;
    if (!validatePerfusedActiveJejunumSpec(spec, &validationError)) {
        throw std::invalid_argument(validationError);
    }
    auto sortedMaterialIndices = materialIndices;
    std::ranges::sort(sortedMaterialIndices);
    if (std::ranges::any_of(materialIndices, [](const std::uint32_t index) {
            return index == NM_INVALID_INDEX;
        }) || std::ranges::adjacent_find(sortedMaterialIndices) !=
            sortedMaterialIndices.end()) {
        throw std::invalid_argument(
            "perfused-active jejunum requires four distinct material indices"
        );
    }

    std::array<const PerfusedJejunalLayerSpec*, 4> canonicalLayers{};
    for (const auto& layer : spec.layers) {
        canonicalLayers[static_cast<std::size_t>(layer.layer)] = &layer;
    }
    const std::uint32_t nz = spec.throughThicknessCells;
    std::array<std::uint32_t, 4> cells{};
    std::array<double, 4> target{};
    std::uint32_t assigned = 0u;
    for (std::size_t layer = 0u; layer < cells.size(); ++layer) {
        target[layer] = canonicalLayers[layer]->thicknessFraction.value * nz;
        cells[layer] = std::max(
            1u, static_cast<std::uint32_t>(std::floor(target[layer]))
        );
        assigned += cells[layer];
    }
    while (assigned > nz) {
        std::size_t selected = cells.size();
        double greatestExcess = -std::numeric_limits<double>::infinity();
        for (std::size_t layer = 0u; layer < cells.size(); ++layer) {
            const double excess = cells[layer] - target[layer];
            if (cells[layer] > 1u && excess > greatestExcess) {
                greatestExcess = excess;
                selected = layer;
            }
        }
        if (selected == cells.size()) {
            throw std::logic_error(
                "four-layer jejunum cell allocation cannot satisfy its mesh budget"
            );
        }
        --cells[selected];
        --assigned;
    }
    while (assigned < nz) {
        std::size_t selected = 0u;
        double greatestDeficit = -std::numeric_limits<double>::infinity();
        for (std::size_t layer = 0u; layer < cells.size(); ++layer) {
            const double deficit = target[layer] - cells[layer];
            if (deficit > greatestDeficit) {
                greatestDeficit = deficit;
                selected = layer;
            }
        }
        ++cells[selected];
        ++assigned;
    }

    std::vector<double> zCoordinates;
    std::vector<std::uint32_t> slabMaterials;
    zCoordinates.reserve(static_cast<std::size_t>(nz) + 1u);
    slabMaterials.reserve(nz);
    double z = -0.5 * spec.thicknessM.value;
    zCoordinates.push_back(z);
    for (std::size_t layer = 0u; layer < cells.size(); ++layer) {
        const double layerThickness = spec.thicknessM.value *
            canonicalLayers[layer]->thicknessFraction.value;
        const double dz = layerThickness / cells[layer];
        for (std::uint32_t cell = 0u; cell < cells[layer]; ++cell) {
            z += dz;
            zCoordinates.push_back(z);
            slabMaterials.push_back(materialIndices[layer]);
        }
    }
    zCoordinates.back() = 0.5 * spec.thicknessM.value;

    PorcineJejunumFungSpec geometry;
    geometry.lengthM = {spec.lengthM.value, JejunalValueBasis::derivedGeometry};
    geometry.widthM = {spec.widthM.value, JejunalValueBasis::derivedGeometry};
    geometry.thicknessM = {
        spec.thicknessM.value, JejunalValueBasis::derivedGeometry
    };
    geometry.incisionLengthM = {
        spec.incisionLengthM.value, JejunalValueBasis::derivedGeometry
    };
    geometry.incisionGapM = {
        spec.incisionGapM.value, JejunalValueBasis::derivedGeometry
    };
    geometry.densityKgPerM3 = {
        canonicalLayers[0]->densityKgPerM3.value,
        JejunalValueBasis::derivedGeometry
    };
    geometry.longitudinalCells = spec.longitudinalCells;
    geometry.circumferentialCells = spec.circumferentialCells;
    geometry.throughThicknessCells = spec.throughThicknessCells;
    geometry.fixLongitudinalEnds = false;

    const std::uint32_t serosa = static_cast<std::uint32_t>(
        JejunalLayer::serosa
    );
    auto built = makeJejunumClosureCoupon(
        materialIndices[serosa],
        geometry,
        zCoordinates,
        slabMaterials,
        false,
        "perfused_active_jejunal_enterotomy_coupon",
        punctureChannels
            ? "Evidence-owned four-layer jejunal geometry, element materials, "
              "mixed transport, activation and mass-conserving puncture "
              "channels; object-level serosa contact; no cohesive mutation "
              "or measured validation."
            : "Evidence-owned four-layer jejunal geometry, element materials, "
              "mixed transport and activation; object-level serosa contact; "
              "no mutation or measured validation."
    );
    built.object.multiphysics.enabled = true;
    built.object.multiphysics.initialTemperature =
        spec.perfusionTemperatureK.value;
    built.object.multiphysics.initialMechanicalPressure = 0.0;
    built.object.multiphysics.initialPorePressure =
        spec.initialPorePressurePa.value;
    built.object.multiphysics.initialElectricPotential = 0.0;
    built.object.multiphysics.initialActivation = 0.0;
    if (punctureChannels) {
        built.object.mutationPolicy.enabled = true;
        built.object.mutationPolicy.cohesiveFracture = false;
        built.object.femCapacity.punctureChannels = 64u;
    }

    PerfusedActiveJejunumClosureCoupon result;
    result.object = std::move(built.object);
    result.spec = spec;
    result.metadata = std::move(built.metadata);
    result.materialIndices = materialIndices;
    const std::uint32_t tetrahedraPerSlab =
        6u * spec.longitudinalCells * spec.circumferentialCells;
    for (std::size_t layer = 0u; layer < cells.size(); ++layer) {
        result.tetrahedronCounts[layer] = cells[layer] * tetrahedraPerSlab;
    }
    return result;
}

} // namespace numi::matter
