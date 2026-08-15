#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace coupled {

inline constexpr std::uint32_t kTrialManifestSchemaVersion = 1u;
inline constexpr std::uint32_t kPhysicsEvidenceSchemaVersion = 1u;

enum class CohortSplit : std::uint8_t {
    calibration,
    modelSelection,
    blindTransfer,
};

enum class SpecimenRole : std::uint8_t {
    benchIdentification,
    dvrkClosure,
};

struct ContentReference {
    std::string path;
    std::string sha256;
    std::uint64_t byteCount = 0u;
};

struct TrialStream {
    std::string name;
    ContentReference content;
    std::string scalarType;
    std::string clockId;
    std::uint32_t channelCount = 0u;
    std::uint32_t strideBytes = 0u;
    std::uint64_t sampleCount = 0u;
    std::int64_t startTimeNs = 0;
    double nominalRateHz = 0.0;
};

// Versioned host-side record for one physical capture. Identifiers are
// pseudonymous and contain no animal-owner, operator, or facility metadata.
// Every referenced payload is immutable and content-addressed.
struct CoupledTrialManifest {
    std::uint32_t schemaVersion = kTrialManifestSchemaVersion;
    std::string trialId;
    std::string animalId;
    std::string specimenId;
    std::string pairedSpecimenId;
    CohortSplit cohortSplit = CohortSplit::calibration;
    SpecimenRole specimenRole = SpecimenRole::benchIdentification;
    std::int64_t captureStartUnixNs = 0;
    std::string sourceRevision;
    std::uint64_t worldFingerprint = 0u;
    std::uint64_t deviceProgramFingerprint = 0u;
    std::vector<ContentReference> calibrations;
    std::vector<TrialStream> streams;
};

struct PhysicalMetric {
    std::string name;
    std::string unit;
    double observed = 0.0;
    double reference = 0.0;
    double absoluteError = 0.0;
    double acceptanceBound = 0.0;
    bool accepted = false;
};

// Promotion artifact. A benchmark may be fast and still be rejected: all
// absolute gates are represented independently from throughput telemetry.
struct PhysicsEvidence {
    std::uint32_t schemaVersion = kPhysicsEvidenceSchemaVersion;
    std::string evidenceId;
    std::string trialManifestSha256;
    std::string sourceRevision;
    std::string deviceName;
    std::string operatingSystem;
    std::string metallibSha256;
    std::uint64_t worldFingerprint = 0u;
    std::uint64_t deviceProgramFingerprint = 0u;
    std::uint64_t replayFingerprint = 0u;
    std::uint64_t failedSteps = 0u;
    std::uint64_t environmentSteps = 0u;
    std::uint64_t physicsSubsteps = 0u;
    std::uint64_t retainedBytes = 0u;
    std::uint64_t peakBytes = 0u;
    double elapsedSeconds = 0.0;
    double environmentStepsPerSecond = 0.0;
    double physicsSubstepsPerSecond = 0.0;
    bool exactReplay = false;
    bool fp64Parity = false;
    bool physicalOutcomesAccepted = false;
    bool zeroSwapGrowth = false;
    bool timelineCaptured = false;
    bool gpuCountersCaptured = false;
    std::string gpuCounterStatus;
    std::vector<PhysicalMetric> physicalMetrics;
};

struct EvidenceDiagnostics {
    bool accepted = false;
    std::string message;
};

[[nodiscard]] std::string_view cohortSplitName(CohortSplit value) noexcept;
[[nodiscard]] std::string_view specimenRoleName(SpecimenRole value) noexcept;

[[nodiscard]] EvidenceDiagnostics validateTrialManifest(
    const CoupledTrialManifest& manifest
);

[[nodiscard]] std::string canonicalTrialManifestJSON(
    const CoupledTrialManifest& manifest
);

[[nodiscard]] EvidenceDiagnostics writeTrialManifest(
    const std::filesystem::path& path,
    const CoupledTrialManifest& manifest
);

[[nodiscard]] EvidenceDiagnostics readTrialManifest(
    const std::filesystem::path& path,
    CoupledTrialManifest& manifest
);

[[nodiscard]] EvidenceDiagnostics verifyTrialPayloads(
    const std::filesystem::path& manifestDirectory,
    const CoupledTrialManifest& manifest
);

[[nodiscard]] std::string sha256File(
    const std::filesystem::path& path,
    std::string* error = nullptr
);

[[nodiscard]] std::uint64_t trialManifestFingerprint(
    const CoupledTrialManifest& manifest
);

[[nodiscard]] EvidenceDiagnostics validatePhysicsEvidence(
    const PhysicsEvidence& evidence
);

[[nodiscard]] std::string canonicalPhysicsEvidenceJSON(
    const PhysicsEvidence& evidence
);

} // namespace coupled
