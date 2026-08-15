#include "coupled/TrialEvidence.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <unistd.h>

namespace {

void require(const bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void writeBytes(
    const std::filesystem::path& path,
    const void* bytes,
    const std::size_t count
) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    require(output.good(), "failed to create evidence probe payload");
    output.write(static_cast<const char*>(bytes), static_cast<std::streamsize>(count));
    require(output.good(), "failed to write evidence probe payload");
}

coupled::ContentReference referenceFor(
    const std::filesystem::path& root,
    const std::string& name
) {
    const auto path = root / name;
    std::string error;
    const std::string hash = coupled::sha256File(path, &error);
    require(error.empty() && hash.size() == 64u, "failed to hash evidence probe payload");
    return {
        .path = name,
        .sha256 = hash,
        .byteCount = std::filesystem::file_size(path),
    };
}

} // namespace

int main() {
    try {
        const auto root = std::filesystem::temp_directory_path() /
            ("coupled-trial-evidence-" + std::to_string(
                static_cast<unsigned long long>(::getpid())));
        std::filesystem::create_directories(root);

        const std::array<std::uint8_t, 8u> calibration{
            0x43u, 0x41u, 0x4cu, 0x49u, 0x42u, 0x52u, 0x41u, 0x54u,
        };
        const std::array<float, 12u> forceTorque{
            1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f,
            1.5f, 2.5f, 3.5f, 4.5f, 5.5f, 6.5f,
        };
        writeBytes(root / "rig.bin", calibration.data(), calibration.size());
        writeBytes(
            root / "force_torque.f32le",
            forceTorque.data(),
            forceTorque.size() * sizeof(float)
        );

        coupled::CoupledTrialManifest manifest;
        manifest.trialId = "trial-0001";
        manifest.animalId = "animal-001";
        manifest.specimenId = "animal-001-a";
        manifest.pairedSpecimenId = "animal-001-b";
        manifest.cohortSplit = coupled::CohortSplit::blindTransfer;
        manifest.specimenRole = coupled::SpecimenRole::benchIdentification;
        manifest.captureStartUnixNs = 1'786'752'000'000'000'000ll;
        manifest.sourceRevision = "fac8d45ef202cf79e865b719becb20abd4b64a1b";
        manifest.worldFingerprint = 0x1234u;
        manifest.deviceProgramFingerprint = 0x5678u;
        manifest.calibrations.push_back(referenceFor(root, "rig.bin"));
        manifest.streams.push_back({
            .name = "force-torque",
            .content = referenceFor(root, "force_torque.f32le"),
            .scalarType = "f32le",
            .clockId = "capture-clock",
            .channelCount = 6u,
            .strideBytes = 6u * sizeof(float),
            .sampleCount = 2u,
            .startTimeNs = 0,
            .nominalRateHz = 2000.0,
        });

        const auto validation = coupled::validateTrialManifest(manifest);
        require(validation.accepted, validation.message);
        const auto manifestPath = root / "manifest.json";
        const auto write = coupled::writeTrialManifest(manifestPath, manifest);
        require(write.accepted, write.message);

        coupled::CoupledTrialManifest decoded;
        const auto read = coupled::readTrialManifest(manifestPath, decoded);
        require(read.accepted, read.message);
        require(
            coupled::canonicalTrialManifestJSON(decoded) ==
                coupled::canonicalTrialManifestJSON(manifest),
            "canonical trial manifest changed after round trip"
        );
        require(
            coupled::trialManifestFingerprint(decoded) ==
                coupled::trialManifestFingerprint(manifest),
            "trial manifest fingerprint changed after round trip"
        );
        const auto payloads = coupled::verifyTrialPayloads(root, decoded);
        require(payloads.accepted, payloads.message);

        auto traversal = manifest;
        traversal.streams.front().content.path = "../force_torque.f32le";
        require(
            !coupled::validateTrialManifest(traversal).accepted,
            "parent path traversal was accepted"
        );

        coupled::PhysicsEvidence evidence;
        evidence.evidenceId = "evidence-0001";
        evidence.trialManifestSha256 = referenceFor(root, "manifest.json").sha256;
        evidence.sourceRevision = manifest.sourceRevision;
        evidence.deviceName = "Apple M4 Pro";
        evidence.operatingSystem = "macOS 26.6";
        evidence.metallibSha256 = std::string(64u, 'a');
        evidence.worldFingerprint = manifest.worldFingerprint;
        evidence.deviceProgramFingerprint = manifest.deviceProgramFingerprint;
        evidence.replayFingerprint = coupled::trialManifestFingerprint(manifest);
        evidence.environmentSteps = 128u;
        evidence.physicsSubsteps = 1024u;
        evidence.retainedBytes = 1024u;
        evidence.peakBytes = 2048u;
        evidence.elapsedSeconds = 1.0;
        evidence.environmentStepsPerSecond = 128.0;
        evidence.physicsSubstepsPerSecond = 1024.0;
        evidence.exactReplay = true;
        evidence.fp64Parity = true;
        evidence.physicalOutcomesAccepted = true;
        evidence.zeroSwapGrowth = true;
        evidence.timelineCaptured = true;
        evidence.gpuCountersCaptured = false;
        evidence.gpuCounterStatus = "unsupported-by-capture";
        evidence.physicalMetrics.push_back({
            .name = "puncture-force",
            .unit = "N",
            .observed = 1.0,
            .reference = 1.01,
            .absoluteError = 0.01,
            .acceptanceBound = 0.02,
            .accepted = true,
        });
        const auto evidenceValidation = coupled::validatePhysicsEvidence(evidence);
        require(evidenceValidation.accepted, evidenceValidation.message);
        const std::string evidenceJSON = coupled::canonicalPhysicsEvidenceJSON(evidence);
        require(
            evidenceJSON.find("\"exactReplay\":true") != std::string::npos &&
                evidenceJSON.find("\"failedSteps\":0") != std::string::npos,
            "canonical physics evidence omitted absolute gates"
        );

        std::filesystem::remove_all(root);
        std::cout << "trial evidence contract passed: fingerprint="
                  << coupled::trialManifestFingerprint(manifest) << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "trial evidence contract failed: " << error.what() << '\n';
        return 1;
    }
}
