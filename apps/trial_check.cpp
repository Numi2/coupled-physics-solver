#include "coupled/TrialEvidence.hpp"

#include <filesystem>
#include <iostream>

int main(const int argc, const char* const* argv) {
    if (argc != 2) {
        std::cerr << "usage: coupled_trial_check MANIFEST.json\n";
        return 64;
    }
    const std::filesystem::path manifestPath(argv[1]);
    coupled::CoupledTrialManifest manifest;
    const auto loaded = coupled::readTrialManifest(manifestPath, manifest);
    if (!loaded.accepted) {
        std::cerr << "trial manifest rejected: " << loaded.message << '\n';
        return 2;
    }
    const auto verified = coupled::verifyTrialPayloads(
        manifestPath.parent_path(),
        manifest
    );
    if (!verified.accepted) {
        std::cerr << "trial payload rejected: " << verified.message << '\n';
        return 3;
    }
    std::cout << "trial accepted: id=" << manifest.trialId
              << " split=" << coupled::cohortSplitName(manifest.cohortSplit)
              << " role=" << coupled::specimenRoleName(manifest.specimenRole)
              << " streams=" << manifest.streams.size()
              << " fingerprint="
              << coupled::trialManifestFingerprint(manifest) << '\n';
    return 0;
}
