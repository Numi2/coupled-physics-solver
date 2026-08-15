#include "metalrobo/InteractionPack.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string_view>
#include <vector>

namespace metalrobo {
namespace {

bool countFits(const std::size_t count) noexcept {
    return count <= std::numeric_limits<std::uint32_t>::max();
}

bool stringFits(const std::string_view value) noexcept {
    return countFits(value.size());
}

} // namespace

bool validInteractionPack(const InteractionPack& pack) noexcept {
    try {
        if (pack.id.empty() || pack.sourceRepository.empty() ||
            pack.sourceRevision.empty() || pack.license.empty() ||
            pack.coordinateFrame != kInteractionCoordinateFrame ||
            pack.jointNames.empty() || pack.clips.empty() ||
            !stringFits(pack.id) || !stringFits(pack.sourceRepository) ||
            !stringFits(pack.sourceRevision) || !stringFits(pack.license) ||
            !countFits(pack.jointNames.size()) ||
            !countFits(pack.contactTracks.size()) ||
            !countFits(pack.clips.size())) {
            return false;
        }
        const auto uniqueStrings = [](const auto& values) {
            for (std::size_t index = 0u; index < values.size(); ++index) {
                if (values[index].empty() || !stringFits(values[index]) ||
                    std::find(values.begin(), values.begin() + index,
                        values[index]) != values.begin() + index) {
                    return false;
                }
            }
            return true;
        };
        if (!uniqueStrings(pack.jointNames)) {
            return false;
        }
        std::vector<std::string> contactIds;
        contactIds.reserve(pack.contactTracks.size());
        for (const InteractionContactTrack& track : pack.contactTracks) {
            if (track.id.empty() || track.taskContactGroup.empty() ||
                track.counterpart.empty() || !stringFits(track.id) ||
                !stringFits(track.taskContactGroup) ||
                !stringFits(track.counterpart) ||
                std::find(contactIds.begin(), contactIds.end(), track.id) !=
                    contactIds.end()) {
                return false;
            }
            contactIds.push_back(track.id);
        }
        const auto finiteValues = [](const std::span<const float> values) {
            return std::all_of(values.begin(), values.end(),
                [](const float value) { return std::isfinite(value); });
        };
        std::vector<std::string> clipIds;
        constexpr std::uint32_t validFlags =
            interactionSamplePredicted | interactionSamplePhysicsCertified;
        for (const InteractionClip& clip : pack.clips) {
            const std::uint64_t frames = clip.frameCount;
            const std::uint64_t samples = frames * pack.contactTracks.size();
            const std::uint64_t features =
                samples * kInteractionContactFeatureCount;
            if (clip.id.empty() || clip.desiredOutcome.empty() ||
                !stringFits(clip.id) || !stringFits(clip.desiredOutcome) ||
                std::find(clipIds.begin(), clipIds.end(), clip.id) !=
                    clipIds.end() || !std::isfinite(clip.framesPerSecond) ||
                !(clip.framesPerSecond > 0.0f) || clip.frameCount < 2u ||
                frames * kInteractionRootTargetCount != clip.rootTargets.size() ||
                frames * pack.jointNames.size() != clip.jointTargets.size() ||
                samples != clip.contactModes.size() ||
                samples != clip.contactFeatureMasks.size() ||
                samples != clip.contactSampleFlags.size() ||
                samples != clip.contactConfidence.size() ||
                features != clip.contactTargets.size() ||
                features != clip.contactTolerances.size() ||
                !finiteValues(clip.rootTargets) ||
                !finiteValues(clip.jointTargets) ||
                !finiteValues(clip.contactConfidence) ||
                !finiteValues(clip.contactTargets) ||
                !finiteValues(clip.contactTolerances)) {
                return false;
            }
            for (std::uint32_t frame = 0u; frame < clip.frameCount; ++frame) {
                const std::size_t root =
                    static_cast<std::size_t>(frame) * kInteractionRootTargetCount;
                const float x = clip.rootTargets[root + 3u];
                const float y = clip.rootTargets[root + 4u];
                const float z = clip.rootTargets[root + 5u];
                const float w = clip.rootTargets[root + 6u];
                if (std::abs(x * x + y * y + z * z + w * w - 1.0f) > 1.0e-3f) {
                    return false;
                }
            }
            for (std::size_t sample = 0u; sample < samples; ++sample) {
                if (clip.contactModes[sample] >
                        static_cast<std::uint32_t>(InteractionContactMode::release) ||
                    (clip.contactFeatureMasks[sample] &
                        ~kInteractionContactFeatureMask) != 0u ||
                    clip.contactSampleFlags[sample] == 0u ||
                    (clip.contactSampleFlags[sample] & ~validFlags) != 0u ||
                    clip.contactConfidence[sample] < 0.0f ||
                    clip.contactConfidence[sample] > 1.0f) {
                    return false;
                }
                const std::size_t base = sample * kInteractionContactFeatureCount;
                for (std::uint32_t feature = 0u;
                     feature < kInteractionContactFeatureCount;
                     ++feature) {
                    const bool valid =
                        (clip.contactFeatureMasks[sample] & (1u << feature)) != 0u;
                    const float target = clip.contactTargets[base + feature];
                    const float tolerance = clip.contactTolerances[base + feature];
                    if ((valid && !(tolerance > 0.0f)) ||
                        (!valid && tolerance < 0.0f) ||
                        (valid && feature >= 8u && target < 0.0f)) {
                        return false;
                    }
                }
            }
            clipIds.push_back(clip.id);
        }
        return true;
    } catch (...) {
        return false;
    }
}

} // namespace metalrobo
