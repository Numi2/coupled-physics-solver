#pragma once

#include "metalrobo/ArticulatedDynamics.hpp"

#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace metalrobo {

// Converts generalized and free-body state into one environment-major global
// body stream. Articulated bodies use the same FP64 kinematics as contact
// validation; this is solver state composition, not a rendering dependency.
[[nodiscard]] bool composeBodyStates(
    const EngineModel& model,
    std::uint32_t environmentCount,
    std::span<const float> q,
    std::span<const float> v,
    std::span<const MRBodyStateGPU> sceneBodies,
    std::vector<MRBodyStateGPU>& output,
    std::string* reason = nullptr
);

} // namespace metalrobo
