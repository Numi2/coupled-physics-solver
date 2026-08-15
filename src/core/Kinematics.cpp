#include "metalrobo/Kinematics.hpp"

#include <cmath>
#include <utility>

namespace metalrobo {
namespace {

bool fail(std::string* reason, std::string message) {
    if (reason != nullptr) {
        *reason = std::move(message);
    }
    return false;
}

bool finite(const mr_float4 value) {
    return std::isfinite(value.x) && std::isfinite(value.y) &&
        std::isfinite(value.z) && std::isfinite(value.w);
}

bool unitQuaternion(const mr_float4 value) {
    if (!finite(value)) {
        return false;
    }
    const double squared =
        static_cast<double>(value.x) * value.x +
        static_cast<double>(value.y) * value.y +
        static_cast<double>(value.z) * value.z +
        static_cast<double>(value.w) * value.w;
    return std::abs(squared - 1.0) <= 1.0e-4;
}

bool finiteValues(const std::span<const float> values) {
    for (const float value : values) {
        if (!std::isfinite(value)) {
            return false;
        }
    }
    return true;
}

} // namespace

bool composeBodyStates(
    const EngineModel& model,
    const std::uint32_t environmentCount,
    const std::span<const float> q,
    const std::span<const float> v,
    const std::span<const MRBodyStateGPU> sceneBodies,
    std::vector<MRBodyStateGPU>& output,
    std::string* reason
) {
    std::string modelReason;
    if (!model.valid(&modelReason) || environmentCount == 0u) {
        return fail(reason, "body composition requires a valid model and batch");
    }
    const std::size_t nq = model.world.nq;
    const std::size_t nv = model.world.nv;
    if (q.size() != environmentCount * nq ||
        v.size() != environmentCount * nv ||
        !finiteValues(q) || !finiteValues(v)) {
        return fail(reason, "body q/v streams have invalid dimensions or values");
    }

    std::vector<std::uint32_t> sceneIndices;
    sceneIndices.reserve(model.bodies.size());
    for (std::uint32_t body = 0u; body < model.bodies.size(); ++body) {
        if (model.bodies[body].articulationIndex == MR_INVALID_INDEX) {
            sceneIndices.push_back(body);
        }
    }
    if (sceneBodies.size() != environmentCount * sceneIndices.size()) {
        return fail(reason, "scene-body stream has invalid dimensions");
    }

    struct ArticulationScratch {
        std::vector<double> q;
        std::vector<double> v;
        std::vector<ArticulatedBodyKinematics> bodies;
    };
    std::vector<ArticulationScratch> scratch(model.articulations.size());
    for (std::size_t index = 0u; index < model.articulations.size(); ++index) {
        const MRArticulationGPU& articulation = model.articulations[index];
        scratch[index].q.resize(articulation.nq);
        scratch[index].v.resize(articulation.nv);
        scratch[index].bodies.resize(articulation.bodyCount);
    }

    std::vector<MRBodyStateGPU> candidate(
        static_cast<std::size_t>(environmentCount) * model.bodies.size());
    for (std::uint32_t environment = 0u;
         environment < environmentCount;
         ++environment) {
        for (std::uint32_t index = 0u; index < model.articulations.size(); ++index) {
            const MRArticulationGPU& articulation = model.articulations[index];
            ArticulationScratch& local = scratch[index];
            for (std::uint32_t coordinate = 0u; coordinate < articulation.nq; ++coordinate) {
                local.q[coordinate] =
                    q[environment * nq + articulation.qOffset + coordinate];
            }
            for (std::uint32_t coordinate = 0u; coordinate < articulation.nv; ++coordinate) {
                local.v[coordinate] =
                    v[environment * nv + articulation.vOffset + coordinate];
            }
            const ArticulatedDynamicsDiagnostics diagnostics =
                computeArticulatedBodyKinematics(
                    model, index, local.q, local.v, local.bodies);
            if (!diagnostics.succeeded()) {
                return fail(reason, "articulated body kinematics failed");
            }
            for (const ArticulatedBodyKinematics& body : local.bodies) {
                MRBodyStateGPU state{};
                state.position = {
                    static_cast<float>(body.centerOfMassPosition[0]),
                    static_cast<float>(body.centerOfMassPosition[1]),
                    static_cast<float>(body.centerOfMassPosition[2]), 1.0f};
                state.orientation = {
                    static_cast<float>(body.orientation[0]),
                    static_cast<float>(body.orientation[1]),
                    static_cast<float>(body.orientation[2]),
                    static_cast<float>(body.orientation[3])};
                state.linearVelocityAndInverseMass = {
                    static_cast<float>(body.linearVelocity[0]),
                    static_cast<float>(body.linearVelocity[1]),
                    static_cast<float>(body.linearVelocity[2]),
                    model.bodies[body.bodyIndex].massAndInverseMass.y};
                state.angularVelocity = {
                    static_cast<float>(body.angularVelocity[0]),
                    static_cast<float>(body.angularVelocity[1]),
                    static_cast<float>(body.angularVelocity[2]), 0.0f};
                state.flagsAndIndices[0] = model.bodies[body.bodyIndex].motionType;
                state.flagsAndIndices[1] =
                    model.bodies[body.bodyIndex].articulationIndex;
                state.flagsAndIndices[2] = body.bodyIndex;
                candidate[environment * model.bodies.size() + body.bodyIndex] = state;
            }
        }
        for (std::size_t local = 0u; local < sceneIndices.size(); ++local) {
            MRBodyStateGPU state = sceneBodies[environment * sceneIndices.size() + local];
            if (!finite(state.position) || !unitQuaternion(state.orientation) ||
                !finite(state.linearVelocityAndInverseMass) ||
                !finite(state.angularVelocity)) {
                return fail(reason, "scene-body state is nonfinite");
            }
            state.flagsAndIndices[2] = sceneIndices[local];
            candidate[environment * model.bodies.size() + sceneIndices[local]] = state;
        }
    }
    output = std::move(candidate);
    return true;
}

} // namespace metalrobo
