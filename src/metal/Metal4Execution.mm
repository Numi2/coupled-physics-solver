#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "coupled/Metal4Execution.hpp"

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <limits>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

namespace coupled {
namespace {

Metal4ExecutionDiagnostics reject(std::string message) {
    return {.accepted = false, .message = std::move(message)};
}

} // namespace

struct Metal4ExecutionContext::State {
    struct CommitFeedbackState {
        std::mutex mutex;
        std::condition_variable condition;
        bool received = false;
        std::string error;
        double gpuStartSeconds = 0.0;
        double gpuEndSeconds = 0.0;
    };

    struct Slot {
        id<MTL4CommandAllocator> allocator = nil;
        id<MTL4CommandBuffer> commandBuffer = nil;
        std::uint64_t completionValue = 0u;
        bool encoding = false;
        bool inFlight = false;
        std::shared_ptr<CommitFeedbackState> feedback;
    };

    id<MTLDevice> device = nil;
    id<MTL4CommandQueue> queue = nil;
    id<MTLResidencySet> residency = nil;
    id<MTLSharedEvent> completionEvent = nil;
    std::vector<Slot> slots;
    std::uint32_t nextSlot = 0u;
    std::uint64_t nextCompletionValue = 1u;
    mutable std::mutex mutex;

    ~State() {
        for (auto& slot : slots) {
            if (slot.inFlight && completionEvent != nil) {
                [completionEvent
                    waitUntilSignaledValue:slot.completionValue
                    timeoutMS:30000u];
            }
        }
        if (residency != nil) [residency endResidency];
    }

    [[nodiscard]] bool busy() const noexcept {
        return std::ranges::any_of(slots, [](const Slot& slot) {
            return slot.encoding || slot.inFlight;
        });
    }

    [[nodiscard]] Metal4ExecutionDiagnostics snapshot(
        const bool accepted,
        std::string message
    ) const {
        Metal4ExecutionDiagnostics result;
        result.accepted = accepted;
        result.message = std::move(message);
        result.deviceName = device == nil || device.name.UTF8String == nullptr
            ? ""
            : device.name.UTF8String;
        result.allocatorCount = static_cast<std::uint32_t>(slots.size());
        for (const auto& slot : slots) {
            if (slot.allocator != nil) {
                result.allocatorBytes += slot.allocator.allocatedSize;
            }
        }
        result.residencyBytes = residency == nil ? 0u : residency.allocatedSize;
        return result;
    }
};

Metal4ExecutionContext::Metal4ExecutionContext() = default;

Metal4ExecutionContext::~Metal4ExecutionContext() = default;

Metal4ExecutionContext::Metal4ExecutionContext(
    Metal4ExecutionContext&&
) noexcept = default;

Metal4ExecutionContext& Metal4ExecutionContext::operator=(
    Metal4ExecutionContext&&
) noexcept = default;

Metal4ExecutionDiagnostics Metal4ExecutionContext::initialize(
    const std::uint32_t allocatorCount
) {
    if (allocatorCount == 0u || allocatorCount > 8u) {
        return reject("Metal 4 allocator count must be in [1, 8]");
    }
    if (state_) {
        return reject("Metal 4 execution context is already initialized");
    }
    auto staged = std::make_unique<State>();
    staged->device = MTLCreateSystemDefaultDevice();
    if (staged->device == nil) {
        return reject("no Apple Metal device is available");
    }
    staged->queue = [staged->device newMTL4CommandQueue];
    if (staged->queue == nil) {
        return reject("device does not provide a Metal 4 command queue");
    }
    staged->completionEvent = [staged->device newSharedEvent];
    if (staged->completionEvent == nil) {
        return reject("failed to create Metal 4 completion event");
    }
    staged->completionEvent.label = @"coupled-metal4-completion";

    MTLResidencySetDescriptor* residencyDescriptor =
        [[MTLResidencySetDescriptor alloc] init];
    residencyDescriptor.label = @"coupled-persistent-residency";
    residencyDescriptor.initialCapacity = 64u;
    NSError* residencyError = nil;
    staged->residency = [staged->device
        newResidencySetWithDescriptor:residencyDescriptor
        error:&residencyError];
    if (staged->residency == nil) {
        const char* message = residencyError == nil
            ? nullptr
            : residencyError.localizedDescription.UTF8String;
        return reject(message == nullptr
            ? "failed to create Metal residency set"
            : message);
    }
    [staged->residency commit];
    [staged->residency requestResidency];
    [staged->queue addResidencySet:staged->residency];

    staged->slots.resize(allocatorCount);
    for (std::uint32_t index = 0u; index < allocatorCount; ++index) {
        MTL4CommandAllocatorDescriptor* descriptor =
            [[MTL4CommandAllocatorDescriptor alloc] init];
        descriptor.label = [NSString stringWithFormat:
            @"coupled-command-allocator-%u", index];
        NSError* allocatorError = nil;
        staged->slots[index].allocator = [staged->device
            newCommandAllocatorWithDescriptor:descriptor
            error:&allocatorError];
        if (staged->slots[index].allocator == nil) {
            const char* message = allocatorError == nil
                ? nullptr
                : allocatorError.localizedDescription.UTF8String;
            return reject(message == nullptr
                ? "failed to create Metal 4 command allocator"
                : message);
        }
    }
    state_ = std::move(staged);
    return state_->snapshot(true, "Metal 4 execution context initialized");
}

bool Metal4ExecutionContext::valid() const noexcept {
    return state_ && state_->device != nil && state_->queue != nil &&
        state_->residency != nil && state_->completionEvent != nil &&
        !state_->slots.empty();
}

void* Metal4ExecutionContext::device() const noexcept {
    return valid() ? (__bridge void*)state_->device : nullptr;
}

void* Metal4ExecutionContext::queue() const noexcept {
    return valid() ? (__bridge void*)state_->queue : nullptr;
}

Metal4ExecutionDiagnostics Metal4ExecutionContext::addPersistentAllocation(
    void* allocation
) {
    if (!valid() || allocation == nullptr) {
        return reject("Metal 4 residency requires a valid context and allocation");
    }
    std::scoped_lock lock(state_->mutex);
    if (state_->busy()) {
        return state_->snapshot(
            false,
            "Metal 4 residency cannot change while a submission is active"
        );
    }
    id<MTLAllocation> resource = (__bridge id<MTLAllocation>)allocation;
    if ([state_->residency containsAllocation:resource]) {
        return state_->snapshot(true, "allocation is already resident");
    }
    [state_->residency addAllocation:resource];
    [state_->residency commit];
    [state_->residency requestResidency];
    return state_->snapshot(true, "persistent allocation added to queue residency");
}

Metal4ExecutionDiagnostics Metal4ExecutionContext::beginSubmission(
    std::string label,
    Metal4Submission& submission
) {
    if (!valid() || label.empty()) {
        return reject("Metal 4 submission requires a valid context and label");
    }
    std::scoped_lock lock(state_->mutex);
    for (std::uint32_t offset = 0u; offset < state_->slots.size(); ++offset) {
        const std::uint32_t index =
            (state_->nextSlot + offset) % state_->slots.size();
        auto& slot = state_->slots[index];
        if (slot.encoding || slot.inFlight) continue;
        slot.commandBuffer = [state_->device newCommandBuffer];
        if (slot.commandBuffer == nil) {
            return state_->snapshot(false, "failed to create Metal 4 command buffer");
        }
        [slot.commandBuffer beginCommandBufferWithAllocator:slot.allocator];
        slot.commandBuffer.label = [NSString stringWithUTF8String:label.c_str()];
        slot.completionValue = state_->nextCompletionValue++;
        slot.encoding = true;
        slot.feedback = std::make_shared<State::CommitFeedbackState>();
        state_->nextSlot = (index + 1u) % state_->slots.size();
        submission = {
            .commandBuffer = (__bridge void*)slot.commandBuffer,
            .allocatorSlot = index,
            .completionValue = slot.completionValue,
        };
        return state_->snapshot(true, "Metal 4 command buffer opened");
    }
    return state_->snapshot(false, "all Metal 4 allocator slots are in flight");
}

Metal4ExecutionDiagnostics Metal4ExecutionContext::endAndCommit(
    Metal4Submission& submission
) {
    if (!valid() || submission.commandBuffer == nullptr ||
        submission.allocatorSlot >= state_->slots.size()) {
        return reject("Metal 4 commit ticket is invalid");
    }
    std::scoped_lock lock(state_->mutex);
    auto& slot = state_->slots[submission.allocatorSlot];
    if (!slot.encoding || slot.inFlight ||
        slot.completionValue != submission.completionValue ||
        (__bridge void*)slot.commandBuffer != submission.commandBuffer) {
        return state_->snapshot(false, "Metal 4 commit ticket is stale");
    }
    [slot.commandBuffer endCommandBuffer];
    slot.encoding = false;
    slot.inFlight = true;

    MTL4CommitOptions* options = [[MTL4CommitOptions alloc] init];
    const auto feedbackState = slot.feedback;
    [options addFeedbackHandler:^(id<MTL4CommitFeedback> feedback) {
        std::scoped_lock feedbackLock(feedbackState->mutex);
        feedbackState->received = true;
        feedbackState->gpuStartSeconds = feedback.GPUStartTime;
        feedbackState->gpuEndSeconds = feedback.GPUEndTime;
        if (feedback.error != nil) {
            const char* text = feedback.error.localizedDescription.UTF8String;
            feedbackState->error = text == nullptr
                ? "Metal 4 commit failed"
                : text;
        }
        feedbackState->condition.notify_all();
    }];
    id<MTL4CommandBuffer> buffers[] = {slot.commandBuffer};
    [state_->queue commit:buffers count:1u options:options];
    [state_->queue
        signalEvent:state_->completionEvent
        value:submission.completionValue];
    return state_->snapshot(true, "Metal 4 command buffer committed");
}

Metal4ExecutionDiagnostics Metal4ExecutionContext::wait(
    Metal4Submission& submission,
    const std::uint64_t timeoutMilliseconds
) {
    if (!valid() || submission.commandBuffer == nullptr ||
        submission.allocatorSlot >= state_->slots.size() ||
        timeoutMilliseconds == 0u) {
        return reject("Metal 4 wait ticket or timeout is invalid");
    }
    const auto timeout = std::chrono::milliseconds(timeoutMilliseconds);
    const auto deadline = std::chrono::steady_clock::now() + timeout;
    if (![state_->completionEvent
            waitUntilSignaledValue:submission.completionValue
            timeoutMS:timeoutMilliseconds]) {
        std::scoped_lock lock(state_->mutex);
        return state_->snapshot(false, "Metal 4 submission timed out");
    }

    std::shared_ptr<State::CommitFeedbackState> feedbackState;
    {
        std::scoped_lock lock(state_->mutex);
        auto& slot = state_->slots[submission.allocatorSlot];
        if (!slot.inFlight || slot.completionValue != submission.completionValue ||
            (__bridge void*)slot.commandBuffer != submission.commandBuffer) {
            return state_->snapshot(false, "Metal 4 wait ticket is stale");
        }
        feedbackState = slot.feedback;
    }
    std::unique_lock feedbackLock(feedbackState->mutex);
    if (!feedbackState->condition.wait_until(
            feedbackLock,
            deadline,
            [&feedbackState] { return feedbackState->received; })) {
        feedbackLock.unlock();
        std::scoped_lock lock(state_->mutex);
        return state_->snapshot(false, "Metal 4 commit feedback timed out");
    }
    const std::string feedbackError = feedbackState->error;
    feedbackLock.unlock();
    std::scoped_lock lock(state_->mutex);
    auto& slot = state_->slots[submission.allocatorSlot];
    if (!slot.inFlight || slot.completionValue != submission.completionValue ||
        (__bridge void*)slot.commandBuffer != submission.commandBuffer) {
        return state_->snapshot(false, "Metal 4 wait ticket became stale");
    }
    [slot.allocator reset];
    slot.commandBuffer = nil;
    slot.encoding = false;
    slot.inFlight = false;
    slot.feedback.reset();
    submission = {};
    return state_->snapshot(
        feedbackError.empty(),
        feedbackError.empty()
            ? "Metal 4 submission completed"
            : feedbackError
    );
}

Metal4ExecutionDiagnostics Metal4ExecutionContext::diagnostics() const {
    if (!valid()) return reject("Metal 4 execution context is not initialized");
    std::scoped_lock lock(state_->mutex);
    return state_->snapshot(true, "Metal 4 execution context is valid");
}

} // namespace coupled
