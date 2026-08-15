#pragma once

#include <cstdint>
#include <memory>
#include <string>

namespace coupled {

struct Metal4ExecutionDiagnostics {
    bool accepted = false;
    std::string deviceName;
    std::string message;
    std::uint32_t allocatorCount = 0u;
    std::uint64_t allocatorBytes = 0u;
    std::uint64_t residencyBytes = 0u;
};

// Opaque ticket for one borrowed Metal 4 command buffer. The context owns the
// command buffer and allocator; callers may encode into commandBuffer but may
// not end, commit, wait, retain, or replace it.
struct Metal4Submission {
    void* commandBuffer = nullptr; // id<MTL4CommandBuffer>
    std::uint32_t allocatorSlot = 0u;
    std::uint64_t completionValue = 0u;
};

// Single Metal 4 submission authority for coupled physics. It owns one queue,
// a bounded allocator ring, queue-level residency and completion feedback.
// This is intentionally generic Metal 4: it performs no GPU-family selection.
class Metal4ExecutionContext {
public:
    Metal4ExecutionContext();
    ~Metal4ExecutionContext();
    Metal4ExecutionContext(Metal4ExecutionContext&&) noexcept;
    Metal4ExecutionContext& operator=(Metal4ExecutionContext&&) noexcept;
    Metal4ExecutionContext(const Metal4ExecutionContext&) = delete;
    Metal4ExecutionContext& operator=(const Metal4ExecutionContext&) = delete;

    [[nodiscard]] Metal4ExecutionDiagnostics initialize(
        std::uint32_t allocatorCount = 3u
    );
    [[nodiscard]] bool valid() const noexcept;
    [[nodiscard]] void* device() const noexcept; // id<MTLDevice>
    [[nodiscard]] void* queue() const noexcept; // id<MTL4CommandQueue>

    // Residency may change only when no ticket is encoding or in flight.
    [[nodiscard]] Metal4ExecutionDiagnostics addPersistentAllocation(
        void* allocation // id<MTLAllocation>
    );

    [[nodiscard]] Metal4ExecutionDiagnostics beginSubmission(
        std::string label,
        Metal4Submission& submission
    );
    [[nodiscard]] Metal4ExecutionDiagnostics endAndCommit(
        Metal4Submission& submission
    );
    [[nodiscard]] Metal4ExecutionDiagnostics wait(
        Metal4Submission& submission,
        std::uint64_t timeoutMilliseconds
    );
    [[nodiscard]] Metal4ExecutionDiagnostics diagnostics() const;

private:
    struct State;
    std::unique_ptr<State> state_;
};

} // namespace coupled
