#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "coupled/Metal4Execution.hpp"

#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef COUPLED_METALLIB
#define COUPLED_METALLIB ""
#endif

namespace {

void require(const bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

} // namespace

int main() {
    try {
        @autoreleasepool {
            coupled::Metal4ExecutionContext context;
            auto diagnostics = context.initialize(3u);
            require(diagnostics.accepted && context.valid(), diagnostics.message);
            id<MTLDevice> device = (__bridge id<MTLDevice>)context.device();

            constexpr std::uint32_t kCount = 32u;
            id<MTLBuffer> input = [device
                newBufferWithLength:kCount * sizeof(std::uint32_t)
                options:MTLResourceStorageModeShared];
            id<MTLBuffer> output = [device
                newBufferWithLength:kCount * sizeof(std::uint32_t)
                options:MTLResourceStorageModeShared];
            require(input != nil && output != nil, "failed to allocate probe buffers");
            auto* inputValues = static_cast<std::uint32_t*>(input.contents);
            auto* outputValues = static_cast<std::uint32_t*>(output.contents);
            for (std::uint32_t index = 0u; index < kCount; ++index) {
                inputValues[index] = index * 17u + 3u;
                outputValues[index] = 0u;
            }
            diagnostics = context.addPersistentAllocation((__bridge void*)input);
            require(diagnostics.accepted, diagnostics.message);
            diagnostics = context.addPersistentAllocation((__bridge void*)output);
            require(diagnostics.accepted && diagnostics.residencyBytes > 0u,
                    diagnostics.message);

            NSError* error = nil;
            NSURL* libraryURL = [NSURL fileURLWithPath:
                [NSString stringWithUTF8String:COUPLED_METALLIB]];
            id<MTLLibrary> library = [device newLibraryWithURL:libraryURL error:&error];
            require(library != nil, error == nil
                ? "failed to load coupled metallib"
                : error.localizedDescription.UTF8String);
            id<MTLFunction> function = [library
                newFunctionWithName:@"coupled_metal4_contract_probe"];
            require(function != nil, "Metal 4 contract kernel is missing");
            id<MTLComputePipelineState> pipeline = [device
                newComputePipelineStateWithFunction:function
                error:&error];
            require(pipeline != nil, error == nil
                ? "failed to create Metal 4 probe pipeline"
                : error.localizedDescription.UTF8String);

            MTL4ArgumentTableDescriptor* argumentDescriptor =
                [[MTL4ArgumentTableDescriptor alloc] init];
            argumentDescriptor.label = @"coupled-metal4-probe-arguments";
            argumentDescriptor.maxBufferBindCount = 2u;
            argumentDescriptor.initializeBindings = YES;
            id<MTL4ArgumentTable> arguments = [device
                newArgumentTableWithDescriptor:argumentDescriptor
                error:&error];
            require(arguments != nil, error == nil
                ? "failed to create Metal 4 argument table"
                : error.localizedDescription.UTF8String);
            [arguments setAddress:input.gpuAddress atIndex:0u];
            [arguments setAddress:output.gpuAddress atIndex:1u];

            coupled::Metal4Submission submission;
            diagnostics = context.beginSubmission(
                "coupled-metal4-contract",
                submission
            );
            require(diagnostics.accepted, diagnostics.message);
            id<MTL4CommandBuffer> commandBuffer =
                (__bridge id<MTL4CommandBuffer>)submission.commandBuffer;
            id<MTL4ComputeCommandEncoder> encoder =
                [commandBuffer computeCommandEncoder];
            require(encoder != nil, "failed to create Metal 4 compute encoder");
            encoder.label = @"coupled-metal4-contract-dispatch";
            [encoder setComputePipelineState:pipeline];
            [encoder setArgumentTable:arguments];
            [encoder dispatchThreads:MTLSizeMake(kCount, 1u, 1u)
                threadsPerThreadgroup:MTLSizeMake(kCount, 1u, 1u)];
            [encoder endEncoding];
            diagnostics = context.endAndCommit(submission);
            require(diagnostics.accepted, diagnostics.message);
            diagnostics = context.wait(submission, 30000u);
            require(diagnostics.accepted, diagnostics.message);

            for (std::uint32_t index = 0u; index < kCount; ++index) {
                const std::uint32_t expected = inputValues[index] ^
                    (0x9e3779b9u + index * 0x85ebca6bu);
                require(outputValues[index] == expected,
                        "Metal 4 argument-table output mismatch");
            }
            require(
                diagnostics.allocatorCount == 3u &&
                    diagnostics.allocatorBytes > 0u &&
                    diagnostics.residencyBytes > 0u,
                "Metal 4 ownership diagnostics are incomplete"
            );
            std::cout << "Metal 4 execution contract passed: device="
                      << diagnostics.deviceName
                      << " allocators=" << diagnostics.allocatorCount
                      << " allocator_bytes=" << diagnostics.allocatorBytes
                      << " residency_bytes=" << diagnostics.residencyBytes
                      << '\n';
        }
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "Metal 4 execution contract failed: "
                  << exception.what() << '\n';
        return 1;
    }
}
