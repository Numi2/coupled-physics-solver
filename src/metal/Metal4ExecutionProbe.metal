#include <metal_stdlib>

using namespace metal;

kernel void coupled_metal4_contract_probe(
    device const uint* input [[buffer(0)]],
    device uint* output [[buffer(1)]],
    const uint index [[thread_position_in_grid]]
) {
    output[index] = input[index] ^ (0x9e3779b9u + index * 0x85ebca6bu);
}
