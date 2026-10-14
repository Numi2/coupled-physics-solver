#pragma once
#include "shared.h"
// A condensed, tension-only elastic thread sliding through authored eyelets.
// Corner order: lower top, upper top, upper bottom, lower bottom.
typedef struct NM_ALIGN16 NMSutureAnchorGPU {
    nm_uint4 identity; // FEM node, loop, corner, reserved
    nm_float4 weights; // interpolation weight, reserved
} NMSutureAnchorGPU;
typedef struct NM_ALIGN16 NMSutureLoopGPU {
    nm_float4 material; // initial loop rest length, total rest length, EA, maximum pull
    nm_float4 state; // current loop rest length, tail rest length, material pulled, tension
    nm_float4 diagnostics; // geometric loop length, elastic energy, stretch, reserved
} NMSutureLoopGPU;
