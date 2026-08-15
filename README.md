# Coupled Physics Solver

Apple-native coupled physics for articulated bodies, rigid contact,
heterogeneous Matter FEM/MPM, and discrete-elastic-rod sutures. The live step
runs through MetalWorld and lets Matter borrow the owning Metal command buffer;
there is no second simulation queue in the normal coupled path.

This is a focused extraction from Numi Lab. It contains the solver, material
compiler/runtime, surgical tissue and suture mechanics, and executable physics
probes. It does not contain training, MLX, rendering, app UI, or unrelated robot
products. See [PROVENANCE.md](PROVENANCE.md) for the exact source revision.

## Requirements

- Apple silicon Mac
- macOS with Xcode and the Metal 4 toolchain
- CMake 3.28 or newer

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The focused executables are written to `build/bin`:

- `coupled_matter_physics_probe`
- `coupled_surgical_tissue_probe`
- `coupled_suture_handoff_probe`

Promotion remains gated by exact replay, FP64 parity, physical outcomes, zero
failed steps, and same-device Metal timeline/counter evidence. A successful
build or probe alone is not a physical-fidelity claim.
