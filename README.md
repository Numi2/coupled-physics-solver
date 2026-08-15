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

The first Metal 4 ownership path requires macOS 26 and Xcode 26. The existing
Matter/MetalWorld solver still uses the original Metal encoder API while its
dispatch graph is migrated; the Metal 4 contract probe is not evidence that
the complete solver has moved to that path.

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
- `coupled_trial_check`
- `coupled_perfused_active_tissue_contract_probe`
- `coupled_perfused_active_tissue_probe`
- `coupled_metal4_execution_probe`

`coupled_trial_check MANIFEST.json` validates the versioned manifest, stream
shapes, normalized relative paths, byte counts, and SHA-256 payload identity.
Perfused-active material authoring has no numerical defaults: every geometry,
layer, perfusion, activation, damage, and contact scalar must carry bounded
evidence provenance before it can be configured.

The workspace command below runs exact replay and captures Metal evidence. GPU
counters are strict by default; use `--timeline-only` only when deliberately
recording timeline evidence without making counter claims.

```sh
numi coupled-profile --mode metal4
numi coupled-profile --mode tissue
numi coupled-profile --mode perfused
numi coupled-profile --mode coupled-perfused
numi coupled-profile --mode heterogeneous
numi coupled-profile --mode heterogeneous-mixed
```

Heterogeneous FEM tetrahedra may own distinct passive constitutive materials
and densities while the object material owns external contact. Mixed-field
heterogeneous elements additionally own their transport, activation, fibre,
and active-stress coefficients. They must share the nodal pressure scale (bulk
modulus and thermal expansion). The compiler still rejects heterogeneous
mutable topology until every mutation transaction preserves material identity.
The perfused probe executes the four-layer path on Metal with synthetic,
provenance-shaped fixtures. It also gates every state-free layer's compiled
passive and mixed pressure/pore-pressure/active-fibre stress and mechanical
tangent against an independent FP64 scalar evaluator, with a separate FP64
directional-difference check. This is a constitutive-point parity boundary,
not a full-step FP64 oracle. The `coupled-perfused` mode connects that same
heterogeneous, active-field tissue directly to live DER strand contact through
the hard needle swage. It deliberately disables puncture mutation until
mutation can preserve per-tetrahedron material identity. Both modes prove
execution and replay, not ex-vivo calibration or indistinguishable
visual/physical fidelity.

Promotion remains gated by exact replay, FP64 parity, physical outcomes, zero
failed steps, and same-device Metal timeline/counter evidence. A successful
build or probe alone is not a physical-fidelity claim.
