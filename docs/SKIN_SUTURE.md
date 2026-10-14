# Skin stitch bench

Run `./suture` from this repository on an Apple silicon Mac with Xcode's
Metal 4 toolchain and CMake. It builds the scoped native target, pulls five
prethreaded stitches, checks the exported states, and opens an AppKit/SceneKit
replay. Drag to inspect; scrub or pause to inspect a particular accepted frame.
Toggle “See thread inside” to inspect the buried loop through translucent skin.
Display tube radii and opacity are viewing choices, not additional physics.

```
./suture                       # simulate and view
./suture run --headless         # simulate and check without a window
./suture view                   # open the last successful replay
./suture check --headless       # thread, no-thread control, half-timestep check
```

A rejected nonlinear/linear frame may retry once with 24 Newton and 256 Krylov
iterations, only after a byte-exact accepted-state rollback check. Tolerances
remain fixed; both attempts are recorded in `attempts.csv`.

Every run gets a unique directory under `build/suture/runs/`. Rejected and
interrupted runs retain their accepted states and logs. `build/suture/latest`
points only to the last successful, independently checked replay. Ctrl-C stops
the owned native process. Closing the replay window ends the viewer.

The bag has its own implementation and entry point in
[Numi2/numi-solver](https://github.com/Numi2/numi-solver). It shares no scene,
launch script, evidence directory, or material with this skin bench.

## What executes

The specimen is a 30 × 24 × 1.5 mm homogeneous synthetic coupon, with a
16 mm incision initially open by 0.600 mm. The 24 × 8 × 2 mesh has 2,304
active tetrahedra and fixed longitudinal ends. Its compressible neo-Hookean
material uses authored density 1,050 kg/m³, shear modulus 25 kPa and lambda
250 kPa. These are scene values, not a measured skin fit.

Five elastic, tension-only threads pass through two ideal, frictionless
material-attached eyelets per stitch, 1 mm behind the incision lips. The
eyelet paths start prethreaded. The condensed thread model has EA = 0.5 N.
A 60 ms cosine tail pull transfers 1.2 mm of material per stitch from the
loop into the tail; total reference thread length is conserved. Fixed material
markers move through the eyelets. The tail boundary holds the authored pull for a further 20 ms. Reference
material flux is prescribed by this boundary condition; it is not a predicted
insertion or abrasion rate.

Thread reaction is the negative gradient of its elastic energy. Its exact
geometric/material tangent is assembled into the same native Newton/FGMRES
operator as FEM. Every candidate recomputes the live thread reaction on Metal;
there is no lagged external spring and no prescribed skin position. Native
swept contact and positive-separation gates retain lip nonpenetration.

This bench demonstrates material pull-through and mechanical approximation of
a wound. It does not execute needle insertion, tissue channel creation,
frictional abrasion, a free DER strand, knot formation, a robot procedure,
healing, or a calibrated clinical result. The earlier needle/robot work remains
separate from this prethreaded closure bench.

## Loaded failure and verification

The original mixed-pressure case rejected its eighth attempted frame. It was
reproduced here, including with a larger Krylov basis. A nodal mixed pressure
approximation did not satisfy the existing per-element volume constraint under
this load. The working bench uses the material's direct compressible
neo-Hookean formulation; it retains volumetric stiffness instead of loosening
that failed mixed constraint. The mixed formulation remains an archived failed
experiment, not a newly qualified solver path.

FEM momentum rows are explicitly scaled, with the same scale in the residual,
operator and preconditioner. This prevents tiny nodal impulses from passing a
unit-sized residual floor without a tissue response. The default scale in
other scenes is one. Equilibrium tolerance remains 1e-4, determinant must exceed
0.20, mass drift must be below 1e-5, and all 2,304 elements must remain active.

The contact search uses relative displacement bounds, terminates after a clear
endpoint, and returns a conservative safe-prefix bound for an unresolved
grazing event. It never treats exhausted search as a proven clear sweep.
Candidate barriers and final positive-separation certificates still govern
acceptance. The contact mollifier uses a fourth-power length floor; the old
squared-length floor incorrectly erased contact on small incision triangles.
Closest-point denominators also use fourth-power floors. The feasibility bound
protects the collision boundary rather than freezing motion at barrier activation.
This scene uses continuous regularized Coulomb friction: retaining an old tangent
impulse at low speed while assembling the current velocity derivative stalled
closure. Its authored skin friction remains 0.45 static / 0.35 dynamic.

This scene retains all 15 vertex/triangle and edge/edge features per triangle
pair, with an explicit bounded capacity. Other scenes retain their existing
manifold setting.

`tools/check_skin_suture.py` independently reconstructs gaps and thread forces
from exported coordinates and interpolation weights using float64 arithmetic.
It checks the full incision surfaces for triangle intersection, as well as
energy gradients, balanced reaction, conserved material length,
material markers exiting the loop, native certificates and mass. The closure
bench requires at least 70% mean gap reduction, a final mean gap below 0.180 mm,
and actual material passage. These are software bench criteria, not clinical
acceptance criteria. A no-thread ablation must retain the open wound; halving
the timestep must change the endpoint gap by less than 0.015 mm.

## History

The old local source was preserved in the immutable Git tag
[`archive/skin-suture-local-20261014`](https://github.com/Numi2/numi-lab/tree/archive/skin-suture-local-20261014).
The complete local directory and the Mac mini's separate source, binaries,
inputs, profiles and incomplete results are preserved in the
[GitHub history release](https://github.com/Numi2/numi-lab/releases/tag/skin-suture-history-20261014),
with file manifests and independently restored archive receipts. The old
failed traction and incomplete passage evidence retain their original status.

Current source and measured run evidence are recorded in this repository's
[skin suture result](SKIN_SUTURE_RESULT_20261014.md).
