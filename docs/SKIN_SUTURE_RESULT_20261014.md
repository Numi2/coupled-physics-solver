# Native Metal skin stitch result — 2026-10-14

The prethreaded synthetic wound bench completes a 60 ms pull and a 20 ms hold.
Five elastic loops transfer 1.2 mm of reference material into their tails.
Mean top-lip gap changes from 0.600000 mm to 0.070125 mm
(88.313% reduction). All 80 native frames pass the fixed
certificates; the independent observer finds no intersecting incision faces.
The boundary prescribes reference material transfer; native elastic reaction
and its exact tangent move the tissue. This is an ideal-eyelet, tension-only
thread bench, with no needle insertion, free DER strand, knot, calibrated skin
material, or clinical qualification.

Run `./suture`. The separate bag remains in Numi2/numi-solver with `./bag`.
[Download the run evidence, including rejected attempts](https://github.com/Numi2/coupled-physics-solver/releases/tag/skin-suture-metal-20261014).

![Measured native gap, stitch tension and material transfer](https://github.com/Numi2/coupled-physics-solver/releases/download/skin-suture-metal-20261014/skin-suture-measured.png)

## What fixed the failures

- The old mixed nodal-pressure loaded case still rejects the eighth attempt.
  The working skin scene uses direct compressible neo-Hookean FEM with the
  same authored volumetric material energy. That is an explicit formulation
  change, not a qualification of the failed mixed path.
- A borrowed external-force field participates in Newton assembly and prevents
  the unloaded static-rest shortcut from hiding tiny loaded nodal impulses.
  Consistent momentum row scaling retains the existing 1e-4 equilibrium gate.
- Thread forces and the exact elastic/geometric tangent are evaluated on Metal
  inside the coupled implicit solve rather than lagged between frames.
- Small-feature contact mollifier and closest-point denominators now use
  fourth-power geometric floors. The scene retains all 15 vertex/triangle and
  edge/edge features per triangle pair, within its explicit capacity.
- Conservative CCD uses relative motion, terminates at a clear endpoint, and
  retains only a proven safe prefix after an unresolved grazing search.
  Contact feasibility protects the collision boundary, not barrier activation.
- Complete-manifold friction evaluates the regularized current-velocity force.
  Previously a stale low-speed history impulse was paired with a different
  derivative, stalling contact equilibrium. Authored friction remains 0.45 / 0.35.

Determinant stays above 0.20; mass drift remains below 1e-5; every one of the
2,304 tetrahedra stays active. No publication threshold was loosened. Optional
bounded solver retries journal their rejection and prove byte-exact rollback.
The successful 80-frame, 1 ms thread run required no retry.

## Evidence

The native tangent/CCD contract passes nine geometry cases. Thread tangent
relative finite-difference error is 0.000175719. Independent float64 checks
reconstruct gaps, elastic energy gradients, balanced reactions, material length,
and the actual incision triangle geometry. The 1 ms run's maximum exported
force reconstruction error is 3.03e-08 N; energy-gradient
error is 2.93e-12 N. Twenty-five fixed material markers
exit the loops. Minimum J is 0.982004344.

The one-entry-point `./suture check --headless` completes the native contract,
thread, no-thread and half-timestep arms. Only complete accepted trajectories
are used for the result. All rejected, cancelled and diagnostic attempts retain
their status in the release archive.

| Arm | Accepted frames | Final gap (mm) | Minimum J |
|---|---:|---:|---:|
| thread | 80 | 0.070124638 | 0.982004344 |
| no thread | 80 | 0.600000101 | 0.999999702 |
| half timestep | 160 | 0.069873003 | 0.981972039 |

Held endpoint timestep difference: **0.000251635 mm**,
below the preregistered 0.015 mm bound. The control stays open, with gap
drift below 1e-10 m. Maximum equilibrium residual is 9.99174154e-05.
The 1 ms thread trajectory takes approximately 81.7 GPU seconds
on the local Apple M4; this is an offline benchmark, not a real-time claim.

Owner regressions pass: FEM contact, heterogeneous FEM, device failure accounting,
solver budget overrides, and the bounded surgical tissue probe. The latter
retains its own exact replay and live FEM progress; this does not extend the
skin scene to a qualified needle procedure.

The 60-frame pull prefix was rerun exactly: tissue coordinates, wound-lip
coordinates, thread coordinates, thread state and thread forces match the
independent earlier successful run byte-for-byte after CSV parsing.

## Preserved history

The original local source snapshot is
[`b6dfd25c`](https://github.com/Numi2/numi-lab/tree/archive/skin-suture-local-20261014).
The [GitHub history release](https://github.com/Numi2/numi-lab/releases/tag/skin-suture-history-20261014)
contains the full old local and Mac mini directories, including source,
binaries, inputs, failed/incomplete results, file manifests, and separate Git
index/staged/unstaged recovery. Archives were restored and every file verified;
independent GitHub downloads matched SHA-256. The admitted old directories and
archive staging copies have been retired from disk.

The native AppKit/SceneKit bundle rendered all 80 measured frames and exported
its final image. Translucent skin and display tubes expose the buried loop;
these are viewing choices. Cancellation starts the native worker, stops only
that process group, records exit 130, and retains the successful latest replay.
No UI automation was used for these checks.

Qualified native executable SHA-256:
`cb3cdb5566b332bb093bdde5e29b36b9fc0a16a1c4accadd84028033fd6cbe31`.
Matter metallib SHA-256:
`472fa86938aa514e091ca45a945ded0019ff32daed2d59d2e5d907eb5dcbea7b`.

New source is in this repository. Native run receipts bind the source files,
input, executable and metallib hashes. Seed is zero and deterministic execution
is enabled. Frame CSVs, geometry, reactions, attempt journals, native replay,
regression logs and original rejected screens accompany the result release.
This remains a bounded software/mechanics benchmark; mesh convergence,
frictional channel damage, knot retention and physical calibration are unqualified.
