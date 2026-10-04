# Perfused suture entry prefix — 2026-10-04

[![Puncture-accepted state from the two-frame clip](perfused-suture-entry-prefix-20261004-poster.png)](perfused-suture-entry-prefix-20261004.mp4)

[Watch the two-frame clip](perfused-suture-entry-prefix-20261004.mp4). It shows the saved pre-entry state followed by the first accepted puncture state from the native Metal probe. The clip holds only those exact snapshots at 1 fps; it does not interpolate or exaggerate motion.

## Recorded result

- The synthetic four-layer coupon published 13,068 tissue-node positions, 8,384 boundary triangles, and 128 discrete-elastic-rod thread nodes in both snapshots.
- At accepted step 1 (0.03125 ms), the saved Matter state contains one active puncture channel. Comparing the two saved tissue meshes gives a maximum node displacement of 5.783 µm and an RMS displacement of 0.0641 µm.
- The curved-passage stage submitted its first 32-step device batch. It was still waiting for the Metal command buffer after 5 min 50 s, with no post-entry accepted snapshot returned. The run was then stopped at the time cap. This is a censored engineering attempt, not a solver rejection or a completed needle passage.
- On the same build, the focused `coupled.suture.perfused_tissue_puncture` and `coupled.suture.perfused_tissue_coupling` CTest gates passed in 13.72 s and 20.99 s. These isolate channel admission and live tissue/thread coupling; they do not qualify the unfinished curved-passage batch.
- The run did not establish thread transport through the wall, distal clearance, wound closure, physical calibration, or clinical capability. The fixture parameters are synthetic.

The native probe was built and run on an Apple M4 Pro Mac mini running macOS 26.6, from base revision `1ddeb40d78be1d7d7eb262adf001adea2376faf4` plus [`simulation-worktree.patch`](simulation-worktree.patch). The binary and changed-source hashes are in [`fingerprints.txt`](fingerprints.txt); [`manifest.txt`](manifest.txt) records the device and build context. The saved native states, run output, and renderer are retained beside this record so the clip can be regenerated.

## Reproduction

From the patched checkout, build `coupled_suture_handoff_probe`, then run:

```sh
build/bin/coupled_suture_handoff_probe \
  --perfused-tissue-curved-pull-through-only \
  --state-output-dir /tmp/perfused-suture-entry-states
```

The clip was rendered from the two captured TSV states with `render_perfused_suture_pullthrough.py` at 1 fps. The exact video and poster hashes appear in the repository media checksum list and the Numi Lab website's `public/media/SHA256SUMS`.
