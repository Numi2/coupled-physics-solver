#!/usr/bin/env python3
"""Summarize displacement in the saved accepted suture-passage states."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from render_perfused_suture_pullthrough import read_state


state_dir = Path(__file__).resolve().parent / "states"
states = sorted(
    (read_state(path) for path in state_dir.glob("*.tsv")),
    key=lambda state: (state.step, state.phase),
)
if len(states) < 2:
    raise SystemExit("need at least two accepted states")

base = states[0]
needle_reference = next((state for state in states if state.needles), None)
base_needles = (
    np.asarray(needle_reference.needles)
    if needle_reference is not None
    else np.asarray([])
)
for state in states:
    tissue_delta = np.linalg.norm(state.nodes - base.nodes, axis=1)
    thread_delta = np.linalg.norm(state.thread - base.thread, axis=1)
    state_needles = np.asarray(state.needles)
    if base_needles.size and state_needles.shape == base_needles.shape:
        needle_proxy_max_mm = (
            np.linalg.norm(state_needles - base_needles, axis=2).max()
            * 1000.0
        )
        needle_metric = f"\tneedle_proxy_max_mm={needle_proxy_max_mm:.6f}"
    else:
        needle_metric = "\tneedle_proxy_max_mm=unavailable"
    print(
        f"{state.phase}\tstep={state.step}"
        f"\ttime_ms={state.simulation_time_s * 1000.0:.6f}"
        f"\tchannels={len(state.channels)}"
        f"\ttissue_max_um={tissue_delta.max() * 1.0e6:.3f}"
        f"\ttissue_rms_um={np.sqrt(np.mean(tissue_delta**2)) * 1.0e6:.3f}"
        f"\tthread_max_mm={thread_delta.max() * 1000.0:.6f}"
        f"{needle_metric}"
    )
if needle_reference is not None:
    print(f"needle_proxy_reference={needle_reference.phase}")
