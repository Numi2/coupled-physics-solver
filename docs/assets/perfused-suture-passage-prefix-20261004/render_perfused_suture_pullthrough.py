#!/usr/bin/env python3
"""Render accepted native states from the perfused suture pull-through probe."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation
from matplotlib.collections import PolyCollection
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


COLORS = ("#8ed0d4", "#efbd78", "#e77f91", "#b69ae5")
LAYER_NAMES = ("Layer 1", "Layer 2", "Layer 3", "Layer 4")
NAVY = "#0a1220"
PANEL = "#111d30"
TEXT = "#edf4ff"
MUTED = "#9db0c8"
THREAD = "#48c5f0"
NEEDLE = "#f5d36c"


@dataclass
class State:
    path: Path
    phase: str
    step: int
    simulation_time_s: float
    timestep_s: float
    nodes: np.ndarray
    triangles: np.ndarray
    materials: np.ndarray
    thread: np.ndarray
    channels: list[tuple[np.ndarray, np.ndarray]]
    needles: list[tuple[np.ndarray, np.ndarray]]
    longitudinal: np.ndarray
    thickness: np.ndarray


def read_state(path: Path) -> State:
    values: dict[str, str] = {}
    nodes: list[list[float]] = []
    triangles: list[list[int]] = []
    materials: list[int] = []
    thread: list[list[float]] = []
    channels: list[tuple[np.ndarray, np.ndarray]] = []
    needles: list[tuple[np.ndarray, np.ndarray]] = []
    with path.open(newline="") as stream:
        for row in csv.reader(stream, delimiter="\t"):
            if not row:
                continue
            if row[0] == "schema":
                if len(row) < 2 or row[1] != "numi.perfused-suture-entry-state.v2":
                    raise ValueError(f"unsupported state schema in {path}")
            elif row[0] in {"phase", "step", "simulation_time_s", "timestep_s", "longitudinal_axis", "thickness_axis"}:
                values[row[0]] = "\t".join(row[1:])
            elif row[0] == "tissue_node":
                nodes.append([float(x) for x in row[2:5]])
            elif row[0] == "tissue_surface_triangle":
                triangles.append([int(x) for x in row[2:5]])
                materials.append(int(row[5]))
            elif row[0] == "thread_node":
                thread.append([float(x) for x in row[2:5]])
            elif row[0] == "puncture_channel":
                origin = np.asarray([float(x) for x in row[2:5]])
                radius = float(row[5])
                axis = np.asarray([float(x) for x in row[6:9]])
                half_length = float(row[9])
                axis_norm = np.linalg.norm(axis)
                if axis_norm > 0.0:
                    axis = axis / axis_norm
                    channels.append((origin - axis * half_length, origin + axis * half_length))
            elif row[0] == "rigid_proxy":
                flags = int(row[2])
                if flags & (1 << 2):
                    start = np.asarray([float(x) for x in row[4:7]])
                    end = np.asarray([float(x) for x in row[8:11]])
                    needles.append((start, end))

    if not nodes or not triangles or not thread:
        raise ValueError(f"incomplete accepted geometry state: {path}")
    longitudinal = np.asarray([float(x) for x in values["longitudinal_axis"].split("\t")])
    thickness = np.asarray([float(x) for x in values["thickness_axis"].split("\t")])
    longitudinal /= np.linalg.norm(longitudinal)
    thickness /= np.linalg.norm(thickness)
    return State(
        path=path,
        phase=values["phase"],
        step=int(values["step"]),
        simulation_time_s=float(values["simulation_time_s"]),
        timestep_s=float(values["timestep_s"]),
        nodes=np.asarray(nodes, dtype=float),
        triangles=np.asarray(triangles, dtype=np.int64),
        materials=np.asarray(materials, dtype=np.int64),
        thread=np.asarray(thread, dtype=float),
        channels=channels,
        needles=needles,
        longitudinal=longitudinal,
        thickness=thickness,
    )


def project(points: np.ndarray, along: np.ndarray, through: np.ndarray) -> np.ndarray:
    return np.column_stack((points @ along, points @ through)) * 1000.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("state_directory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--poster", type=Path)
    parser.add_argument("--fps", type=int, default=2)
    args = parser.parse_args()
    if args.fps < 1 or args.fps > 12:
        parser.error("--fps must be between 1 and 12")

    states = sorted(
        (read_state(path) for path in args.state_directory.glob("*.tsv")),
        key=lambda state: (state.step, state.phase),
    )
    if len(states) < 2:
        raise SystemExit(f"need at least two accepted states; found {len(states)}")
    for state in states[1:]:
        if not np.array_equal(state.triangles, states[0].triangles):
            raise SystemExit("tissue surface topology changed between saved states")
        if not np.array_equal(state.materials, states[0].materials):
            raise SystemExit("tissue material ownership changed between saved states")

    all_points = np.vstack([
        point
        for state in states
        for point in (state.nodes, state.thread,
                      *(np.vstack(pair) for pair in state.needles))
    ]) * 1000.0
    low = all_points.min(axis=0)
    high = all_points.max(axis=0)
    span = np.maximum(high - low, 1.0)
    low -= span * 0.08
    high += span * 0.08

    first = states[0]
    fig = plt.figure(figsize=(12.8, 7.2), facecolor=NAVY)
    grid = fig.add_gridspec(1, 2, width_ratios=(1.35, 1.0), left=0.055,
                            right=0.985, top=0.86, bottom=0.18, wspace=0.13)
    ax3d = fig.add_subplot(grid[0, 0], projection="3d", facecolor=PANEL)
    ax2d = fig.add_subplot(grid[0, 1], facecolor=PANEL)
    for ax in (ax3d, ax2d):
        ax.tick_params(colors=MUTED, labelsize=8)
        for spine in ax.spines.values():
            spine.set_color("#34445d")

    triangle_colors = [COLORS[int(i) % len(COLORS)] for i in first.materials]
    surface3d = Poly3DCollection(
        first.nodes[first.triangles] * 1000.0,
        facecolors=triangle_colors,
        edgecolors=(0.05, 0.08, 0.13, 0.22),
        linewidths=0.12,
        alpha=0.69,
    )
    ax3d.add_collection3d(surface3d)
    ax3d.set_xlim(low[0], high[0])
    ax3d.set_ylim(low[1], high[1])
    ax3d.set_zlim(low[2], high[2])
    ax3d.set_box_aspect(np.maximum(high - low, 1.0))
    ax3d.view_init(elev=22, azim=-54)
    ax3d.set_xlabel("x (mm)", color=MUTED, labelpad=4)
    ax3d.set_ylabel("y (mm)", color=MUTED, labelpad=4)
    ax3d.set_zlabel("z (mm)", color=MUTED, labelpad=4)
    ax3d.set_title("ACCEPTED 3D TISSUE STATE", color=TEXT, fontsize=10, pad=10,
                   loc="left", fontweight="bold")
    for axis in (ax3d.xaxis, ax3d.yaxis, ax3d.zaxis):
        axis.pane.set_facecolor((*to_rgb(PANEL), 1.0))
        axis.pane.set_edgecolor("#34445d")

    thread3d, = ax3d.plot([], [], [], color=THREAD, linewidth=2.3,
                          solid_capstyle="round", label="DER thread")
    max_needle_lines = max(len(state.needles) for state in states)
    needle3d = [ax3d.plot([], [], [], color=NEEDLE, linewidth=2.5,
                          solid_capstyle="round")[0]
                for _ in range(max_needle_lines)]
    max_channel_lines = max(len(state.channels) for state in states)
    channel3d = [ax3d.plot([], [], [], color="#ff8c72", linewidth=1.5,
                           linestyle="--", alpha=0.8)[0]
                 for _ in range(max_channel_lines)]

    projected_triangles = project(first.nodes[first.triangles].reshape(-1, 3),
                                  first.longitudinal, first.thickness).reshape(-1, 3, 2)
    surface2d = PolyCollection(
        projected_triangles,
        facecolors=triangle_colors,
        edgecolors=(0.05, 0.08, 0.13, 0.28),
        linewidths=0.13,
        alpha=0.78,
    )
    ax2d.add_collection(surface2d)
    projected_all = np.vstack([
        project(np.vstack((state.nodes, state.thread)), state.longitudinal,
                state.thickness)
        for state in states
    ])
    lo2, hi2 = projected_all.min(axis=0), projected_all.max(axis=0)
    span2 = np.maximum(hi2 - lo2, 1.0)
    ax2d.set_xlim(lo2[0] - span2[0] * 0.08, hi2[0] + span2[0] * 0.08)
    ax2d.set_ylim(lo2[1] - span2[1] * 0.14, hi2[1] + span2[1] * 0.14)
    ax2d.set_aspect("equal", adjustable="box")
    ax2d.set_xlabel("along specimen (mm)", color=MUTED)
    ax2d.set_ylabel("through thickness (mm)", color=MUTED)
    ax2d.set_title("THROUGH-THICKNESS VIEW", color=TEXT, fontsize=10, pad=10,
                   loc="left", fontweight="bold")
    thread2d, = ax2d.plot([], [], color=THREAD, linewidth=2.2,
                          solid_capstyle="round", label="DER thread")
    needle2d = [ax2d.plot([], [], color=NEEDLE, linewidth=2.4,
                          solid_capstyle="round")[0]
                for _ in range(max_needle_lines)]
    channel2d = [ax2d.plot([], [], color="#ff8c72", linewidth=1.6,
                           linestyle="--", alpha=0.8)[0]
                 for _ in range(max_channel_lines)]
    phase_text = fig.text(0.057, 0.91, "", color="#8ed0d4", fontsize=9,
                          fontweight="bold")
    step_text = fig.text(0.985, 0.91, "", color=MUTED, fontsize=9,
                         ha="right")
    legend_handles = [
        Patch(facecolor=COLORS[index], edgecolor="none", label=LAYER_NAMES[index])
        for index in range(len(COLORS))
    ] + [
        Line2D([0], [0], color=THREAD, linewidth=2.3, label="DER thread"),
        Line2D([0], [0], color=NEEDLE, linewidth=2.5, label="Needle proxies"),
        Line2D([0], [0], color="#ff8c72", linewidth=1.6, linestyle="--",
               label="Puncture channels"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=7,
               bbox_to_anchor=(0.52, 0.09), frameon=False,
               labelcolor=MUTED, fontsize=8, handlelength=1.8,
               columnspacing=1.4)
    last_phase = states[-1].phase
    if last_phase == "pull-through-certified":
        footer = (
            "SYNTHETIC FOUR-LAYER PERFUSED TISSUE  ·  SAVED ACCEPTED STATES ONLY  ·  NOT CLINICAL EVIDENCE"
        )
    elif last_phase.startswith("passage-step-"):
        footer = (
            "CURVED NEEDLE PASSAGE PREFIX ONLY  ·  THROUGH-WALL CLEARANCE AND THREAD PULL-THROUGH NOT SHOWN  ·  NOT CLINICAL EVIDENCE"
        )
    elif last_phase.startswith("pull-"):
        footer = (
            "PULL-THROUGH PREFIX ONLY  ·  COMPLETION GATES NOT SHOWN  ·  NOT CLINICAL EVIDENCE"
        )
    else:
        footer = (
            "ENTRY PREFIX ONLY  ·  NEEDLE PASSAGE AND THREAD PULL-THROUGH INCOMPLETE  ·  NOT CLINICAL EVIDENCE"
        )
    fig.text(0.057, 0.035,
             footer,
             color=MUTED, fontsize=8, family="sans-serif")

    def update(index: int):
        state = states[index]
        surface3d.set_verts(state.nodes[state.triangles] * 1000.0)
        surface2d.set_verts(project(
            state.nodes[state.triangles].reshape(-1, 3),
            state.longitudinal,
            state.thickness,
        ).reshape(-1, 3, 2))
        thread3d.set_data_3d(*((state.thread * 1000.0).T))
        thread2d.set_data(*project(state.thread, state.longitudinal,
                                   state.thickness).T)
        for artists, curves, scale in (
            (needle3d, state.needles, 1000.0),
            (channel3d, state.channels, 1000.0),
        ):
            for artist, (start, end) in zip(artists, curves):
                points = np.vstack((start, end)) * scale
                artist.set_data_3d(points[:, 0], points[:, 1], points[:, 2])
            for artist in artists[len(curves):]:
                artist.set_data_3d([], [], [])
        for artists, curves in ((needle2d, state.needles),
                                (channel2d, state.channels)):
            for artist, (start, end) in zip(artists, curves):
                points = project(np.vstack((start, end)),
                                 state.longitudinal, state.thickness)
                artist.set_data(points[:, 0], points[:, 1])
            for artist in artists[len(curves):]:
                artist.set_data([], [])
        phase_text.set_text(state.phase.replace("-", " ").upper())
        step_text.set_text(
            f"ACCEPTED STEP {state.step}  ·  {state.simulation_time_s * 1000.0:.2f} ms"
        )
        return (surface3d, surface2d, thread3d, thread2d, *needle3d,
                *needle2d, *channel3d, *channel2d, phase_text, step_text)

    poster_index = len(states) - 1
    update(poster_index)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.poster:
        args.poster.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.poster, dpi=150, facecolor=fig.get_facecolor())
    animation = FuncAnimation(fig, update, frames=len(states), interval=1000 // args.fps,
                              blit=False, repeat=False)
    animation.save(args.output, writer=FFMpegWriter(
        fps=args.fps, codec="libx264", bitrate=6500,
        extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"],
        metadata={"title": (
                      "Synthetic perfused tissue curved needle passage prefix"
                      if last_phase.startswith("passage-step-")
                      else "Synthetic perfused tissue suture entry attempt"
                  ),
                  "comment": "Accepted native states; no frame interpolation"},
    ), dpi=100)
    plt.close(fig)
    print(f"states={len(states)}")
    print(f"fps={args.fps}")
    print(f"video={args.output}")
    if args.poster:
        print(f"poster={args.poster}")


if __name__ == "__main__":
    main()
