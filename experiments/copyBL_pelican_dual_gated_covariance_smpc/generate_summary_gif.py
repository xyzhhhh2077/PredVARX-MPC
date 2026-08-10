"""Create a BJ-layout synchronized GIF with copyBL dual-gate evidence."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
MEDIA = HERE / "results_media"
DATA_PATH = RESULTS / "copyBL_dual_gated_run.npz"
EVENT_PATH = RESULTS / "copyBL_dual_gate_events.json"
SUMMARY_PATH = RESULTS / "copyBL_dual_gated_summary.json"
GIF_PATH = MEDIA / "copyBL_dual_gate_summary_200frames.gif"
FRAME_COUNT = 200
FRAME_DURATION_MS = 90
DISPLAY_POINTS = 1800

BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
RED = "#D55E00"
PURPLE = "#CC79A7"
SKY = "#56B4E9"
BLACK = "#111827"
GRAY = "#6B7280"
LIGHT_GRAY = "#D1D5DB"
AXIS_COLORS = (BLUE, ORANGE, GREEN)
INNOVATION_COLORS = (BLUE, ORANGE, GREEN, PURPLE, SKY)
CONTROL_COLORS = (BLUE, ORANGE, GREEN, PURPLE)


def _style_axis(axis):
    axis.grid(True, color="#E5E7EB", linewidth=0.6, alpha=0.75)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(labelsize=7)


def _box_edges(bound):
    corners = np.array(
        [[x, y, z] for x in (-bound, bound) for y in (-bound, bound) for z in (-bound, bound)]
    )
    return [
        (a, b)
        for index, a in enumerate(corners)
        for b in corners[index + 1:]
        if np.count_nonzero(a != b) == 1
    ]


def main():
    data = np.load(DATA_PATH)
    events = json.loads(EVENT_PATH.read_text(encoding="utf-8"))
    summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
    task = data["trajectory_standardized"]
    reference = data["reference_standardized"]
    control = data["control"]
    innovation = data["residual_innovation_latent"]
    wind = data["wind_velocity_mps"]
    cost = np.maximum(data["stage_cost"], 1e-5)
    dt = float(summary["sample_time_seconds"])
    time_s = np.arange(len(task), dtype=float) * dt
    hard_bound = float(summary["hard_bound_standardized"])
    input_bound = 7.0
    if not all(len(array) == len(time_s) for array in (reference, task, control, innovation, wind, cost)):
        raise ValueError("All animation arrays must share the same sample count")
    if control.ndim != 2 or control.shape[1] != 4:
        raise ValueError("The Pelican control history must contain four channels")

    display_idx = np.unique(
        np.linspace(0, len(time_s) - 1, min(DISPLAY_POINTS, len(time_s)), dtype=int)
    )
    frame_end = np.linspace(1, len(display_idx), FRAME_COUNT, dtype=int)
    td = time_s[display_idx]
    rd = reference[display_idx]
    yd = task[display_idx]
    ud = control[display_idx]
    ed = innovation[display_idx]
    wd = wind[display_idx]
    cd = cost[display_idx]
    event_steps = np.asarray([event["step"] for event in events], dtype=int)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial"],
        "axes.titlesize": 10,
        "axes.labelsize": 8,
        "legend.fontsize": 7,
    })
    figure = plt.figure(figsize=(12.8, 11.5), dpi=90, facecolor="white")
    outer = figure.add_gridspec(
        4, 2, height_ratios=(1.25, 0.62, 0.95, 0.72), hspace=0.38, wspace=0.19,
    )
    xyz_grid = outer[0, 0].subgridspec(3, 1, hspace=0.08)
    noise_grid = outer[2, :].subgridspec(2, 1, hspace=0.10)
    xyz_axes = [figure.add_subplot(xyz_grid[index, 0]) for index in range(3)]
    axis_3d = figure.add_subplot(outer[0, 1], projection="3d")
    axis_control = figure.add_subplot(outer[1, :])
    axis_innovation = figure.add_subplot(noise_grid[0, 0])
    axis_wind = figure.add_subplot(noise_grid[1, 0], sharex=axis_innovation)
    axis_cost = figure.add_subplot(outer[3, :])

    figure.suptitle(
        "copyBL dual-gated covariance SMPC: trajectory, disturbances, stage cost",
        fontsize=14, fontweight="bold", y=0.985,
    )
    progress_text = figure.text(
        0.5, 0.952, "", ha="center", va="center", fontsize=9, color=BLACK,
    )

    xyz_lines, xyz_markers, reference_markers, xyz_cursors = [], [], [], []
    for index, (axis, label, color) in enumerate(zip(xyz_axes, "xyz", AXIS_COLORS)):
        axis.plot(td, rd[:, index], "--", color=BLACK, linewidth=1.8, zorder=6)
        axis.plot(td, yd[:, index], color=color, linewidth=0.9, alpha=0.12, zorder=2)
        line, = axis.plot([], [], color=color, linewidth=1.7, label=f"SMPC {label}", zorder=5)
        marker, = axis.plot([], [], "o", color=color, markersize=4.5, zorder=7)
        reference_marker, = axis.plot(
            [], [], marker="X", color=BLACK, markerfacecolor="white",
            markeredgewidth=1.1, markersize=5.5, zorder=8,
        )
        axis.axhline(hard_bound, color=RED, linestyle=":", linewidth=0.9)
        axis.axhline(-hard_bound, color=RED, linestyle=":", linewidth=0.9)
        cursor = axis.axvline(td[0], color=BLACK, linewidth=0.8, alpha=0.7)
        axis.set_ylim(-4.2, 4.2)
        axis.set_ylabel(label)
        _style_axis(axis)
        xyz_lines.append(line)
        xyz_markers.append(marker)
        reference_markers.append(reference_marker)
        xyz_cursors.append(cursor)
    xyz_axes[0].set_title("a  Task coordinates and references", loc="left", fontweight="bold")
    xyz_axes[0].plot([], [], "--", color=BLACK, linewidth=1.8, label="reference")
    xyz_axes[0].plot([], [], ":", color=RED, linewidth=1.0, label="hard bound +/-3.8")
    xyz_axes[0].legend(loc="upper right", ncol=3, frameon=False)
    xyz_axes[0].tick_params(labelbottom=False)
    xyz_axes[1].tick_params(labelbottom=False)
    xyz_axes[2].set_xlabel("Time (s)")

    axis_3d.plot(
        rd[:, 0], rd[:, 1], rd[:, 2], "--", color=ORANGE, linewidth=3.2,
        marker="X", markevery=120, markersize=4.5, markeredgecolor=BLACK,
        markeredgewidth=0.6, label="reference",
    )
    axis_3d.plot(yd[:, 0], yd[:, 1], yd[:, 2], color=BLUE, linewidth=0.8, alpha=0.10)
    line_3d, = axis_3d.plot([], [], [], color=BLUE, linewidth=1.8, label="dual-gated SMPC")
    marker_3d, = axis_3d.plot([], [], [], "o", color=RED, markersize=5)
    reference_marker_3d, = axis_3d.plot(
        [], [], [], marker="X", color=ORANGE, markeredgecolor=BLACK, markersize=6,
    )
    for start, end in _box_edges(hard_bound):
        axis_3d.plot(*zip(start, end), color=RED, linewidth=0.65, alpha=0.32)
    axis_3d.set_title("b  Three-dimensional trajectory", loc="left", fontweight="bold")
    axis_3d.set_xlabel("x", labelpad=1)
    axis_3d.set_ylabel("y", labelpad=1)
    axis_3d.set_zlabel("z", labelpad=1)
    axis_3d.set_xlim(-4.15, 4.15)
    axis_3d.set_ylim(-4.15, 4.15)
    axis_3d.set_zlim(-4.15, 4.15)
    axis_3d.set_box_aspect((1, 1, 1))
    axis_3d.view_init(elev=24, azim=-53)
    axis_3d.legend(loc="upper right", frameon=False)

    control_lines = []
    for index, color in enumerate(CONTROL_COLORS):
        axis_control.plot(td, ud[:, index], color=color, linewidth=0.55, alpha=0.10)
        line, = axis_control.plot([], [], color=color, linewidth=1.05, label=f"u{index + 1}")
        control_lines.append(line)
    axis_control.axhline(input_bound, color=RED, linestyle=":", linewidth=1.1, label="input bounds +/-7")
    axis_control.axhline(-input_bound, color=RED, linestyle=":", linewidth=1.1)
    control_cursor = axis_control.axvline(td[0], color=BLACK, linewidth=0.8)
    axis_control.set_xlim(td[0], td[-1])
    axis_control.set_ylim(-1.12 * input_bound, 1.12 * input_bound)
    axis_control.set_title("c  Standardized control inputs and box constraints", loc="left", fontweight="bold")
    axis_control.set_ylabel("Control u")
    axis_control.set_xlabel("Time (s)")
    axis_control.legend(loc="upper right", ncol=5, frameon=False)
    _style_axis(axis_control)

    innovation_lines = []
    for index, color in enumerate(INNOVATION_COLORS):
        axis_innovation.plot(td, ed[:, index], color=color, linewidth=0.45, alpha=0.10)
        line, = axis_innovation.plot([], [], color=color, linewidth=0.85, label=f"eps{index + 1}")
        innovation_lines.append(line)
    innovation_cursor = axis_innovation.axvline(td[0], color=BLACK, linewidth=0.8)
    axis_innovation.set_title(
        "d  Composite latent innovation and injected physical wind", loc="left", fontweight="bold",
    )
    axis_innovation.set_ylabel("Innovation")
    axis_innovation.legend(loc="upper right", ncol=5, frameon=False)
    axis_innovation.tick_params(labelbottom=False)
    _style_axis(axis_innovation)

    wind_lines = []
    for index, (label, color) in enumerate(zip("xyz", AXIS_COLORS)):
        axis_wind.plot(td, wd[:, index], color=color, linewidth=0.45, alpha=0.10)
        line, = axis_wind.plot([], [], color=color, linewidth=0.85, label=f"w{label}")
        wind_lines.append(line)
    wind_cursor = axis_wind.axvline(td[0], color=BLACK, linewidth=0.8)
    axis_wind.set_ylabel("Wind (m/s)")
    axis_wind.set_xlabel("Time (s)")
    axis_wind.legend(loc="upper right", ncol=3, frameon=False)
    _style_axis(axis_wind)

    axis_cost.plot(td, cd, color=RED, linewidth=0.7, alpha=0.12)
    cost_line, = axis_cost.plot([], [], color=RED, linewidth=1.35)
    cost_marker, = axis_cost.plot([], [], "o", color=BLACK, markersize=4)
    cost_cursor = axis_cost.axvline(td[0], color=BLACK, linewidth=0.8)
    gate_cursor = axis_cost.axvline(td[0], color=GRAY, linewidth=1.25, alpha=0.0)
    axis_cost.set_yscale("log")
    axis_cost.set_xlim(td[0], td[-1])
    axis_cost.set_ylim(max(1e-4, cd.min() * 0.65), cd.max() * 1.45)
    axis_cost.set_title("e  Realized stage cost and latest dual-gate decision", loc="left", fontweight="bold")
    axis_cost.set_xlabel("Time (s)")
    axis_cost.set_ylabel("Stage cost")
    axis_cost.text(
        0.03, 0.95, "stage cost = Q||s-r||^2 + R||u-u_mean||^2\n"
        "green: commit | purple: value reject | red: feasibility reject",
        transform=axis_cost.transAxes, va="top", fontsize=7.5, color=BLACK,
        bbox={"facecolor": "white", "edgecolor": LIGHT_GRAY, "alpha": 0.90, "pad": 4},
    )
    _style_axis(axis_cost)

    def update(frame):
        end = int(frame_end[frame])
        shown = slice(0, end)
        current = end - 1
        now = td[current]
        source_step = int(display_idx[current])
        prior = np.flatnonzero(event_steps <= source_step)
        if prior.size:
            event = events[int(prior[-1])]
            if event["update_committed"]:
                gate_label, gate_color = "COMMITTED", GREEN
            elif event["value_triggered"]:
                gate_label, gate_color = "VALUE REJECTED", PURPLE
            elif not event["feasibility_accepted"]:
                gate_label, gate_color = "FEASIBILITY REJECTED", RED
            else:
                gate_label, gate_color = "FEASIBILITY PASSED / VALUE NOT SCHEDULED", GRAY
            gate_cursor.set_xdata([event["time_seconds"], event["time_seconds"]])
            gate_cursor.set_color(gate_color)
            gate_cursor.set_alpha(0.95)
        else:
            gate_label = "COLD START"
        progress_text.set_text(
            f"t = {now:6.2f} s   |   frame {frame + 1:03d}/{FRAME_COUNT}   |   "
            f"hard bound +/-{hard_bound:.1f}   |   latest gate: {gate_label}"
        )
        for index in range(3):
            xyz_lines[index].set_data(td[shown], yd[shown, index])
            xyz_markers[index].set_data([now], [yd[current, index]])
            reference_markers[index].set_data([now], [rd[current, index]])
            xyz_cursors[index].set_xdata([now, now])
        line_3d.set_data_3d(yd[shown, 0], yd[shown, 1], yd[shown, 2])
        marker_3d.set_data_3d([yd[current, 0]], [yd[current, 1]], [yd[current, 2]])
        reference_marker_3d.set_data_3d(
            [rd[current, 0]], [rd[current, 1]], [rd[current, 2]],
        )
        for index, line in enumerate(control_lines):
            line.set_data(td[shown], ud[shown, index])
        control_cursor.set_xdata([now, now])
        for index, line in enumerate(innovation_lines):
            line.set_data(td[shown], ed[shown, index])
        innovation_cursor.set_xdata([now, now])
        for index, line in enumerate(wind_lines):
            line.set_data(td[shown], wd[shown, index])
        wind_cursor.set_xdata([now, now])
        cost_line.set_data(td[shown], cd[shown])
        cost_marker.set_data([now], [cd[current]])
        cost_cursor.set_xdata([now, now])
        return []

    animation_object = animation.FuncAnimation(
        figure, update, frames=FRAME_COUNT, interval=FRAME_DURATION_MS, blit=False,
    )
    MEDIA.mkdir(parents=True, exist_ok=True)
    animation_object.save(
        GIF_PATH, writer=animation.PillowWriter(fps=1000 / FRAME_DURATION_MS), dpi=90,
    )
    plt.close(figure)
    print(GIF_PATH)


if __name__ == "__main__":
    main()
