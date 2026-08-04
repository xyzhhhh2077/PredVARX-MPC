"""Create a synchronized multi-panel GIF from the verified copyBJ SMPC figures."""

from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results_smpc"
DATA_PATH = RESULTS / "copyBJ_smpc_figure_data.npz"
GIF_PATH = RESULTS / "copyBJ_smpc_summary_200frames.gif"
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


def _style_axis(ax):
    ax.grid(True, color="#E5E7EB", linewidth=0.6, alpha=0.75)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=7)


def _box_edges(bound):
    corners = np.array(
        [
            [x, y, z]
            for x in (-bound, bound)
            for y in (-bound, bound)
            for z in (-bound, bound)
        ]
    )
    edges = []
    for i, a in enumerate(corners):
        for b in corners[i + 1 :]:
            if np.count_nonzero(a != b) == 1:
                edges.append((a, b))
    return edges


def main():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Run generate_smpc_figures.py first: {DATA_PATH}")

    data = np.load(DATA_PATH)
    t = data["time_seconds"]
    ref = data["reference_standardized"]
    task = data["smpc_task_standardized"]
    control = data["smpc_control"]
    innovation = data["composite_latent_innovation"]
    wind = data["physical_wind_velocity_mps"]
    cost = data["realized_stage_cost"]
    hard_bound = float(data["hard_bound_standardized"][0])
    input_bound = float(data["input_bound_standardized"][0])
    online_weight = float(data["online_weight"][0]) if "online_weight" in data.files else 0.8

    n = len(t)
    if not all(len(a) == n for a in (ref, task, control, innovation, wind, cost)):
        raise ValueError("All animation arrays must share the same sample count")
    if control.ndim != 2 or control.shape[1] != 4:
        raise ValueError("The Pelican control history must contain four input channels")
    if FRAME_COUNT != 200:
        raise ValueError("This deliverable is contracted to exactly 200 source frames")

    display_idx = np.unique(np.linspace(0, n - 1, min(DISPLAY_POINTS, n), dtype=int))
    frame_end = np.linspace(1, len(display_idx), FRAME_COUNT, dtype=int)
    td = t[display_idx]
    rd = ref[display_idx]
    yd = task[display_idx]
    ud = control[display_idx]
    ed = innovation[display_idx]
    wd = wind[display_idx]
    cd = np.maximum(cost[display_idx], 1e-5)

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial"],
            "axes.titlesize": 10,
            "axes.labelsize": 8,
            "legend.fontsize": 7,
        }
    )
    fig = plt.figure(figsize=(12.8, 11.5), dpi=90, facecolor="white")
    outer = fig.add_gridspec(
        4,
        2,
        height_ratios=(1.25, 0.62, 0.95, 0.72),
        hspace=0.38,
        wspace=0.19,
    )
    xyz_grid = outer[0, 0].subgridspec(3, 1, hspace=0.08)
    noise_grid = outer[2, :].subgridspec(2, 1, hspace=0.10)
    xyz_axes = [fig.add_subplot(xyz_grid[i, 0]) for i in range(3)]
    ax3d = fig.add_subplot(outer[0, 1], projection="3d")
    ax_control = fig.add_subplot(outer[1, :])
    ax_innov = fig.add_subplot(noise_grid[0, 0])
    ax_wind = fig.add_subplot(noise_grid[1, 0], sharex=ax_innov)
    ax_cost = fig.add_subplot(outer[3, :])

    fig.suptitle(
        f"copyBJ online SMPC (beta={online_weight:.1f}): trajectory, disturbances, stage cost",
        fontsize=14,
        fontweight="bold",
        y=0.985,
    )
    progress_text = fig.text(
        0.5, 0.952, "", ha="center", va="center", fontsize=9, color=BLACK,
    )

    xyz_lines = []
    xyz_markers = []
    xyz_reference_markers = []
    xyz_cursors = []
    labels = ("x", "y", "z")
    for j, (ax, label, color) in enumerate(zip(xyz_axes, labels, AXIS_COLORS)):
        ax.plot(td, rd[:, j], "--", color=BLACK, linewidth=1.8, alpha=1.0, zorder=6)
        ax.plot(td, yd[:, j], color=color, linewidth=0.9, alpha=0.12, zorder=2)
        line, = ax.plot([], [], color=color, linewidth=1.7, label=f"SMPC {label}", zorder=5)
        marker, = ax.plot([], [], "o", color=color, markersize=4.5, zorder=7)
        reference_marker, = ax.plot(
            [],
            [],
            marker="X",
            color=BLACK,
            markerfacecolor="white",
            markeredgewidth=1.1,
            markersize=5.5,
            zorder=8,
        )
        ax.axhline(hard_bound, color=RED, linestyle=":", linewidth=0.9)
        ax.axhline(-hard_bound, color=RED, linestyle=":", linewidth=0.9)
        cursor = ax.axvline(td[0], color=BLACK, linewidth=0.8, alpha=0.7)
        ax.set_ylim(-4.2, 4.2)
        ax.set_ylabel(label)
        _style_axis(ax)
        xyz_lines.append(line)
        xyz_markers.append(marker)
        xyz_reference_markers.append(reference_marker)
        xyz_cursors.append(cursor)
    xyz_axes[0].set_title("a  Task coordinates and references", loc="left", fontweight="bold")
    xyz_axes[0].plot([], [], "--", color=BLACK, linewidth=1.8, label="reference")
    xyz_axes[0].plot([], [], ":", color=RED, linewidth=1.0, label="hard bound ±3.8")
    xyz_axes[0].legend(loc="upper right", ncol=3, frameon=False)
    xyz_axes[0].tick_params(labelbottom=False)
    xyz_axes[1].tick_params(labelbottom=False)
    xyz_axes[2].set_xlabel("Time (s)")

    ax3d.plot(
        rd[:, 0],
        rd[:, 1],
        rd[:, 2],
        "--",
        color=ORANGE,
        linewidth=3.2,
        alpha=1.0,
        marker="X",
        markevery=120,
        markersize=4.5,
        markeredgecolor=BLACK,
        markeredgewidth=0.6,
        label="reference",
    )
    ax3d.plot(yd[:, 0], yd[:, 1], yd[:, 2], color=BLUE, linewidth=0.8, alpha=0.10)
    line3d, = ax3d.plot([], [], [], color=BLUE, linewidth=1.8, label="online SMPC")
    marker3d, = ax3d.plot([], [], [], "o", color=RED, markersize=5)
    reference_marker3d, = ax3d.plot(
        [], [], [], marker="X", color=ORANGE, markeredgecolor=BLACK, markersize=6,
    )
    for a, b in _box_edges(hard_bound):
        ax3d.plot(*zip(a, b), color=RED, linewidth=0.65, alpha=0.32)
    ax3d.set_title("b  Three-dimensional trajectory", loc="left", fontweight="bold")
    ax3d.set_xlabel("x", labelpad=1)
    ax3d.set_ylabel("y", labelpad=1)
    ax3d.set_zlabel("z", labelpad=1)
    ax3d.set_xlim(-4.15, 4.15)
    ax3d.set_ylim(-4.15, 4.15)
    ax3d.set_zlim(-4.15, 4.15)
    ax3d.set_box_aspect((1, 1, 1))
    ax3d.view_init(elev=24, azim=-53)
    ax3d.legend(loc="upper right", frameon=False)

    control_lines = []
    for j, color in enumerate(CONTROL_COLORS):
        ax_control.plot(td, ud[:, j], color=color, linewidth=0.55, alpha=0.10)
        line, = ax_control.plot([], [], color=color, linewidth=1.05, label=f"u{j + 1}")
        control_lines.append(line)
    ax_control.axhline(
        input_bound,
        color=RED,
        linestyle=":",
        linewidth=1.1,
        label=f"input bounds ±{input_bound:.0f}",
    )
    ax_control.axhline(-input_bound, color=RED, linestyle=":", linewidth=1.1)
    control_cursor = ax_control.axvline(td[0], color=BLACK, linewidth=0.8)
    ax_control.set_xlim(td[0], td[-1])
    ax_control.set_ylim(-1.12 * input_bound, 1.12 * input_bound)
    ax_control.set_title(
        "c  Standardized control inputs and box constraints",
        loc="left",
        fontweight="bold",
    )
    ax_control.set_ylabel("Control u")
    ax_control.set_xlabel("Time (s)")
    ax_control.legend(loc="upper right", ncol=5, frameon=False)
    _style_axis(ax_control)

    innov_lines = []
    for j, color in enumerate(INNOVATION_COLORS):
        ax_innov.plot(td, ed[:, j], color=color, linewidth=0.45, alpha=0.10)
        line, = ax_innov.plot([], [], color=color, linewidth=0.85, label=f"ε{j + 1}")
        innov_lines.append(line)
    innov_cursor = ax_innov.axvline(td[0], color=BLACK, linewidth=0.8)
    ax_innov.set_title(
        "d  Composite latent innovation and injected physical wind",
        loc="left",
        fontweight="bold",
    )
    ax_innov.set_ylabel("Innovation")
    ax_innov.legend(loc="upper right", ncol=5, frameon=False)
    ax_innov.tick_params(labelbottom=False)
    _style_axis(ax_innov)

    wind_lines = []
    for j, (label, color) in enumerate(zip(labels, AXIS_COLORS)):
        ax_wind.plot(td, wd[:, j], color=color, linewidth=0.45, alpha=0.10)
        line, = ax_wind.plot([], [], color=color, linewidth=0.85, label=f"w{label}")
        wind_lines.append(line)
    wind_cursor = ax_wind.axvline(td[0], color=BLACK, linewidth=0.8)
    ax_wind.set_ylabel("Wind (m/s)")
    ax_wind.set_xlabel("Time (s)")
    ax_wind.legend(loc="upper right", ncol=3, frameon=False)
    _style_axis(ax_wind)

    ax_cost.plot(td, cd, color=RED, linewidth=0.7, alpha=0.12)
    cost_line, = ax_cost.plot([], [], color=RED, linewidth=1.35)
    cost_marker, = ax_cost.plot([], [], "o", color=BLACK, markersize=4)
    cost_cursor = ax_cost.axvline(td[0], color=BLACK, linewidth=0.8)
    ax_cost.set_yscale("log")
    ax_cost.set_xlim(td[0], td[-1])
    ax_cost.set_ylim(max(1e-4, cd.min() * 0.65), cd.max() * 1.45)
    ax_cost.set_title("e  Realized stage cost (log scale)", loc="left", fontweight="bold")
    ax_cost.set_xlabel("Time (s)")
    ax_cost.set_ylabel("Stage cost ℓk")
    ax_cost.text(
        0.03,
        0.95,
        "ℓk = Q‖sk − rk‖² + R‖uk − ū‖²\n(not the horizon QP objective)",
        transform=ax_cost.transAxes,
        va="top",
        fontsize=7.5,
        color=BLACK,
        bbox={"facecolor": "white", "edgecolor": LIGHT_GRAY, "alpha": 0.90, "pad": 4},
    )
    _style_axis(ax_cost)

    def update(frame):
        end = int(frame_end[frame])
        sl = slice(0, end)
        current = end - 1
        now = td[current]
        progress_text.set_text(
            f"t = {now:6.2f} s   |   frame {frame + 1:03d}/{FRAME_COUNT}   |   "
            f"hard bound ±{hard_bound:.1f}   |   beta={online_weight:.1f}"
        )
        for j in range(3):
            xyz_lines[j].set_data(td[sl], yd[sl, j])
            xyz_markers[j].set_data([now], [yd[current, j]])
            xyz_reference_markers[j].set_data([now], [rd[current, j]])
            xyz_cursors[j].set_xdata([now, now])
        line3d.set_data_3d(yd[sl, 0], yd[sl, 1], yd[sl, 2])
        marker3d.set_data_3d([yd[current, 0]], [yd[current, 1]], [yd[current, 2]])
        reference_marker3d.set_data_3d(
            [rd[current, 0]], [rd[current, 1]], [rd[current, 2]]
        )
        for j, line in enumerate(control_lines):
            line.set_data(td[sl], ud[sl, j])
        control_cursor.set_xdata([now, now])
        for j, line in enumerate(innov_lines):
            line.set_data(td[sl], ed[sl, j])
        innov_cursor.set_xdata([now, now])
        for j, line in enumerate(wind_lines):
            line.set_data(td[sl], wd[sl, j])
        wind_cursor.set_xdata([now, now])
        cost_line.set_data(td[sl], cd[sl])
        cost_marker.set_data([now], [cd[current]])
        cost_cursor.set_xdata([now, now])
        return []

    ani = animation.FuncAnimation(
        fig, update, frames=FRAME_COUNT, interval=FRAME_DURATION_MS, blit=False,
    )
    RESULTS.mkdir(parents=True, exist_ok=True)
    ani.save(GIF_PATH, writer=animation.PillowWriter(fps=1000 / FRAME_DURATION_MS), dpi=90)
    plt.close(fig)
    print(GIF_PATH)


if __name__ == "__main__":
    main()
