"""Plot the three copyBF boundary scenarios in three-dimensional task space."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def _draw_bound_box(ax, limits, color, linestyle, linewidth, label):
    lx, ly, lz = np.asarray(limits, dtype=float)
    corners = np.array([
        [sx * lx, sy * ly, sz * lz]
        for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)
    ])
    edges = []
    for i, a in enumerate(corners):
        for j, b in enumerate(corners):
            if j > i and np.sum(a != b) == 1:
                edges.append((a, b))
    for edge_index, (a, b) in enumerate(edges):
        ax.plot(
            [a[0], b[0]], [a[1], b[1]], [a[2], b[2]],
            color=color, ls=linestyle, lw=linewidth, alpha=0.65,
            label=label if edge_index == 0 else None,
        )


def create_figure():
    arrays = np.load(RESULTS / "copyBF_long_3d_boundary_disturbance.npz")
    summary = json.loads(
        (RESULTS / "copyBF_long_3d_boundary_disturbance.json").read_text(
            encoding="ascii"
        )
    )
    hard = float(summary["hard_bound"])
    tightened = np.asarray(summary["terminal_tightened_bound"], dtype=float)
    reset_steps = [0, *arrays["scenario_reset_steps"].tolist(), len(arrays["reference"])]
    names = [
        "x reachable-limit pressure (0-40 s)",
        "y near-boundary pressure (40-80 s)",
        "z near-boundary pressure (80-120 s)",
    ]

    fig = plt.figure(figsize=(18, 6.4))
    for scenario, (start, stop, title) in enumerate(
        zip(reset_steps[:-1], reset_steps[1:], names), start=1
    ):
        ax = fig.add_subplot(1, 3, scenario, projection="3d")
        ref = arrays["reference"][start:stop]
        mpc = arrays["deterministic_mpc"][start:stop]
        smpc = arrays["smpc"][start:stop]

        _draw_bound_box(
            ax, [hard, hard, hard], "#1a1a1a", "-.", 0.8, "hard bound",
        )
        _draw_bound_box(
            ax, tightened, "#7b2cbf", "--", 0.75,
            "terminal tightened mean bound",
        )
        ax.plot(
            ref[:, 0], ref[:, 1], ref[:, 2],
            color="#1a1a1a", ls="--", lw=1.4, label="reference", zorder=2,
        )
        ax.plot(
            mpc[:, 0], mpc[:, 1], mpc[:, 2],
            color="#c44536", ls=":", lw=1.8,
            label="deterministic MPC", zorder=3,
        )
        ax.plot(
            smpc[:, 0], smpc[:, 1], smpc[:, 2],
            color="#087e8b", lw=1.2, label="SMPC", zorder=4,
        )
        ax.scatter(
            smpc[0, 0], smpc[0, 1], smpc[0, 2],
            color="#087e8b", marker="o", s=28, label="scenario start", zorder=5,
        )
        ax.scatter(
            smpc[-1, 0], smpc[-1, 1], smpc[-1, 2],
            color="#087e8b", marker="x", s=36, label="scenario end", zorder=5,
        )
        ax.set_xlim(-hard, hard)
        ax.set_ylim(-hard, hard)
        ax.set_zlim(-hard, hard)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("standardized x")
        ax.set_ylabel("standardized y")
        ax.set_zlabel("standardized z")
        ax.set_title(title, fontsize=11)
        ax.view_init(elev=24, azim=-58)
        ax.grid(True, alpha=0.22)
        if scenario == 1:
            ax.legend(loc="upper left", fontsize=7)

    fig.suptitle(
        "copyBF three-dimensional trajectories: three independent 40 s scenarios",
        fontsize=14,
    )
    fig.text(
        0.5, 0.015,
        "The state is reset at 40 s and 80 s; panels are not a continuous flight path.",
        ha="center", fontsize=10,
    )
    fig.subplots_adjust(left=0.02, right=0.99, bottom=0.08, top=0.88, wspace=0.02)
    output = RESULTS / "copyBF_long_3d_trajectories.png"
    fig.savefig(output, dpi=190)
    plt.close(fig)
    return output


if __name__ == "__main__":
    print(create_figure())
