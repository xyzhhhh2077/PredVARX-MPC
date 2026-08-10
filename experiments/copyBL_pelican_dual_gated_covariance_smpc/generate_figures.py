"""Generate BJ-style four evidence figures for the completed copyBL run."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
MEDIA = HERE / "results_media"
DATA = RESULTS / "copyBL_dual_gated_run.npz"
EVENTS = RESULTS / "copyBL_dual_gate_events.json"
SUMMARY = RESULTS / "copyBL_dual_gated_summary.json"
BJ_FIGURES = HERE.parent / "copyBJ_pelican_online_covariance_smpc" / "generate_smpc_figures.py"

_spec = importlib.util.spec_from_file_location("copybl_bj_figure_style", BJ_FIGURES)
BJ = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = BJ
_spec.loader.exec_module(BJ)

BLUE = BJ.BLUE
GREEN = BJ.GREEN
RED = BJ.RED
BLACK = BJ.BLACK
GRAY = BJ.GRAY
PURPLE = "#CC79A7"


def save_xyz(path, time_s, reference, task, hard_bound):
    figure, axes = plt.subplots(
        3, 1, figsize=(12.0, 8.0), sharex=True, constrained_layout=True,
    )
    for index, (axis, label, color) in enumerate(zip(axes, "xyz", BJ.AXIS_COLORS)):
        axis.plot(
            time_s, reference[:, index], color=BLACK, linestyle="--", linewidth=1.15,
            label=f"reference {label} = +/-3.46",
        )
        axis.plot(
            time_s, task[:, index], color=color, linewidth=1.0,
            label=f"dual-gated SMPC {label}",
        )
        axis.axhline(
            hard_bound, color=RED, linestyle=":", linewidth=1.0,
            label="hard bound +/-3.8" if index == 0 else None,
        )
        axis.axhline(-hard_bound, color=RED, linestyle=":", linewidth=1.0)
        axis.set_ylabel(f"{label} (standardized)")
        axis.set_ylim(-4.1, 4.1)
        axis.legend(loc="upper right", ncol=3, fontsize=8, frameon=False)
        BJ._style_axis(axis)
    axes[-1].set_xlabel("time (s)")
    figure.suptitle(
        "copyBL dual-gated SMPC: task trajectory and reference",
        fontsize=13, fontweight="bold",
    )
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def save_3d(path, reference, task, hard_bound):
    figure = plt.figure(figsize=(9.0, 7.5), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    axis.plot(*reference.T, color=BLACK, linestyle="--", linewidth=1.1, label="reference")
    axis.plot(*task.T, color=BLUE, linewidth=1.35, label="dual-gated SMPC")
    BJ.BI.BH._draw_bound_box(axis, hard_bound, RED, ":", "hard box +/-3.8")
    axis.scatter(*task[0], color=GREEN, s=35, label="start", depthshade=False)
    axis.set_xlabel("x (standardized)")
    axis.set_ylabel("y (standardized)")
    axis.set_zlabel("z (standardized)")
    axis.set_xlim(-4.0, 4.0)
    axis.set_ylim(-4.0, 4.0)
    axis.set_zlim(-4.0, 4.0)
    axis.set_box_aspect((1, 1, 1))
    axis.view_init(elev=23, azim=-52)
    axis.set_title(
        "copyBL dual-gated SMPC six-face trajectory in the hard output box",
        fontsize=13, fontweight="bold",
    )
    axis.legend(loc="upper left", fontsize=9, frameon=False)
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def save_noise(path, time_s, innovation, wind_velocity, innovation_std, sigma_wind):
    figure, axes = plt.subplots(
        2, 1, figsize=(12.0, 7.0), sharex=True, constrained_layout=True,
    )
    stride = max(1, len(time_s) // 4000)
    shown = slice(None, None, stride)
    for index in range(innovation.shape[1]):
        axes[0].plot(
            time_s[shown], innovation[shown, index], linewidth=0.55, alpha=0.82,
            color=BJ.AXIS_COLORS[index % 3], label=f"latent mode {index + 1}",
        )
    axes[0].text(
        0.01, 0.96,
        "injected plant innovation std (sample): "
        + ", ".join(f"z{index + 1}={value:.4f}" for index, value in enumerate(innovation_std)),
        transform=axes[0].transAxes, va="top", fontsize=8,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )
    axes[0].axhline(0.0, color=GRAY, linewidth=0.6)
    axes[0].set_ylabel("latent innovation")
    axes[0].set_title(
        "(a) Sampled composite latent innovation epsilon[k] injected into the plant",
        loc="left", fontweight="bold",
    )
    axes[0].legend(loc="upper right", ncol=3, fontsize=8, frameon=False)
    BJ._style_axis(axes[0])
    for index, label in enumerate("xyz"):
        axes[1].plot(
            time_s[shown], wind_velocity[shown, index], linewidth=0.6, alpha=0.82,
            color=BJ.AXIS_COLORS[index], label=f"wind {label}",
        )
    axes[1].axhline(
        sigma_wind, color=RED, linestyle=":", linewidth=1.0,
        label=f"long-run +/-sigma = {sigma_wind:.3f} m/s",
    )
    axes[1].axhline(-sigma_wind, color=RED, linestyle=":", linewidth=1.0)
    axes[1].axhline(0.0, color=GRAY, linewidth=0.6)
    axes[1].set_ylabel("wind velocity (m/s)")
    axes[1].set_xlabel("time (s)")
    axes[1].set_title("(b) Injected physical-wind realization w[k]", loc="left", fontweight="bold")
    axes[1].legend(loc="upper right", ncol=4, fontsize=8, frameon=False)
    BJ._style_axis(axes[1])
    figure.suptitle(
        "copyBL: two disturbance objects used by the simulation",
        fontsize=13, fontweight="bold",
    )
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def save_cost_with_gates(path, time_s, stage_cost, events, deadband):
    """Keep BJ's single cost panel and add a compact gate event rail."""
    figure, axis = plt.subplots(figsize=(12.0, 4.2), constrained_layout=True)
    axis.plot(time_s, stage_cost, color=BLUE, linewidth=0.85)
    axis.fill_between(time_s, 0.0, stage_cost, color=BLUE, alpha=0.12)
    axis.set_xlabel("time (s)")
    axis.set_ylabel("realized stage cost")
    axis.set_title(
        "copyBL realized stage cost and dual-gate decisions: "
        "Q||s[k]-r[k]||^2 + R||u[k]-u_mean||^2",
        fontweight="bold",
    )
    axis.text(
        0.99,
        0.95,
        f"mean = {np.mean(stage_cost):.4f}\nmax = {np.max(stage_cost):.4f}\n"
        f"value deadband = {100 * deadband:.1f}%",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": GRAY},
    )

    event_time = np.asarray([event["time_seconds"] for event in events], dtype=float)
    feasibility_reject = np.asarray(
        [not event["feasibility_accepted"] for event in events], dtype=bool,
    )
    value_events = [event for event in events if event["value_triggered"]]
    value_reject_time = [
        event["time_seconds"] for event in value_events if not event["value_accepted"]
    ]
    commit_time = [event["time_seconds"] for event in value_events if event["update_committed"]]
    ymin, ymax = axis.get_ylim()
    rail_y = ymin + 0.025 * (ymax - ymin)
    axis.scatter(
        event_time[feasibility_reject],
        np.full(np.sum(feasibility_reject), rail_y),
        marker="|",
        s=115,
        linewidths=1.6,
        color=RED,
        label="feasibility rejected",
        zorder=6,
    )
    axis.scatter(
        value_reject_time,
        np.full(len(value_reject_time), rail_y),
        marker="x",
        s=28,
        linewidths=1.0,
        color=PURPLE,
        label="value rejected",
        zorder=7,
    )
    axis.scatter(
        commit_time,
        np.full(len(commit_time), rail_y),
        marker="o",
        s=32,
        facecolor=GREEN,
        edgecolor=BLACK,
        linewidth=0.6,
        label="covariance committed",
        zorder=8,
    )
    axis.legend(loc="upper left", ncol=3, frameon=False, fontsize=8)
    BJ._style_axis(axis)
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def main():
    if not all(path.exists() for path in (DATA, EVENTS, SUMMARY)):
        raise FileNotFoundError("Run the complete copyBL experiment before generating figures")
    data = np.load(DATA)
    events = json.loads(EVENTS.read_text(encoding="utf-8"))
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    run = summary["copyBL"]
    if int(run["completed_steps"]) != int(summary["steps"]):
        raise RuntimeError("Source copyBL run did not complete")
    if int(run["fallback_count"]) != 0 or int(run["hard_violation_steps"]) != 0:
        raise RuntimeError("Source copyBL run used fallback or crossed the hard bound")

    task = np.asarray(data["trajectory_standardized"], dtype=float)
    reference = np.asarray(data["reference_standardized"], dtype=float)
    innovation = np.asarray(data["residual_innovation_latent"], dtype=float)
    wind = np.asarray(data["wind_velocity_mps"], dtype=float)
    cost = np.asarray(data["stage_cost"], dtype=float)
    time_s = np.arange(len(task), dtype=float) * float(summary["sample_time_seconds"])
    bound = float(summary["hard_bound_standardized"])
    if not all(len(array) == len(task) for array in (reference, innovation, wind, cost)):
        raise ValueError("Source arrays do not share a common step count")

    MEDIA.mkdir(parents=True, exist_ok=True)
    paths = {
        "xyz": MEDIA / "copyBL_smpc_xyz_reference.png",
        "trajectory_3d": MEDIA / "copyBL_smpc_3d_trajectory.png",
        "noise": MEDIA / "copyBL_smpc_noise_timeseries.png",
        "cost": MEDIA / "copyBL_smpc_stage_cost_dual_gates.png",
    }
    save_xyz(paths["xyz"], time_s, reference, task, bound)
    save_3d(paths["trajectory_3d"], reference, task, bound)
    save_noise(
        paths["noise"], time_s, innovation, wind,
        np.std(innovation, axis=0, ddof=1), BJ.BI.BH.SIGMA_WIND_MPS,
    )
    save_cost_with_gates(
        paths["cost"], time_s, cost, events, float(summary["value_deadband_fraction"]),
    )

    value_events = [event for event in events if event["value_triggered"]]
    manifest = {
        "style_reference": "copyBJ generate_smpc_figures.py",
        "figures": {key: value.name for key, value in paths.items()},
        "figure_count": 4,
        "feasibility_trigger_count": len(events),
        "feasibility_rejection_count": sum(not event["feasibility_accepted"] for event in events),
        "value_trigger_count": len(value_events),
        "value_accept_count": sum(event["value_accepted"] for event in value_events),
        "value_rejection_count": sum(not event["value_accepted"] for event in value_events),
    }
    (MEDIA / "copyBL_figure_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
