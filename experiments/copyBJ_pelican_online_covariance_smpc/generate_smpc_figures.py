"""Generate BI-style four evidence figures for the verified copyBJ online SMPC.

Uses the already-saved beta=0.8 online trajectory in
results/copyBJ_online_covariance_comparison.npz. Does not re-run the controller.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results_smpc"
SOURCE_NPZ = HERE / "results" / "copyBJ_online_covariance_comparison.npz"
SOURCE_JSON = HERE / "results" / "copyBJ_online_covariance_summary.json"
BI_SCRIPT = (
    HERE.parent
    / "copyBI_pelican_probabilistic_boundary_advantage"
    / "run_probabilistic_boundary_advantage.py"
)

_spec = importlib.util.spec_from_file_location("copybj_bi_figures", BI_SCRIPT)
BI = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = BI
_spec.loader.exec_module(BI)

BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
RED = "#D55E00"
BLACK = "#111827"
GRAY = "#6B7280"
AXIS_COLORS = [BLUE, ORANGE, GREEN]


def _style_axis(axis):
    axis.grid(True, alpha=0.22, linewidth=0.6)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def save_xyz(path, time_s, reference, task, hard_bound):
    figure, axes = plt.subplots(
        3, 1, figsize=(12.0, 8.0), sharex=True, constrained_layout=True,
    )
    labels = ["x", "y", "z"]
    for i, axis in enumerate(axes):
        axis.plot(
            time_s,
            reference[:, i],
            color=BLACK,
            linestyle="--",
            linewidth=1.15,
            label=f"reference {labels[i]} = +/-3.46",
        )
        axis.plot(
            time_s,
            task[:, i],
            color=AXIS_COLORS[i],
            linewidth=1.0,
            label=f"online SMPC {labels[i]}",
        )
        axis.axhline(
            hard_bound,
            color=RED,
            linestyle=":",
            linewidth=1.0,
            label="hard bound +/-3.8" if i == 0 else None,
        )
        axis.axhline(-hard_bound, color=RED, linestyle=":", linewidth=1.0)
        axis.set_ylabel(f"{labels[i]} (standardized)")
        axis.set_ylim(-4.1, 4.1)
        axis.legend(loc="upper right", ncol=3, fontsize=8, frameon=False)
        _style_axis(axis)
    axes[-1].set_xlabel("time (s)")
    figure.suptitle(
        "copyBJ online SMPC (beta=0.8): task trajectory and reference",
        fontsize=13,
        fontweight="bold",
    )
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def save_3d(path, reference, task, hard_bound):
    figure = plt.figure(figsize=(9.0, 7.5), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    axis.plot(*reference.T, color=BLACK, linestyle="--", linewidth=1.1, label="reference")
    axis.plot(*task.T, color=BLUE, linewidth=1.35, label="online SMPC")
    BI.BH._draw_bound_box(axis, hard_bound, RED, ":", "hard box +/-3.8")
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
        "copyBJ online SMPC six-face trajectory in the hard output box",
        fontsize=13,
        fontweight="bold",
    )
    axis.legend(loc="upper left", fontsize=9, frameon=False)
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def save_noise(path, time_s, innovation, wind_velocity, sigma_eps_diag_std, sigma_wind):
    figure, axes = plt.subplots(
        2, 1, figsize=(12.0, 7.0), sharex=True, constrained_layout=True,
    )
    stride = max(1, len(time_s) // 4000)
    shown = slice(None, None, stride)
    for i in range(innovation.shape[1]):
        axes[0].plot(
            time_s[shown],
            innovation[shown, i],
            linewidth=0.55,
            alpha=0.82,
            color=AXIS_COLORS[i % 3],
            label=f"latent mode {i + 1}",
        )
    axes[0].text(
        0.01,
        0.96,
        "injected plant innovation std (sample): "
        + ", ".join(f"z{i + 1}={value:.4f}" for i, value in enumerate(sigma_eps_diag_std)),
        transform=axes[0].transAxes,
        va="top",
        fontsize=8,
        bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
    )
    axes[0].axhline(0.0, color=GRAY, linewidth=0.6)
    axes[0].set_ylabel("latent innovation")
    axes[0].set_title(
        "(a) Sampled composite latent innovation epsilon[k] injected into the plant",
        loc="left",
        fontweight="bold",
    )
    axes[0].legend(loc="upper right", ncol=3, fontsize=8, frameon=False)
    _style_axis(axes[0])

    for i, name in enumerate("xyz"):
        axes[1].plot(
            time_s[shown],
            wind_velocity[shown, i],
            linewidth=0.6,
            alpha=0.82,
            color=AXIS_COLORS[i],
            label=f"wind {name}",
        )
    axes[1].axhline(
        sigma_wind,
        color=RED,
        linestyle=":",
        linewidth=1.0,
        label=f"long-run +/-sigma = {sigma_wind:.3f} m/s",
    )
    axes[1].axhline(-sigma_wind, color=RED, linestyle=":", linewidth=1.0)
    axes[1].axhline(0.0, color=GRAY, linewidth=0.6)
    axes[1].set_ylabel("wind velocity (m/s)")
    axes[1].set_xlabel("time (s)")
    axes[1].set_title(
        "(b) Injected physical-wind realization w[k]",
        loc="left",
        fontweight="bold",
    )
    axes[1].legend(loc="upper right", ncol=4, fontsize=8, frameon=False)
    _style_axis(axes[1])

    figure.suptitle(
        "copyBJ: two disturbance objects used by the simulation",
        fontsize=13,
        fontweight="bold",
    )
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def save_cost(path, time_s, stage_cost):
    figure, axis = plt.subplots(figsize=(12.0, 4.2), constrained_layout=True)
    axis.plot(time_s, stage_cost, color=BLUE, linewidth=0.85)
    axis.fill_between(time_s, 0.0, stage_cost, color=BLUE, alpha=0.12)
    axis.set_xlabel("time (s)")
    axis.set_ylabel("realized stage cost")
    axis.set_title(
        "copyBJ realized stage cost: Q||s[k]-r[k]||^2 + R||u[k]-u_mean||^2",
        fontweight="bold",
    )
    axis.text(
        0.99,
        0.95,
        f"mean = {np.mean(stage_cost):.4f}\nmax = {np.max(stage_cost):.4f}",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": GRAY},
    )
    _style_axis(axis)
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def main():
    if not SOURCE_NPZ.exists():
        raise FileNotFoundError(
            f"Missing source NPZ from the verified online run: {SOURCE_NPZ}"
        )
    summary_src = json.loads(SOURCE_JSON.read_text(encoding="utf-8"))
    online = summary_src["online_copyBJ_smpc"]
    if int(online["completed_steps"]) != int(summary_src["steps"]):
        raise RuntimeError("Source online run did not complete all steps")
    if int(online["fallback_count"]) != 0 or int(online["hard_violation_steps"]) != 0:
        raise RuntimeError("Source online run used fallback or crossed the hard bound")

    data = np.load(SOURCE_NPZ)
    reference = np.asarray(data["reference_standardized"], dtype=float)
    task = np.asarray(data["online_trajectory_standardized"], dtype=float)
    control = np.asarray(data["online_control"], dtype=float)
    innovation = np.asarray(data["residual_innovation_latent"], dtype=float)
    wind_velocity = np.asarray(data["wind_velocity_mps"], dtype=float)
    steps = int(task.shape[0])
    if not all(arr.shape[0] == steps for arr in (reference, control, innovation, wind_velocity)):
        raise ValueError("Source arrays do not share a common step count")

    hard_bound = float(BI.EXPERIMENT_HARD_BOUND)
    model, _E, source_bound, _scales, _noise = BI.BH.load_identification_and_noise_objects()
    u_mean = model["u_mean"].ravel()
    cost = BI.realized_stage_cost(task, reference, control, u_mean)
    time_s = np.arange(steps, dtype=float) * BI.BH.SAMPLE_TIME_SECONDS
    innov_std = np.std(innovation, axis=0, ddof=1)

    RESULTS.mkdir(parents=True, exist_ok=True)
    paths = {
        "xyz": RESULTS / "copyBJ_smpc_xyz_reference.png",
        "trajectory_3d": RESULTS / "copyBJ_smpc_3d_trajectory.png",
        "noise": RESULTS / "copyBJ_smpc_noise_timeseries.png",
        "cost": RESULTS / "copyBJ_smpc_stage_cost.png",
        "data": RESULTS / "copyBJ_smpc_figure_data.npz",
        "summary": RESULTS / "copyBJ_smpc_figure_summary.json",
    }
    save_xyz(paths["xyz"], time_s, reference, task, hard_bound)
    save_3d(paths["trajectory_3d"], reference, task, hard_bound)
    save_noise(
        paths["noise"],
        time_s,
        innovation,
        wind_velocity,
        innov_std,
        BI.BH.SIGMA_WIND_MPS,
    )
    save_cost(paths["cost"], time_s, cost)

    abs_task = np.abs(task)
    np.savez_compressed(
        paths["data"],
        time_seconds=time_s,
        reference_standardized=reference,
        smpc_task_standardized=task,
        smpc_control=control,
        composite_latent_innovation=innovation,
        physical_wind_velocity_mps=wind_velocity,
        realized_stage_cost=cost,
        hard_bound_standardized=np.array([hard_bound]),
        input_bound_standardized=np.array([BI.BH.INPUT_COMMAND_BOUND_STANDARDIZED]),
        online_estimated_sigma_eps_std=np.asarray(
            data["online_estimated_sigma_eps_std"], dtype=float,
        ),
        online_stage_tightening=np.asarray(data["online_stage_tightening"], dtype=float),
        online_weight=np.array([float(summary_src["shrinkage_candidate"]["online_weight"])]),
    )
    summary = {
        "experiment": "copyBJ online covariance-adaptive SMPC (beta=0.8 candidate)",
        "source_npz": SOURCE_NPZ.name,
        "steps": steps,
        "hard_bound_standardized": hard_bound,
        "source_hard_bound_standardized": float(source_bound),
        "reference_peak_standardized": BI.FACE_PRESSURES.tolist(),
        "hard_violation_steps": int(online["hard_violation_steps"]),
        "minimum_hard_margin": float(online["minimum_hard_margin_standardized"]),
        "fallback_count": int(online["fallback_count"]),
        "active_qp_steps": int(online["active_qp_steps"]),
        "online_weight": float(summary_src["shrinkage_candidate"]["online_weight"]),
        "offline_prior_weight": float(
            summary_src["shrinkage_candidate"]["offline_prior_weight"]
        ),
        "task_peak_magnitudes": abs_task.max(axis=0).tolist(),
        "stage_cost_definition": "Q||s-r||^2 + R||u-u_mean||^2",
        "stage_cost_mean": float(np.mean(cost)),
        "stage_cost_maximum": float(np.max(cost)),
        "noise_scope": (
            "composite identified latent innovation and injected physical wind; "
            "not separately identified process-noise Qw and measurement-noise Rv; "
            "controller uses online composite innovation covariance with offline prior shrinkage"
        ),
        "claim_boundary": summary_src["claim_boundary"],
        "artifacts": {key: value.name for key, value in paths.items()},
    }
    paths["summary"].write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
