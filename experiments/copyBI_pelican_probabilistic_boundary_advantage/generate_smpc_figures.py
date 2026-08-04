"""Run the fixed copyBI SMPC case and generate its four evidence figures."""

import importlib.util
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results_smpc"
SCRIPT = HERE / "run_probabilistic_boundary_advantage.py"
spec = importlib.util.spec_from_file_location("copybi_figure_source", SCRIPT)
BI = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = BI
spec.loader.exec_module(BI)

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
    figure, axes = plt.subplots(3, 1, figsize=(12.0, 8.0), sharex=True,
                                constrained_layout=True)
    labels = ["x", "y", "z"]
    for i, axis in enumerate(axes):
        axis.plot(time_s, reference[:, i], color=BLACK, linestyle="--",
                  linewidth=1.15, label=f"reference {labels[i]} = +/-3.46")
        axis.plot(time_s, task[:, i], color=AXIS_COLORS[i], linewidth=1.0,
                  label=f"SMPC {labels[i]}")
        axis.axhline(hard_bound, color=RED, linestyle=":", linewidth=1.0,
                     label="hard bound +/-3.8" if i == 0 else None)
        axis.axhline(-hard_bound, color=RED, linestyle=":", linewidth=1.0)
        axis.set_ylabel(f"{labels[i]} (standardized)")
        axis.set_ylim(-4.1, 4.1)
        axis.legend(loc="upper right", ncol=3, fontsize=8, frameon=False)
        _style_axis(axis)
    axes[-1].set_xlabel("time (s)")
    figure.suptitle("SMPC task trajectory and reference", fontsize=13,
                    fontweight="bold")
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def save_3d(path, reference, task, hard_bound):
    figure = plt.figure(figsize=(9.0, 7.5), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    axis.plot(*reference.T, color=BLACK, linestyle="--", linewidth=1.1,
              label="reference")
    axis.plot(*task.T, color=BLUE, linewidth=1.35, label="SMPC")
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
    axis.set_title("SMPC six-face trajectory in the hard output box",
                   fontsize=13, fontweight="bold")
    axis.legend(loc="upper left", fontsize=9, frameon=False)
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def save_noise(path, time_s, innovation, wind_velocity, sigma_eps, sigma_wind):
    figure, axes = plt.subplots(2, 1, figsize=(12.0, 7.0), sharex=True,
                                constrained_layout=True)
    stride = max(1, len(time_s) // 4000)
    shown = slice(None, None, stride)
    for i in range(innovation.shape[1]):
        axes[0].plot(time_s[shown], innovation[shown, i], linewidth=0.55,
                     alpha=0.82, color=AXIS_COLORS[i % 3], label=f"latent mode {i + 1}")
    eps_std = np.sqrt(np.maximum(np.diag(sigma_eps), 0.0))
    axes[0].text(0.01, 0.96, "estimated std: " + ", ".join(
        f"z{i + 1}={value:.4f}" for i, value in enumerate(eps_std)
    ), transform=axes[0].transAxes, va="top", fontsize=8,
                 bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})
    axes[0].axhline(0.0, color=GRAY, linewidth=0.6)
    axes[0].set_ylabel("latent innovation")
    axes[0].set_title("(a) Sampled composite latent innovation epsilon[k]",
                      loc="left", fontweight="bold")
    axes[0].legend(loc="upper right", ncol=3, fontsize=8, frameon=False)
    _style_axis(axes[0])

    for i, name in enumerate("xyz"):
        axes[1].plot(time_s[shown], wind_velocity[shown, i], linewidth=0.6,
                     alpha=0.82, color=AXIS_COLORS[i], label=f"wind {name}")
    axes[1].axhline(sigma_wind, color=RED, linestyle=":", linewidth=1.0,
                    label=f"long-run +/-sigma = {sigma_wind:.3f} m/s")
    axes[1].axhline(-sigma_wind, color=RED, linestyle=":", linewidth=1.0)
    axes[1].axhline(0.0, color=GRAY, linewidth=0.6)
    axes[1].set_ylabel("wind velocity (m/s)")
    axes[1].set_xlabel("time (s)")
    axes[1].set_title("(b) Injected physical-wind realization w[k]",
                      loc="left", fontweight="bold")
    axes[1].legend(loc="upper right", ncol=4, fontsize=8, frameon=False)
    _style_axis(axes[1])

    figure.suptitle("Two disturbance objects used by the simulation",
                    fontsize=13, fontweight="bold")
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def save_cost(path, time_s, stage_cost):
    figure, axis = plt.subplots(figsize=(12.0, 4.2), constrained_layout=True)
    axis.plot(time_s, stage_cost, color=BLUE, linewidth=0.85)
    axis.fill_between(time_s, 0.0, stage_cost, color=BLUE, alpha=0.12)
    axis.set_xlabel("time (s)")
    axis.set_ylabel("realized stage cost")
    axis.set_title(
        "Realized stage cost: Q||s[k]-r[k]||^2 + R||u[k]-u_mean||^2",
        fontweight="bold",
    )
    axis.text(0.99, 0.95,
              f"mean = {np.mean(stage_cost):.4f}\nmax = {np.max(stage_cost):.4f}",
              transform=axis.transAxes, ha="right", va="top", fontsize=9,
              bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": GRAY})
    _style_axis(axis)
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def main():
    steps = BI.DEFAULT_STEPS
    reference = BI.hard_edge_face_petal_reference(steps)
    model, E, source_bound, scales, noise = BI.BH.load_identification_and_noise_objects()
    hard_bound = float(BI.EXPERIMENT_HARD_BOUND)
    covariance = BI.BH.build_process_covariances(
        model, E, scales["y_scale"][:3], BI.BH.SIGMA_WIND_MPS,
    )
    innovation = BI.innovation_sequence(model, steps, BI.INNOVATION_SEED)
    wind_velocity = BI.physical_wind_velocity(
        steps, BI.BH.SIGMA_WIND_MPS, BI.WIND_SEED,
    )
    task_map = E.T @ model["P"][:, :model["R"].shape[1]]
    wind_latent, wind_map = BI.BH.wind_velocity_to_latent_increment(
        wind_velocity, task_map, scales["y_scale"][:3],
    )
    if not np.allclose(wind_map, covariance["G_w"], atol=1e-14):
        raise RuntimeError("Wind realization and covariance use different maps")

    smpc = BI.BH.simulate_controller(
        model, E, hard_bound, reference, innovation, wind_latent,
        covariance["Sigma_eps_aug"], BI.BH.ALPHA_JOINT,
        initial_task=np.zeros(3),
    )
    if smpc["completed_steps"] != steps or smpc["qp_failure_step"] is not None:
        raise RuntimeError(f"SMPC failed at step {smpc['qp_failure_step']}")
    if smpc["fallback_count"] != 0 or smpc["violation_steps"] != 0:
        raise RuntimeError("SMPC run used fallback or crossed the hard bound")

    cost = BI.realized_stage_cost(
        smpc["S"], reference, smpc["U"], model["u_mean"].ravel(),
    )
    time_s = np.arange(steps, dtype=float) * BI.BH.SAMPLE_TIME_SECONDS
    RESULTS.mkdir(parents=True, exist_ok=True)
    paths = {
        "xyz": RESULTS / "copyBI_smpc_xyz_reference.png",
        "trajectory_3d": RESULTS / "copyBI_smpc_3d_trajectory.png",
        "noise": RESULTS / "copyBI_smpc_noise_timeseries.png",
        "cost": RESULTS / "copyBI_smpc_stage_cost.png",
        "data": RESULTS / "copyBI_smpc_figure_data.npz",
        "summary": RESULTS / "copyBI_smpc_figure_summary.json",
    }
    save_xyz(paths["xyz"], time_s, reference, smpc["S"], hard_bound)
    save_3d(paths["trajectory_3d"], reference, smpc["S"], hard_bound)
    save_noise(paths["noise"], time_s, innovation, wind_velocity,
               noise["Sigma_eps"], BI.BH.SIGMA_WIND_MPS)
    save_cost(paths["cost"], time_s, cost)
    np.savez_compressed(
        paths["data"], time_seconds=time_s, reference_standardized=reference,
        smpc_task_standardized=smpc["S"], smpc_control=smpc["U"],
        composite_latent_innovation=innovation,
        physical_wind_velocity_mps=wind_velocity,
        realized_stage_cost=cost, hard_bound_standardized=np.array([hard_bound]),
        input_bound_standardized=np.array([
            BI.BH.INPUT_COMMAND_BOUND_STANDARDIZED
        ]),
    )
    summary = {
        "steps": steps,
        "hard_bound_standardized": hard_bound,
        "source_hard_bound_standardized": float(source_bound),
        "reference_peak_standardized": BI.FACE_PRESSURES.tolist(),
        "hard_violation_steps": smpc["violation_steps"],
        "minimum_hard_margin": smpc["minimum_hard_margin"],
        "fallback_count": smpc["fallback_count"],
        "active_qp_steps": smpc["active_qp_steps"],
        "stage_cost_definition": "Q||s-r||^2 + R||u-u_mean||^2",
        "stage_cost_mean": float(np.mean(cost)),
        "stage_cost_maximum": float(np.max(cost)),
        "noise_scope": (
            "composite identified latent innovation and injected physical wind; "
            "not separately identified process-noise Qw and measurement-noise Rv"
        ),
        "artifacts": {key: value.name for key, value in paths.items()},
    }
    paths["summary"].write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
