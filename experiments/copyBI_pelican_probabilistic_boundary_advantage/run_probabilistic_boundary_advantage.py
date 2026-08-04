"""copyBI: probabilistic hard-boundary advantage demonstration.

This experiment changes only the exogenous reference and random seeds. It reuses
copyBH's frozen identified Pelican model, deterministic hard-constrained MPC,
SMPC chance tightening, input box, covariance object, and physical-wind model.
"""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BH_SCRIPT = (
    HERE.parent / "copyBH_pelican_u_cap_ablation" / "run_u_cap_ablation.py"
)
os.environ.setdefault("COPYBH_INPUT_CAP", "7.0")
_spec = importlib.util.spec_from_file_location("copybi_frozen_copybh", BH_SCRIPT)
BH = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = BH
_spec.loader.exec_module(BH)

# Explicit aliases make the no-controller-change contract auditable.
CachedBoundaryController = BH.CachedBoundaryController
run_controller_pair = BH.run_controller_pair
innovation_sequence = BH.innovation_sequence
physical_wind_velocity = BH.physical_wind_velocity

DEFAULT_STEPS = 18_000
FACE_PRESSURES = np.array([3.46, 3.46, 3.46], dtype=float)
EXPERIMENT_HARD_BOUND = 3.80
WIND_SEED = 7016
INNOVATION_SEED = 16


def hard_edge_face_petal_reference(steps=DEFAULT_STEPS):
    """Frozen six-face petal with all axes equally close to the hard boundary."""
    return BH.face_petal_stress_reference(
        steps=int(steps),
        pressures=FACE_PRESSURES,
        slide=BH.FACE_PETAL_SLIDE,
    )


def realized_stage_cost(task, reference, control, control_mean,
                        q_weight=BH.Q_WEIGHT, ru=BH.RU):
    """Evaluate the realized one-step tracking and input penalty."""
    task = np.asarray(task, dtype=float)
    reference = np.asarray(reference, dtype=float)
    control = np.asarray(control, dtype=float)
    control_mean = np.asarray(control_mean, dtype=float).reshape(1, -1)
    if task.shape != reference.shape:
        raise ValueError("task and reference must have the same shape")
    if control.ndim != 2 or control.shape[0] != task.shape[0]:
        raise ValueError("control and task must have the same number of steps")
    tracking = float(q_weight) * np.sum((task - reference) ** 2, axis=1)
    effort = float(ru) * np.sum((control - control_mean) ** 2, axis=1)
    return tracking + effort


def validate_sampled_advantage(smpc, deterministic_mpc, expected_steps):
    """Fail closed unless the sampled event is a valid like-for-like contrast."""
    if smpc["disturbance_sha256"] != deterministic_mpc["disturbance_sha256"]:
        raise RuntimeError("Controllers did not use an identical disturbance")
    for name, result in (("SMPC", smpc), ("MPC", deterministic_mpc)):
        if result["completed_steps"] != int(expected_steps):
            raise RuntimeError(f"{name} did not complete the requested trajectory")
        if result["qp_failure_step"] is not None:
            raise RuntimeError(f"{name} QP failed during the sampled trajectory")
        if result["fallback_count"]:
            raise RuntimeError(f"{name} used a forbidden fallback")
    if smpc["violation_steps"] != 0:
        raise RuntimeError("SMPC crossed the hard bound in the sampled trajectory")
    if deterministic_mpc["violation_steps"] <= 0:
        raise RuntimeError("MPC did not cross the hard bound in the sampled trajectory")


def slow_z_hard_edge_reference(steps=12_000, edge=3.88):
    """Slow center→+z dwell→center→−z dwell→center hard-edge reference."""
    steps = int(steps)
    edge = float(edge)
    if steps < 100:
        raise ValueError("steps must be at least 100")
    if edge <= 0.0:
        raise ValueError("edge must be positive")

    knot_times = np.array([0.0, 0.20, 0.40, 0.50, 0.70, 0.90, 1.0])
    z_anchors = np.array([0.0, edge, edge, 0.0, -edge, -edge, 0.0])
    timeline = np.linspace(0.0, 1.0, steps, endpoint=True)
    segment = np.minimum(
        np.searchsorted(knot_times[1:], timeline, side="right"),
        len(knot_times) - 2,
    )
    left = knot_times[segment]
    right = knot_times[segment + 1]
    local = (timeline - left) / (right - left)
    blend = local ** 3 * (10.0 - 15.0 * local + 6.0 * local ** 2)
    z = (1.0 - blend) * z_anchors[segment] + blend * z_anchors[segment + 1]
    reference = np.zeros((steps, 3), dtype=float)
    reference[:, 2] = z
    return reference


def _controller_summary(result):
    return {
        "completed_steps": int(result["completed_steps"]),
        "qp_failure_step": result["qp_failure_step"],
        "fallback_count": int(result["fallback_count"]),
        "hard_violation_steps": int(result["violation_steps"]),
        "hard_violation_steps_per_axis": np.asarray(
            result["violation_steps_per_axis"], dtype=int,
        ).tolist(),
        "minimum_hard_margin": float(result["minimum_hard_margin"]),
        "peak_absolute_task": np.max(np.abs(result["S"]), axis=0).tolist(),
        "active_qp_steps": int(result["active_qp_steps"]),
        "maximum_qp_constraint_residual": float(
            result["maximum_qp_constraint_residual"]
        ),
    }


def _save_comparison_figure(path, reference, smpc, deterministic_mpc, hard_bound):
    import matplotlib.pyplot as plt

    blue = "#0072B2"
    red = "#D55E00"
    black = "#111827"
    gray = "#6b7280"
    s_smpc = np.asarray(smpc["S"], dtype=float)
    s_mpc = np.asarray(deterministic_mpc["S"], dtype=float)
    time_s = np.arange(len(reference), dtype=float) * BH.SAMPLE_TIME_SECONDS
    violation_by_axis = np.abs(s_mpc) > hard_bound
    violation_indices = np.flatnonzero(np.any(violation_by_axis, axis=1))
    if not violation_indices.size:
        raise RuntimeError("Cannot build crossing figure without an MPC violation")
    axis_counts = np.sum(violation_by_axis, axis=0)
    violation_axis = int(np.argmax(axis_counts))
    axis_name = "xyz"[violation_axis]
    axis_violation = violation_by_axis[:, violation_axis]
    axis_violation_indices = np.flatnonzero(axis_violation)
    center = int(axis_violation_indices[len(axis_violation_indices) // 2])
    zoom = slice(max(0, center - 120), min(len(reference), center + 121))
    zoom_indices = np.arange(len(reference))[zoom]
    zoom_violation = axis_violation[zoom]
    signed_bound = np.sign(s_mpc[center, violation_axis]) * hard_bound
    excess = np.max(np.abs(s_mpc[:, violation_axis])) - hard_bound

    margin_smpc = hard_bound - np.max(np.abs(s_smpc), axis=1)
    margin_mpc = hard_bound - np.max(np.abs(s_mpc), axis=1)

    figure = plt.figure(figsize=(13.2, 8.4), constrained_layout=True)
    grid = figure.add_gridspec(2, 3, height_ratios=[1.05, 0.95], width_ratios=[1.1, 1.5, 0.8])
    ax3d = figure.add_subplot(grid[0, 0], projection="3d")
    ax_zoom = figure.add_subplot(grid[0, 1:])
    ax_margin = figure.add_subplot(grid[1, :2])
    ax_bar = figure.add_subplot(grid[1, 2])

    ax3d.plot(*reference.T, "--", color=gray, linewidth=0.9, label="reference")
    ax3d.plot(*s_smpc.T, color=blue, linewidth=1.4, label="SMPC")
    ax3d.plot(*s_mpc.T, color=red, linewidth=1.0, alpha=0.88, label="hard MPC")
    ax3d.scatter(
        s_mpc[axis_violation_indices, 0], s_mpc[axis_violation_indices, 1],
        s_mpc[axis_violation_indices, 2], marker="x", s=55, linewidths=2.0,
        color="#b91c1c", depthshade=False, label="MPC outside hard",
    )
    BH._draw_bound_box(ax3d, hard_bound, black, ":", "hard box")
    ax3d.set_xlabel("x")
    ax3d.set_ylabel("y")
    ax3d.set_zlabel("z")
    ax3d.set_title("(a) Full six-face trajectory")
    ax3d.view_init(elev=22, azim=-52)
    ax3d.legend(loc="upper left", fontsize=7)

    ax_zoom.axhspan(
        signed_bound, signed_bound + np.sign(signed_bound) * 0.008,
        color="#fee2e2", alpha=0.9, label="outside hard bound",
    )
    ax_zoom.plot(
        time_s[zoom], s_smpc[zoom, violation_axis], color=blue,
        linewidth=2.0, label="SMPC",
    )
    ax_zoom.plot(
        time_s[zoom], s_mpc[zoom, violation_axis], color=red,
        linewidth=1.7, linestyle="--", label="hard MPC",
    )
    ax_zoom.axhline(
        signed_bound, color=black, linestyle=":", linewidth=1.6,
        label=f"hard bound = {signed_bound:.5f}",
    )
    ax_zoom.scatter(
        time_s[zoom_indices[zoom_violation]],
        s_mpc[zoom_indices[zoom_violation], violation_axis],
        marker="x", s=75, linewidths=2.2, color="#b91c1c", zorder=6,
        label=f"MPC violation ({axis_violation_indices.size} samples)",
    )
    local_values = np.concatenate([
        s_smpc[zoom, violation_axis], s_mpc[zoom, violation_axis],
        np.array([signed_bound]),
    ])
    pad = max(0.004, 0.18 * np.ptp(local_values))
    ax_zoom.set_ylim(np.min(local_values) - pad, np.max(local_values) + pad)
    ax_zoom.annotate(
        f"peak excess = {excess:.5f}",
        xy=(time_s[np.argmax(np.abs(s_mpc[:, violation_axis]))],
            s_mpc[np.argmax(np.abs(s_mpc[:, violation_axis])), violation_axis]),
        xytext=(0.03, 0.08), textcoords="axes fraction",
        arrowprops={"arrowstyle": "->", "color": "#b91c1c"},
        color="#991b1b", fontsize=10, fontweight="bold",
    )
    ax_zoom.set_title(f"(b) Magnified {axis_name}-axis boundary crossing")
    ax_zoom.set_xlabel("time (s)")
    ax_zoom.set_ylabel(f"{axis_name} (standardized)")
    ax_zoom.grid(True, alpha=0.25)
    ax_zoom.legend(fontsize=8, ncol=2, loc="best")

    ax_margin.plot(
        time_s[zoom], margin_smpc[zoom], color=blue, linewidth=2.0,
        label=f"SMPC (minimum {smpc['minimum_hard_margin']:.5f})",
    )
    ax_margin.plot(
        time_s[zoom], margin_mpc[zoom], color=red, linewidth=1.7,
        linestyle="--", label=f"hard MPC (minimum {deterministic_mpc['minimum_hard_margin']:.5f})",
    )
    ax_margin.axhline(0.0, color=black, linestyle=":", linewidth=1.5)
    ax_margin.fill_between(
        time_s[zoom], margin_mpc[zoom], 0.0, where=margin_mpc[zoom] < 0.0,
        color="#ef4444", alpha=0.55, interpolate=True, label="hard MPC outside",
    )
    ax_margin.scatter(
        time_s[axis_violation_indices], margin_mpc[axis_violation_indices],
        marker="x", s=55, linewidths=2.0, color="#b91c1c", zorder=5,
    )
    ax_margin.set_title("(c) Magnified instantaneous hard-bound margin")
    ax_margin.set_xlabel("time (s)")
    ax_margin.set_ylabel("hard bound - max |task|\n(negative means violation)")
    ax_margin.grid(True, alpha=0.25)
    ax_margin.legend(fontsize=8, ncol=2, loc="best")

    labels = ["SMPC", "hard MPC"]
    violations = [smpc["violation_steps"], deterministic_mpc["violation_steps"]]
    bars = ax_bar.bar(labels, violations, color=[blue, red], width=0.58)
    for bar, value in zip(bars, violations):
        ax_bar.text(
            bar.get_x() + bar.get_width() / 2.0, value + 0.12,
            str(int(value)), ha="center", va="bottom", fontweight="bold", fontsize=12,
        )
    ax_bar.set_ylim(0.0, max(7.0, max(violations) + 1.0))
    ax_bar.set_ylabel("hard-violation samples")
    ax_bar.set_title("(d) Same disturbance")
    ax_bar.grid(True, axis="y", alpha=0.25)
    ax_bar.text(
        0.5, 0.94,
        f"SMPC active QPs\n{smpc['active_qp_steps']} samples\n\n"
        f"MPC peak excess\n{excess:.5f}",
        transform=ax_bar.transAxes, ha="center", va="top", fontsize=9,
    )

    figure.suptitle(
        "Sampled advantage: hard MPC crosses; SMPC retains a positive margin\n"
        "Identical innovation and physical-wind realization",
        fontsize=13, fontweight="bold",
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def run_experiment(output_dir=None):
    """Reproduce the fixed sampled event and write data, metrics, and figure."""
    steps = DEFAULT_STEPS
    reference = hard_edge_face_petal_reference(steps)
    model, E, source_hard_bound, scales, noise = BH.load_identification_and_noise_objects()
    hard_bound = float(EXPERIMENT_HARD_BOUND)
    if hard_bound >= float(source_hard_bound):
        raise RuntimeError("copyBI hard bound must be stricter than the source bound")
    covariance = BH.build_process_covariances(
        model, E, scales["y_scale"][:3], BH.SIGMA_WIND_MPS,
    )
    residual_innovation = innovation_sequence(model, steps, INNOVATION_SEED)
    wind_velocity_mps = physical_wind_velocity(
        steps, BH.SIGMA_WIND_MPS, WIND_SEED,
    )
    task_map = E.T @ model["P"][:, :model["R"].shape[1]]
    wind_latent_increment, wind_map = BH.wind_velocity_to_latent_increment(
        wind_velocity_mps, task_map, scales["y_scale"][:3],
    )
    if not np.allclose(wind_map, covariance["G_w"], atol=1e-14):
        raise RuntimeError("Wind realization and covariance use different maps")

    smpc, deterministic_mpc = run_controller_pair(
        model, E, hard_bound, reference, residual_innovation,
        wind_latent_increment, covariance["Sigma_eps_aug"],
        initial_task=np.zeros(3),
    )
    validate_sampled_advantage(smpc, deterministic_mpc, steps)
    for name, result in (("SMPC", smpc), ("MPC", deterministic_mpc)):
        if result["maximum_qp_constraint_residual"] > 1e-6:
            raise RuntimeError(f"{name} QP residual exceeds 1e-6")

    results_dir = Path(output_dir) if output_dir else HERE / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    figure_path = results_dir / "copyBI_sampled_boundary_advantage.png"
    npz_path = results_dir / "copyBI_sampled_boundary_advantage.npz"
    json_path = results_dir / "copyBI_sampled_boundary_advantage.json"
    report_path = results_dir / "REPORT.md"

    np.savez_compressed(
        npz_path,
        reference_standardized=reference,
        smpc_trajectory_standardized=smpc["S"],
        deterministic_mpc_trajectory_standardized=deterministic_mpc["S"],
        smpc_control=smpc["U"],
        deterministic_mpc_control=deterministic_mpc["U"],
        smpc_active_history=smpc["active_history"],
        residual_innovation_latent=residual_innovation,
        wind_velocity_mps=wind_velocity_mps,
        wind_latent_increment=wind_latent_increment,
        hard_bound_standardized=np.array([hard_bound]),
        face_pressures_standardized=FACE_PRESSURES,
    )
    _save_comparison_figure(
        figure_path, reference, smpc, deterministic_mpc, hard_bound,
    )

    summary = {
        "experiment": "copyBI sampled probabilistic hard-boundary advantage",
        "scope": (
            "one deliberately selected sampled disturbance event; mechanism demonstration, "
            "not a Monte Carlo violation-probability estimate"
        ),
        "steps": steps,
        "sample_time_seconds": BH.SAMPLE_TIME_SECONDS,
        "hard_bound_standardized": float(hard_bound),
        "source_hard_bound_standardized": float(source_hard_bound),
        "face_pressures_standardized": FACE_PRESSURES.tolist(),
        "input_command_cap_standardized": BH.INPUT_COMMAND_BOUND_STANDARDIZED,
        "sigma_wind_mps": BH.SIGMA_WIND_MPS,
        "wind_seed": WIND_SEED,
        "innovation_seed": INNOVATION_SEED,
        "seed_selection": (
            "selected after an exploratory scan of 10 consecutive seed pairs. "
            "Six of the 10 deterministic-MPC probes ended at a non-optimal QP, "
            "so the scan is not retained or interpreted as a violation-frequency "
            "estimate. This pair was selected because both controllers completed "
            "all 18,000 steps and the sampled MPC crossing occurred."
        ),
        "shared_disturbance_sha256": smpc["disturbance_sha256"],
        "chance_covariance": (
            "identified Sigma_eps_aug only; physical wind is injected into both "
            "plants but is not added to the SMPC chance covariance"
        ),
        "smpc": _controller_summary(smpc),
        "deterministic_mpc": _controller_summary(deterministic_mpc),
        "claim_boundary": (
            "This selected run demonstrates one sampled outward disturbance crossing for "
            "hard-constrained nominal MPC while SMPC remains inside. It does not "
            "establish empirical calibration, recursive feasibility, or a robust "
            "hard guarantee under physical wind."
        ),
        "artifacts": {
            "data": npz_path.name,
            "figure": figure_path.name,
            "report": report_path.name,
        },
    }
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    report_path.write_text(
        "# copyBI — sampled hard-boundary advantage\n\n"
        "## What this run verifies\n\n"
        f"- Same innovation and physical-wind realization: `{smpc['disturbance_sha256']}`.\n"
        f"- Both controllers completed **{steps:,}** plant steps; fallback = 0.\n"
        f"- SMPC: **0** hard-violation steps, minimum margin "
        f"**{smpc['minimum_hard_margin']:.6f}**, active chance QPs "
        f"**{smpc['active_qp_steps']}**.\n"
        f"- Hard-constrained nominal MPC: **{deterministic_mpc['violation_steps']}** "
        f"hard-violation steps, minimum margin "
        f"**{deterministic_mpc['minimum_hard_margin']:.6f}**.\n\n"
        "## Claim boundary\n\n"
        "This seed pair was deliberately selected after an exploratory scan of 10 "
        "consecutive pairs. Six deterministic-MPC probes ended at a non-optimal QP, "
        "so that scan is not a valid violation-frequency estimate. The retained run "
        "is one full-length sampled crossing event, not a Monte Carlo estimate of "
        "the violation probability. The SMPC prediction covariance is the identified "
        "`Sigma_eps_aug`; physical wind is injected into both plants but is not added "
        "to that chance covariance. Therefore this run demonstrates the mechanism and "
        "a same-noise contrast, not probability calibration, recursive feasibility, "
        "or a robust hard guarantee under physical wind.\n",
        encoding="utf-8",
    )
    return summary


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    summary = run_experiment(output_dir=args.output_dir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    _main()
