"""copyBJ: copyBI boundary event with copyAU-style online covariance learning.

The Pelican geometry and dynamics (E, A, B, P, R) remain frozen. At each plant
sample, the newest closed-loop innovation enters a 40-sample window. The window
covariance updates the SMPC chance tightening at the next control decision.
"""

import argparse
from collections import deque
import importlib.util
import json
import os
from pathlib import Path
import sys

import cvxpy as cp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm


HERE = Path(__file__).resolve().parent
BI_SCRIPT = (
    HERE.parent / "copyBI_pelican_probabilistic_boundary_advantage"
    / "run_probabilistic_boundary_advantage.py"
)
os.environ.setdefault("COPYBH_INPUT_CAP", "7.0")
_spec = importlib.util.spec_from_file_location("copybj_copybi", BI_SCRIPT)
BI = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = BI
_spec.loader.exec_module(BI)
BH = BI.BH

WINDOW_SAMPLES = 40
UPDATE_MIN_SAMPLES = 5
COVARIANCE_FLOOR = 1e-8
ONLINE_UPDATED_OBJECTS = ("Sigma_eps", "Sigma_obs_proxy")
FIXED_OBJECTS = ("E", "A", "B", "P", "R")
DEFAULT_STEPS = BI.DEFAULT_STEPS
SHRUNK_ONLINE_WEIGHT = 0.8


class OnlineInnovationCovariance:
    """copyAU-compatible finite-window covariance with an offline prior."""

    def __init__(self, dimension, window, initial_covariance,
                 floor=COVARIANCE_FLOOR, min_samples=UPDATE_MIN_SAMPLES):
        self.dimension = int(dimension)
        self.window = int(window)
        self.floor = float(floor)
        self.min_samples = int(min_samples)
        self.initial_covariance = np.asarray(initial_covariance, dtype=float).copy()
        if self.initial_covariance.shape != (self.dimension, self.dimension):
            raise ValueError("initial_covariance has the wrong shape")
        self._samples = deque(maxlen=self.window)
        self.covariance = self.initial_covariance.copy()

    @property
    def count(self):
        return len(self._samples)

    def update(self, sample):
        sample = np.asarray(sample, dtype=float).reshape(self.dimension)
        self._samples.append(sample.copy())
        if self.count >= self.min_samples:
            values = np.asarray(self._samples, dtype=float)
            centered = values - values.mean(axis=0, keepdims=True)
            covariance = centered.T @ centered / max(self.count - 1, 1)
            self.covariance = BH._symmetrize_psd(
                covariance + self.floor * np.eye(self.dimension),
                floor=self.floor,
            )
        return self.covariance.copy()


def shrink_online_covariance(prior, online, online_weight):
    """Shrink a short-window estimate toward the offline training covariance."""
    weight = float(online_weight)
    if not 0.0 <= weight <= 1.0:
        raise ValueError("online_weight must lie in [0, 1]")
    prior = np.asarray(prior, dtype=float)
    online = np.asarray(online, dtype=float)
    if prior.shape != online.shape or prior.ndim != 2 or prior.shape[0] != prior.shape[1]:
        raise ValueError("prior and online covariance must be equal square matrices")
    return BH._symmetrize_psd((1.0 - weight) * prior + weight * online)


class AdaptiveBoundaryController:
    """Parameterized copyBI QP whose chance bounds follow online covariance."""

    def __init__(self, model, E, task_bound, process_covariance,
                 horizon=BH.HORIZON_STEPS, q_weight=BH.Q_WEIGHT, ru=BH.RU,
                 alpha_joint=BH.ALPHA_JOINT,
                 control_interval=BH.CONTROL_INTERVAL_STEPS):
        self.model = model
        self.E = np.asarray(E, dtype=float)
        self.N = int(horizon)
        self.d = int(control_interval)
        self.nu = model["B"].shape[1]
        self.n = model["A"].shape[0]
        self.task_axes = self.E.shape[1]
        self.task_bound = float(task_bound)
        self.alpha_joint = float(alpha_joint)
        self.u_min = np.full(self.nu, -BH.INPUT_COMMAND_BOUND_STANDARDIZED)
        self.u_max = np.full(self.nu, BH.INPUT_COMMAND_BOUND_STANDARDIZED)
        self.u_mean = model["u_mean"].ravel()
        self.y_mean = model["y_mean"].ravel()
        self.Ad, self.Bd = BH.mstep_lift(model["A"], model["B"], self.d)
        per_face = self.alpha_joint / (2.0 * self.task_axes * self.N)
        self.quantile = float(norm.ppf(1.0 - per_face))

        self.U = cp.Variable((self.N, self.nu), name="adaptive_control_plan")
        self.state_parameter = cp.Parameter(self.n, name="adaptive_current_state")
        self.reference_parameter = cp.Parameter(
            (self.N, self.task_axes), name="adaptive_reference",
        )
        self.stage_bound_parameter = cp.Parameter(
            (self.N, self.task_axes), nonneg=True, name="adaptive_stage_bound",
        )
        constraints = [self.U >= self.u_min, self.U <= self.u_max]
        objective = float(ru) * cp.sum_squares(
            self.U - np.tile(self.u_mean, (self.N, 1))
        )
        state = self.state_parameter
        self.task_expressions = []
        self.task_constraints = []
        for stage in range(self.N):
            state = self.Ad @ state + self.Bd @ (self.U[stage] - self.u_mean)
            task = self.E.T @ (self.y_mean + model["P"] @ state)
            self.task_expressions.append(task)
            objective += float(q_weight) * cp.sum_squares(
                task - self.reference_parameter[stage]
            )
            bound = self.stage_bound_parameter[stage]
            upper = task <= bound
            lower = task >= -bound
            constraints.extend((upper, lower))
            self.task_constraints.append((upper, lower))
        self.problem = cp.Problem(cp.Minimize(objective), constraints)
        if not self.problem.is_dpp():
            raise RuntimeError("Adaptive control problem must satisfy DPP")
        self.update_process_covariance(process_covariance)

    def update_process_covariance(self, process_covariance):
        process_covariance = BH.lift_noise_covariance(
            self.model["A"], BH._symmetrize_psd(process_covariance), self.d,
        )
        covariance = np.zeros((self.n, self.n))
        tightening = []
        for _ in range(self.N):
            covariance = BH._symmetrize_psd(
                self.Ad @ covariance @ self.Ad.T + process_covariance
            )
            task_covariance = (
                self.E.T @ self.model["P"] @ covariance
                @ self.model["P"].T @ self.E
            )
            variances = np.maximum(np.diag(task_covariance), 1e-15)
            tightening.append(self.quantile * np.sqrt(variances))
        self.stage_tightening = np.asarray(tightening)
        self.stage_mean_bound = self.task_bound - self.stage_tightening
        if np.any(self.stage_mean_bound <= 0.0):
            raise ValueError("Online chance tightening leaves an empty task interval")
        self.stage_bound_parameter.value = self.stage_mean_bound

    def step(self, state, reference_horizon):
        self.state_parameter.value = np.asarray(state, dtype=float).reshape(self.n)
        reference_horizon = np.asarray(reference_horizon, dtype=float)
        if reference_horizon.shape != (self.N, self.task_axes):
            raise ValueError("reference_horizon has the wrong shape")
        self.reference_parameter.value = reference_horizon
        retry_count = 0
        try:
            self.problem.solve(
                solver=cp.OSQP, warm_start=True, verbose=False,
                eps_abs=1e-8, eps_rel=1e-8, max_iter=20000, polishing=True,
            )
            if self.problem.status == cp.USER_LIMIT:
                retry_count = 1
                self.problem.solve(
                    solver=cp.OSQP, warm_start=False, verbose=False,
                    eps_abs=1e-8, eps_rel=1e-8, max_iter=100000,
                    polishing=True,
                )
        except cp.error.SolverError as exc:
            return None, {"status": f"solver_error: {exc}", "retry_count": retry_count}
        if self.problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
            return None, {"status": self.problem.status, "retry_count": retry_count}
        planned = np.asarray(self.U.value, dtype=float)
        predicted = np.vstack([
            np.asarray(expr.value, dtype=float).reshape(1, -1)
            for expr in self.task_expressions
        ])
        slacks = self.stage_mean_bound - np.abs(predicted)
        duals = []
        for upper, lower in self.task_constraints:
            duals.extend(np.asarray(upper.dual_value, dtype=float).ravel())
            duals.extend(np.asarray(lower.dual_value, dtype=float).ravel())
        residual = max(
            float(np.max(self.u_min - planned)),
            float(np.max(planned - self.u_max)),
            float(np.max(np.abs(predicted) - self.stage_mean_bound)),
        )
        return np.clip(planned[0], self.u_min, self.u_max), {
            "status": self.problem.status,
            "minimum_predicted_slack": float(np.min(slacks)),
            "maximum_constraint_dual": float(np.max(duals)),
            "maximum_constraint_residual": residual,
            "solve_time_seconds": float(self.problem.solver_stats.solve_time or 0.0),
            "retry_count": retry_count,
        }


def simulate_online_controller(model, E, hard_bound, reference,
                               residual_innovation, wind_latent_increment,
                               initial_covariance, initial_obs_covariance,
                               initial_task=None, online_weight=1.0):
    """Run copyBI with copyAU's 40-sample online residual updates."""
    reference = np.asarray(reference, dtype=float)
    residual_innovation = np.asarray(residual_innovation, dtype=float)
    wind_latent_increment = np.asarray(wind_latent_increment, dtype=float)
    steps = len(reference)
    ell = model["A"].shape[0] // 2
    p = model["P"].shape[0]
    prior_covariance = np.asarray(initial_covariance[:ell, :ell], dtype=float)
    controller = AdaptiveBoundaryController(
        model, E, hard_bound, initial_covariance,
    )
    eps_estimator = OnlineInnovationCovariance(
        ell, WINDOW_SAMPLES, initial_covariance[:ell, :ell]
    )
    obs_estimator = OnlineInnovationCovariance(
        p, WINDOW_SAMPLES, initial_obs_covariance
    )

    if initial_task is None:
        initial_task = reference[0]
    state, control = BH.equilibrium_initial_condition(model, E, initial_task)
    previous_state = None
    previous_control = None
    task_history = []
    control_history = []
    eps_trace = []
    obs_trace = []
    tightening_history = []
    slacks = []
    duals = []
    residuals = []
    solve_times = []
    statuses = {}
    solver_retry_count = 0
    qp_failure_step = None
    IminusPR = np.eye(p) - model["P"][:, :ell] @ model["R"].T

    for step in range(steps):
        output = model["y_mean"] + model["P"] @ state
        task = (E.T @ output).ravel()
        if previous_state is not None:
            predicted = (
                model["A"] @ previous_state
                + model["B"] @ (previous_control.reshape(-1, 1) - model["u_mean"])
            )
            innovation = state[:ell, 0] - predicted[:ell, 0]
            eps_estimator.update(innovation)
        projection_residual = (
            IminusPR @ (output - model["y_mean"])
        ).ravel()
        obs_estimator.update(projection_residual)
        eps_trace.append(float(np.trace(eps_estimator.covariance) / ell) ** 0.5)
        obs_trace.append(float(np.trace(obs_estimator.covariance) / p) ** 0.5)

        if step % BH.CONTROL_INTERVAL_STEPS == 0:
            covariance = np.zeros_like(initial_covariance)
            covariance[:ell, :ell] = shrink_online_covariance(
                prior_covariance, eps_estimator.covariance, online_weight,
            )
            try:
                controller.update_process_covariance(covariance)
            except ValueError as exc:
                qp_failure_step = int(step)
                status = f"empty_chance_interval: {exc}"
                statuses[status] = statuses.get(status, 0) + 1
                break
            refs = BH.future_reference_horizon(
                reference, step, controller.N,
                control_interval=BH.CONTROL_INTERVAL_STEPS,
            )
            candidate, info = controller.step(state, refs)
            solver_retry_count += int(info.get("retry_count", 0))
            status = info["status"]
            statuses[status] = statuses.get(status, 0) + 1
            if candidate is None:
                qp_failure_step = int(step)
                break
            control = candidate
            tightening_history.append(controller.stage_tightening.copy())
            slacks.append(info["minimum_predicted_slack"])
            duals.append(info["maximum_constraint_dual"])
            residuals.append(info["maximum_constraint_residual"])
            solve_times.append(info["solve_time_seconds"])

        task_history.append(task)
        control_history.append(control.copy())
        previous_state = state.copy()
        previous_control = control.copy()
        state = (
            model["A"] @ state
            + model["B"] @ (control.reshape(-1, 1) - model["u_mean"])
        )
        state[:ell, 0] += residual_innovation[step]
        state[:ell, 0] += wind_latent_increment[step]

    S = np.asarray(task_history, dtype=float).reshape(-1, 3)
    U = np.asarray(control_history, dtype=float).reshape(-1, model["B"].shape[1])
    completed = len(S)
    hard_margin = float(hard_bound) - np.abs(S)
    axis_violations = np.abs(S) > float(hard_bound)
    step_violations = np.any(axis_violations, axis=1)
    slack_array = np.asarray(slacks)
    dual_array = np.asarray(duals)
    active = (slack_array <= 5e-4) | (dual_array > 1e-6)
    saturation = np.any(
        np.abs(U) >= BH.INPUT_COMMAND_BOUND_STANDARDIZED - 1e-4, axis=1
    )
    return {
        "S": S,
        "U": U,
        "completed_steps": completed,
        "qp_count": len(slacks) + int(qp_failure_step is not None),
        "qp_failure_count": int(qp_failure_step is not None),
        "qp_failure_step": qp_failure_step,
        "qp_failure_status": "" if qp_failure_step is None else status,
        "fallback_count": 0,
        "qp_status_counts": statuses,
        "violation_steps": int(np.sum(step_violations)),
        "minimum_hard_margin": float(np.min(hard_margin)) if completed else float("nan"),
        "minimum_hard_margin_per_axis": np.min(hard_margin, axis=0),
        "rmse": np.sqrt(np.mean((S - reference[:completed]) ** 2, axis=0)),
        "input_saturation_steps": int(np.sum(saturation)),
        "active_qp_steps": int(np.sum(active)),
        "positive_dual_qp_steps": int(np.sum(dual_array > 1e-6)),
        "maximum_qp_constraint_residual": (
            float(np.max(residuals)) if residuals else float("nan")
        ),
        "mean_qp_solve_time_seconds": (
            float(np.mean(solve_times)) if solve_times else float("nan")
        ),
        "disturbance_sha256": BH._disturbance_sha256(
            residual_innovation, wind_latent_increment
        ),
        "estimated_sigma_eps_std": np.asarray(eps_trace),
        "estimated_sigma_obs_std": np.asarray(obs_trace),
        "stage_tightening_history": np.asarray(tightening_history),
        "final_Sigma_eps": eps_estimator.covariance,
        "final_Sigma_obs_proxy": obs_estimator.covariance,
        "online_weight": float(online_weight),
        "solver_retry_count": int(solver_retry_count),
    }


def _summary(result):
    return {
        "completed_steps": int(result["completed_steps"]),
        "qp_failure_count": int(result["qp_failure_count"]),
        "fallback_count": int(result["fallback_count"]),
        "hard_violation_steps": int(result["violation_steps"]),
        "minimum_hard_margin_standardized": float(result["minimum_hard_margin"]),
        "rmse_standardized": np.asarray(result["rmse"]).tolist(),
        "active_qp_steps": int(result["active_qp_steps"]),
        "input_saturation_steps": int(result["input_saturation_steps"]),
        "mean_qp_solve_time_seconds": float(result["mean_qp_solve_time_seconds"]),
        "solver_retry_count": int(result.get("solver_retry_count", 0)),
    }


def _save_figure(path, reference, frozen, online, hard_bound):
    time = np.arange(len(reference)) * BH.SAMPLE_TIME_SECONDS
    fig, axes = plt.subplots(3, 1, figsize=(13, 10), constrained_layout=True)
    colors = ("#c62828", "#1565c0", "#2e7d32")
    for axis, color, label in zip(range(3), colors, ("x", "y", "z")):
        axes[0].plot(time, reference[:, axis], color=color, ls=":", lw=1.2,
                     alpha=0.7, label=f"{label} reference")
        axes[0].plot(time, online["S"][:, axis], color=color, lw=1.0,
                     label=f"{label} online")
    axes[0].axhline(hard_bound, color="black", ls="--", lw=1)
    axes[0].axhline(-hard_bound, color="black", ls="--", lw=1)
    axes[0].set_ylabel("Task output (standardized)")
    axes[0].legend(ncol=3, fontsize=8)

    axes[1].plot(time, online["estimated_sigma_eps_std"], color="#6a1b9a",
                 label="online composite innovation std")
    axes[1].axhline(online["estimated_sigma_eps_std"][0], color="#555555",
                    ls="--", label="offline initialization")
    axes[1].set_ylabel("Latent innovation RMS std")
    axes[1].legend()

    decisions = np.arange(len(online["stage_tightening_history"])) \
        * BH.CONTROL_INTERVAL_STEPS * BH.SAMPLE_TIME_SECONDS
    terminal_tightening = online["stage_tightening_history"][:, -1, :]
    for axis, color, label in zip(range(3), colors, ("x", "y", "z")):
        axes[2].plot(decisions, terminal_tightening[:, axis], color=color,
                     label=f"{label} terminal tightening")
    frozen_terminal = np.asarray(frozen["stage_tightening"])[-1]
    for value, color in zip(frozen_terminal, colors):
        axes[2].axhline(value, color=color, ls="--", alpha=0.55)
    axes[2].set_ylabel("18-step tightening")
    axes[2].set_xlabel("Time (s)")
    axes[2].legend(ncol=3, fontsize=8)
    fig.suptitle("copyBJ: online covariance-adaptive SMPC on the copyBI event")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def run_experiment(output_dir=None, steps=DEFAULT_STEPS):
    steps = int(steps)
    reference = BI.hard_edge_face_petal_reference(steps)
    model, E, source_hard_bound, scales, noise = BH.load_identification_and_noise_objects()
    hard_bound = float(BI.EXPERIMENT_HARD_BOUND)
    covariance = BH.build_process_covariances(
        model, E, scales["y_scale"][:3], BH.SIGMA_WIND_MPS,
    )
    residual = BI.innovation_sequence(model, steps, BI.INNOVATION_SEED)
    wind_velocity = BI.physical_wind_velocity(
        steps, BH.SIGMA_WIND_MPS, BI.WIND_SEED,
    )
    task_map = E.T @ model["P"][:, :model["R"].shape[1]]
    wind_latent, _ = BH.wind_velocity_to_latent_increment(
        wind_velocity, task_map, scales["y_scale"][:3],
    )

    frozen = BH.simulate_controller(
        model, E, hard_bound, reference, residual, wind_latent,
        covariance["Sigma_eps_aug"], BH.ALPHA_JOINT,
        initial_task=np.zeros(3),
    )
    raw_online = simulate_online_controller(
        model, E, hard_bound, reference, residual, wind_latent,
        covariance["Sigma_eps_aug"], noise["Sigma_obs_proxy"],
        initial_task=np.zeros(3), online_weight=1.0,
    )
    online = simulate_online_controller(
        model, E, hard_bound, reference, residual, wind_latent,
        covariance["Sigma_eps_aug"], noise["Sigma_obs_proxy"],
        initial_task=np.zeros(3), online_weight=SHRUNK_ONLINE_WEIGHT,
    )
    hashes = {
        frozen["disturbance_sha256"], raw_online["disturbance_sha256"],
        online["disturbance_sha256"],
    }
    if len(hashes) != 1:
        raise RuntimeError("Frozen and online runs did not share disturbances")
    if steps == DEFAULT_STEPS:
        if raw_online["qp_failure_step"] is None:
            raise RuntimeError("Raw 40-sample online run unexpectedly remained feasible")
        if not str(raw_online["qp_failure_status"]).startswith("empty_chance_interval"):
            raise RuntimeError("Raw online run failed for an unexpected reason")
    for name, result in (("frozen", frozen), ("shrunk online", online)):
        if result["completed_steps"] != steps or result["qp_failure_count"]:
            raise RuntimeError(f"{name} controller did not complete all steps")
        if result["fallback_count"]:
            raise RuntimeError(f"{name} controller used a fallback")
        if result["maximum_qp_constraint_residual"] > 1e-6:
            raise RuntimeError(f"{name} QP residual exceeds 1e-6")

    results = Path(output_dir) if output_dir else HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    npz_path = results / "copyBJ_online_covariance_comparison.npz"
    json_path = results / "copyBJ_online_covariance_summary.json"
    figure_path = results / "copyBJ_online_covariance_comparison.png"
    np.savez_compressed(
        npz_path,
        reference_standardized=reference,
        frozen_trajectory_standardized=frozen["S"],
        raw_online_trajectory_standardized=raw_online["S"],
        online_trajectory_standardized=online["S"],
        frozen_control=frozen["U"],
        online_control=online["U"],
        residual_innovation_latent=residual,
        wind_velocity_mps=wind_velocity,
        wind_latent_increment=wind_latent,
        online_estimated_sigma_eps_std=online["estimated_sigma_eps_std"],
        online_estimated_sigma_obs_std=online["estimated_sigma_obs_std"],
        online_stage_tightening=online["stage_tightening_history"],
        online_final_Sigma_eps=online["final_Sigma_eps"],
        online_final_Sigma_obs_proxy=online["final_Sigma_obs_proxy"],
        raw_online_completed_steps=np.array([raw_online["completed_steps"]]),
        raw_online_failure_step=np.array([raw_online["qp_failure_step"]]),
        raw_online_stage_tightening=raw_online["stage_tightening_history"],
    )
    _save_figure(figure_path, reference, frozen, online, hard_bound)
    summary = {
        "experiment": "copyBJ copyBI event with copyAU-style online covariance learning",
        "steps": steps,
        "window_samples": WINDOW_SAMPLES,
        "update_min_samples": UPDATE_MIN_SAMPLES,
        "fixed_objects": list(FIXED_OBJECTS),
        "online_updated_objects": list(ONLINE_UPDATED_OBJECTS),
        "raw_copyAU_style_online": {
            **_summary(raw_online),
            "online_weight": 1.0,
            "failure_step": raw_online["qp_failure_step"],
            "failure_time_seconds": (
                None if raw_online["qp_failure_step"] is None
                else raw_online["qp_failure_step"] * BH.SAMPLE_TIME_SECONDS
            ),
            "failure_status": raw_online["qp_failure_status"],
        },
        "shrinkage_candidate": {
            "online_weight": SHRUNK_ONLINE_WEIGHT,
            "offline_prior_weight": 1.0 - SHRUNK_ONLINE_WEIGHT,
            "status": "experimental regularization; not part of the copyAU theory",
        },
        "online_covariance_semantics": (
            "composite closed-loop latent innovation; in this simulation it absorbs "
            "both sampled identified innovation and unmodelled injected physical wind"
        ),
        "chance_covariance_semantics": (
            "initially offline Sigma_eps; thereafter the 40-sample online composite "
            "innovation covariance; no Oracle wind covariance is supplied"
        ),
        "shared_disturbance_sha256": online["disturbance_sha256"],
        "frozen_copyBI_smpc": _summary(frozen),
        "online_copyBJ_smpc": _summary(online),
        "initial_sigma_eps_rms_std": float(online["estimated_sigma_eps_std"][0]),
        "final_sigma_eps_rms_std": float(online["estimated_sigma_eps_std"][-1]),
        "minimum_sigma_eps_rms_std": float(np.min(online["estimated_sigma_eps_std"])),
        "maximum_sigma_eps_rms_std": float(np.max(online["estimated_sigma_eps_std"])),
        "claim_boundary": (
            "This is an online covariance-adaptation experiment on the frozen identified "
            "Pelican model, not online learning of E/A/B/P/R, not process/measurement "
            "noise separation, and not real-flight validation."
        ),
        "artifacts": {"data": npz_path.name, "figure": figure_path.name},
    }
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--steps", type=int, default=DEFAULT_STEPS)
    args = parser.parse_args()
    print(json.dumps(run_experiment(args.output_dir, args.steps), indent=2))


if __name__ == "__main__":
    _main()
