"""copyBL: dual-gated online covariance adaptation for the Pelican SMPC task.

Candidate covariance updates first pass a backtracking feasibility gate and then
an offline paired closed-loop control-value proxy. Only candidates passing both
gates atomically replace the last accepted online covariance.
"""

import argparse
from dataclasses import dataclass
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from typing import Callable, Iterable

import numpy as np


HERE = Path(__file__).resolve().parent
BJ_SCRIPT = (
    HERE.parent / "copyBJ_pelican_online_covariance_smpc"
    / "run_online_covariance_smpc.py"
)
_spec = importlib.util.spec_from_file_location("copybl_copybj", BJ_SCRIPT)
BJ = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = BJ
_spec.loader.exec_module(BJ)
BI = BJ.BI
BH = BJ.BH

WINDOW_SAMPLES = 40
VALUE_EVALUATION_INTERVAL_STEPS = 400
VALUE_REPLAY_STEPS = 180
VALUE_DEADBAND = 0.01
BACKTRACKING_FACTORS = (1.0, 0.8, 0.64, 0.512, 0.4096, 0.32768, 0.262144)
DEFAULT_STEPS = BI.DEFAULT_STEPS
ONLINE_UPDATED_OBJECTS = ("Sigma_eps",)
FIXED_OBJECTS = ("E", "A", "B", "P", "R")


def _json_safe(value):
    """Recursively convert NumPy values and non-finite floats to strict JSON values."""
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if np.isfinite(number) else None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def _strict_json_dumps(value):
    return json.dumps(_json_safe(value), indent=2, allow_nan=False)


@dataclass(frozen=True)
class FeasibilityGateResult:
    accepted: bool
    covariance: np.ndarray
    factor: float
    attempts: tuple


@dataclass(frozen=True)
class ControlValueGateResult:
    accepted: bool
    relative_improvement: float
    reason: str


def control_value_gate(
    baseline_cost: float,
    candidate_cost: float,
    deadband: float,
    candidate_completed: bool,
    baseline_completed: bool = True,
) -> ControlValueGateResult:
    """Apply the Lin-inspired paired closed-loop cost deadband."""
    baseline = float(baseline_cost)
    candidate = float(candidate_cost)
    threshold = float(deadband)
    if not 0.0 <= threshold < 1.0:
        raise ValueError("deadband must lie in [0, 1)")
    if not baseline_completed:
        return ControlValueGateResult(False, float("nan"), "baseline_replay_failed")
    if not candidate_completed:
        return ControlValueGateResult(False, float("nan"), "candidate_replay_failed")
    if not np.isfinite(baseline) or baseline <= 0.0 or not np.isfinite(candidate):
        raise ValueError("completed replay costs must be finite and baseline cost positive")
    improvement = (baseline - candidate) / baseline
    if improvement > threshold:
        return ControlValueGateResult(
            True, float(improvement), "improvement_exceeds_deadband",
        )
    return ControlValueGateResult(
        False, float(improvement), "improvement_below_deadband",
    )


def backtracking_feasibility_gate(
    accepted_covariance: np.ndarray,
    raw_candidate: np.ndarray,
    evaluator: Callable[[np.ndarray], tuple[bool, dict]],
    factors: Iterable[float],
) -> FeasibilityGateResult:
    """Return the first feasible convex step from accepted to raw candidate."""
    accepted = np.asarray(accepted_covariance, dtype=float)
    raw = np.asarray(raw_candidate, dtype=float)
    if accepted.shape != raw.shape:
        raise ValueError("accepted and candidate covariances must share a shape")

    attempts = []
    for factor in factors:
        factor = float(factor)
        if not 0.0 < factor <= 1.0:
            raise ValueError("backtracking factors must lie in (0, 1]")
        candidate = (1.0 - factor) * accepted + factor * raw
        feasible, details = evaluator(candidate)
        attempts.append({"factor": factor, "feasible": bool(feasible), **details})
        if feasible:
            return FeasibilityGateResult(
                accepted=True,
                covariance=candidate.copy(),
                factor=factor,
                attempts=tuple(attempts),
            )

    return FeasibilityGateResult(
        accepted=False,
        covariance=accepted.copy(),
        factor=0.0,
        attempts=tuple(attempts),
    )


def _augment_covariance(latent_covariance: np.ndarray, state_dimension: int) -> np.ndarray:
    latent = np.asarray(latent_covariance, dtype=float)
    augmented = np.zeros((state_dimension, state_dimension), dtype=float)
    augmented[: latent.shape[0], : latent.shape[1]] = latent
    return augmented


def _paired_replay(
    model: dict,
    E: np.ndarray,
    hard_bound: float,
    full_reference: np.ndarray,
    start_step: int,
    residual: np.ndarray,
    wind: np.ndarray,
    initial_state: np.ndarray,
    initial_control: np.ndarray,
    latent_covariance: np.ndarray,
) -> dict:
    """Replay one global future segment from a copied state and control.

    Keeping the full reference is essential: QPs near the end of the cost segment
    still see the real future reference, rather than a falsely clamped segment end.
    """
    state = np.asarray(initial_state, dtype=float).copy()
    control = np.asarray(initial_control, dtype=float).copy()
    full_reference = np.asarray(full_reference, dtype=float)
    residual = np.asarray(residual, dtype=float)
    wind = np.asarray(wind, dtype=float)
    start_step = int(start_step)
    replay_steps = len(residual)
    if len(wind) != replay_steps or start_step < 0 or start_step + replay_steps > len(full_reference):
        raise ValueError("paired replay arrays do not define one valid global future segment")
    try:
        controller = BJ.AdaptiveBoundaryController(
            model, E, hard_bound,
            _augment_covariance(latent_covariance, model["A"].shape[0]),
        )
    except Exception as exc:
        return {
            "completed": False,
            "total_cost": float("inf"),
            "violation_steps": 0,
            "failure_status": f"controller_setup_failed: {type(exc).__name__}: {exc}",
        }
    tasks = []
    controls = []
    completed = True
    failure_status = ""
    for local_step in range(replay_steps):
        global_step = start_step + local_step
        output = model["y_mean"] + model["P"] @ state
        task = (E.T @ output).ravel()
        if local_step % BH.CONTROL_INTERVAL_STEPS == 0:
            refs = BH.future_reference_horizon(
                full_reference, global_step, controller.N,
                control_interval=BH.CONTROL_INTERVAL_STEPS,
            )
            try:
                candidate, info = controller.step(state, refs)
            except Exception as exc:
                completed = False
                failure_status = f"controller_step_failed: {type(exc).__name__}: {exc}"
                break
            if candidate is None:
                completed = False
                failure_status = str(info["status"])
                break
            control = candidate
        tasks.append(task)
        controls.append(control.copy())
        state = (
            model["A"] @ state
            + model["B"] @ (control.reshape(-1, 1) - model["u_mean"])
        )
        state[: residual.shape[1], 0] += residual[local_step]
        state[: wind.shape[1], 0] += wind[local_step]

    task_array = np.asarray(tasks, dtype=float).reshape(-1, 3)
    control_array = np.asarray(controls, dtype=float).reshape(-1, model["B"].shape[1])
    if len(task_array):
        stage_cost = BI.realized_stage_cost(
            task_array,
            full_reference[start_step : start_step + len(task_array)],
            control_array,
            model["u_mean"].ravel(),
        )
        total_cost = float(np.sum(stage_cost))
        violation_steps = int(np.sum(np.any(np.abs(task_array) > hard_bound, axis=1)))
    else:
        total_cost = float("inf")
        violation_steps = 0
    return {
        "completed": bool(completed and len(task_array) == replay_steps),
        "total_cost": total_cost,
        "violation_steps": violation_steps,
        "failure_status": failure_status,
    }


def simulate_dual_gated_controller(
    model: dict,
    E: np.ndarray,
    hard_bound: float,
    reference: np.ndarray,
    residual_innovation: np.ndarray,
    wind_latent_increment: np.ndarray,
    initial_covariance: np.ndarray,
    initial_task=None,
    value_interval_steps: int = VALUE_EVALUATION_INTERVAL_STEPS,
    value_replay_steps: int = VALUE_REPLAY_STEPS,
    value_deadband: float = VALUE_DEADBAND,
) -> dict:
    """Run transactional feasibility and offline control-value gates."""
    reference = np.asarray(reference, dtype=float)
    residual_innovation = np.asarray(residual_innovation, dtype=float)
    wind_latent_increment = np.asarray(wind_latent_increment, dtype=float)
    steps = len(reference)
    ell = model["A"].shape[0] // 2
    accepted_covariance = np.asarray(initial_covariance[:ell, :ell], dtype=float).copy()
    estimator = BJ.OnlineInnovationCovariance(
        ell, WINDOW_SAMPLES, accepted_covariance,
    )
    controller = BJ.AdaptiveBoundaryController(
        model, E, hard_bound, initial_covariance,
    )
    probe = BJ.AdaptiveBoundaryController(
        model, E, hard_bound, initial_covariance,
    )
    if initial_task is None:
        initial_task = reference[0]
    state, control = BH.equilibrium_initial_condition(model, E, initial_task)
    previous_state = None
    previous_control = None

    task_history = []
    control_history = []
    estimated_std = []
    accepted_std = []
    tightening_history = []
    slacks = []
    duals = []
    residuals = []
    solve_times = []
    events = []
    qp_status_counts = {}
    qp_failure_step = None
    qp_failure_status = ""

    for step in range(steps):
        output = model["y_mean"] + model["P"] @ state
        task = (E.T @ output).ravel()
        if previous_state is not None:
            predicted = (
                model["A"] @ previous_state
                + model["B"] @ (previous_control.reshape(-1, 1) - model["u_mean"])
            )
            estimator.update(state[:ell, 0] - predicted[:ell, 0])
        estimated_std.append(float(np.trace(estimator.covariance) / ell) ** 0.5)

        if step % BH.CONTROL_INTERVAL_STEPS == 0:
            refs = BH.future_reference_horizon(
                reference,
                step,
                controller.N,
                control_interval=BH.CONTROL_INTERVAL_STEPS,
            )

            def evaluate(candidate_covariance):
                try:
                    probe.update_process_covariance(
                        _augment_covariance(candidate_covariance, model["A"].shape[0])
                    )
                except ValueError as exc:
                    return False, {"reason": "empty_chance_interval", "detail": str(exc)}
                candidate_control, candidate_info = probe.step(state, refs)
                feasible = (
                    candidate_control is not None
                    and candidate_info.get("maximum_constraint_residual", float("inf")) <= 1e-6
                    and candidate_info.get("minimum_predicted_slack", -float("inf")) >= -1e-6
                )
                return feasible, {
                    "reason": "feasible" if feasible else str(candidate_info["status"]),
                    "minimum_predicted_slack": float(
                        candidate_info.get("minimum_predicted_slack", float("nan"))
                    ),
                }

            feasibility = backtracking_feasibility_gate(
                accepted_covariance,
                estimator.covariance,
                evaluate,
                BACKTRACKING_FACTORS,
            )
            value_triggered = (
                step >= WINDOW_SAMPLES
                and step % int(value_interval_steps) == 0
                and steps - step >= int(value_replay_steps)
                and feasibility.accepted
                and not np.allclose(feasibility.covariance, accepted_covariance)
            )
            event = {
                "step": int(step),
                "time_seconds": float(step * BH.SAMPLE_TIME_SECONDS),
                "feasibility_triggered": True,
                "feasibility_accepted": bool(feasibility.accepted),
                "backtracking_factor": float(feasibility.factor),
                "backtracking_attempts": [dict(item) for item in feasibility.attempts],
                "value_triggered": bool(value_triggered),
                "value_accepted": False,
                "value_reason": "not_scheduled",
                "relative_cost_improvement": None,
                "update_committed": False,
            }
            if not feasibility.accepted:
                event["value_reason"] = "feasibility_rejected"
            elif value_triggered:
                stop = step + int(value_replay_steps)
                replay_args = (
                    model,
                    E,
                    hard_bound,
                    reference,
                    step,
                    residual_innovation[step:stop],
                    wind_latent_increment[step:stop],
                    state,
                    control,
                )
                baseline_replay = _paired_replay(
                    *replay_args, accepted_covariance,
                )
                candidate_replay = _paired_replay(
                    *replay_args, feasibility.covariance,
                )
                candidate_completed = (
                    candidate_replay["completed"]
                    and candidate_replay["violation_steps"] == 0
                )
                baseline_completed = (
                    baseline_replay["completed"]
                    and baseline_replay["violation_steps"] == 0
                )
                value = control_value_gate(
                    baseline_replay["total_cost"],
                    candidate_replay["total_cost"],
                    value_deadband,
                    candidate_completed,
                    baseline_completed,
                )
                event.update({
                    "value_accepted": bool(value.accepted),
                    "value_reason": value.reason,
                    "relative_cost_improvement": (
                        float(value.relative_improvement)
                        if np.isfinite(value.relative_improvement) else None
                    ),
                    "baseline_replay_cost": _json_safe(baseline_replay["total_cost"]),
                    "baseline_replay_completed": bool(baseline_replay["completed"]),
                    "baseline_replay_violation_steps": int(
                        baseline_replay["violation_steps"]
                    ),
                    "baseline_replay_failure_status": baseline_replay["failure_status"],
                    "candidate_replay_cost": _json_safe(candidate_replay["total_cost"]),
                    "candidate_replay_completed": bool(candidate_replay["completed"]),
                    "candidate_replay_violation_steps": int(
                        candidate_replay["violation_steps"]
                    ),
                    "candidate_replay_failure_status": candidate_replay["failure_status"],
                })
                if value.accepted:
                    accepted_covariance = feasibility.covariance.copy()
                    controller.update_process_covariance(
                        _augment_covariance(accepted_covariance, model["A"].shape[0])
                    )
                    event["update_committed"] = True
            events.append(event)

            candidate, info = controller.step(state, refs)
            status = str(info["status"])
            qp_status_counts[status] = qp_status_counts.get(status, 0) + 1
            if candidate is None:
                qp_failure_step = int(step)
                qp_failure_status = status
                break
            control = candidate
            tightening_history.append(controller.stage_tightening.copy())
            slacks.append(float(info["minimum_predicted_slack"]))
            duals.append(float(info["maximum_constraint_dual"]))
            residuals.append(float(info["maximum_constraint_residual"]))
            solve_times.append(float(info["solve_time_seconds"]))

        accepted_std.append(float(np.trace(accepted_covariance) / ell) ** 0.5)
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
    slack_array = np.asarray(slacks, dtype=float)
    dual_array = np.asarray(duals, dtype=float)
    hard_margin = float(hard_bound) - np.abs(S)
    stage_cost = BI.realized_stage_cost(
        S, reference[:completed], U, model["u_mean"].ravel(),
    )
    feasibility_rejections = sum(not event["feasibility_accepted"] for event in events)
    value_events = [event for event in events if event["value_triggered"]]
    commits = [event for event in events if event["update_committed"]]
    return {
        "S": S,
        "U": U,
        "stage_cost": stage_cost,
        "completed_steps": completed,
        "qp_failure_count": int(qp_failure_step is not None),
        "qp_failure_step": qp_failure_step,
        "qp_failure_status": qp_failure_status,
        "qp_status_counts": qp_status_counts,
        "fallback_count": 0,
        "violation_steps": int(np.sum(np.any(np.abs(S) > hard_bound, axis=1))),
        "minimum_hard_margin": float(np.min(hard_margin)) if completed else float("nan"),
        "rmse": np.sqrt(np.mean((S - reference[:completed]) ** 2, axis=0)),
        "active_qp_steps": int(np.sum((slack_array <= 5e-4) | (dual_array > 1e-6))),
        "maximum_qp_constraint_residual": float(np.max(residuals)),
        "mean_qp_solve_time_seconds": float(np.mean(solve_times)),
        "estimated_sigma_eps_std": np.asarray(estimated_std, dtype=float),
        "accepted_sigma_eps_std": np.asarray(accepted_std, dtype=float),
        "stage_tightening_history": np.asarray(tightening_history, dtype=float),
        "final_accepted_Sigma_eps": accepted_covariance,
        "events": events,
        "feasibility_trigger_count": len(events),
        "feasibility_rejection_count": int(feasibility_rejections),
        "value_trigger_count": len(value_events),
        "value_accept_count": int(sum(event["value_accepted"] for event in value_events)),
        "value_rejection_count": int(sum(not event["value_accepted"] for event in value_events)),
        "commit_count": len(commits),
        "disturbance_sha256": BH._disturbance_sha256(
            residual_innovation, wind_latent_increment,
        ),
    }


def _summary(result: dict) -> dict:
    return {
        "completed_steps": int(result["completed_steps"]),
        "qp_failure_count": int(result["qp_failure_count"]),
        "fallback_count": int(result["fallback_count"]),
        "hard_violation_steps": int(result["violation_steps"]),
        "minimum_hard_margin_standardized": float(result["minimum_hard_margin"]),
        "rmse_standardized": np.asarray(result["rmse"]).tolist(),
        "active_qp_steps": int(result["active_qp_steps"]),
        "maximum_qp_constraint_residual": float(result["maximum_qp_constraint_residual"]),
        "mean_qp_solve_time_seconds": float(result["mean_qp_solve_time_seconds"]),
        "feasibility_trigger_count": int(result["feasibility_trigger_count"]),
        "feasibility_rejection_count": int(result["feasibility_rejection_count"]),
        "value_trigger_count": int(result["value_trigger_count"]),
        "value_accept_count": int(result["value_accept_count"]),
        "value_rejection_count": int(result["value_rejection_count"]),
        "committed_update_count": int(result["commit_count"]),
    }


def run_experiment(output_dir=None, steps=DEFAULT_STEPS) -> dict:
    """Run the fixed copyBI event and persist the copyBL evidence bundle."""
    steps = int(steps)
    reference = BI.hard_edge_face_petal_reference(steps)
    model, E, source_hard_bound, scales, _noise = BH.load_identification_and_noise_objects()
    hard_bound = float(BI.EXPERIMENT_HARD_BOUND)
    covariance = BH.build_process_covariances(
        model, E, scales["y_scale"][:3], BH.SIGMA_WIND_MPS,
    )
    residual = BI.innovation_sequence(model, steps, BI.INNOVATION_SEED)
    wind_velocity = BI.physical_wind_velocity(steps, BH.SIGMA_WIND_MPS, BI.WIND_SEED)
    task_map = E.T @ model["P"][:, : model["R"].shape[1]]
    wind_latent, _ = BH.wind_velocity_to_latent_increment(
        wind_velocity, task_map, scales["y_scale"][:3],
    )
    result = simulate_dual_gated_controller(
        model,
        E,
        hard_bound,
        reference,
        residual,
        wind_latent,
        covariance["Sigma_eps_aug"],
        initial_task=np.zeros(3),
    )
    if result["completed_steps"] != steps or result["qp_failure_count"]:
        raise RuntimeError("copyBL did not complete the requested trajectory")
    if result["fallback_count"] or result["violation_steps"]:
        raise RuntimeError("copyBL used a fallback or crossed the hard bound")
    if result["maximum_qp_constraint_residual"] > 1e-6:
        raise RuntimeError("copyBL QP residual exceeds 1e-6")

    results = Path(output_dir) if output_dir else HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    npz_path = results / "copyBL_dual_gated_run.npz"
    events_path = results / "copyBL_dual_gate_events.json"
    summary_path = results / "copyBL_dual_gated_summary.json"
    np.savez_compressed(
        npz_path,
        reference_standardized=reference,
        trajectory_standardized=result["S"],
        control=result["U"],
        residual_innovation_latent=residual,
        wind_velocity_mps=wind_velocity,
        wind_latent_increment=wind_latent,
        stage_cost=result["stage_cost"],
        estimated_sigma_eps_std=result["estimated_sigma_eps_std"],
        accepted_sigma_eps_std=result["accepted_sigma_eps_std"],
        stage_tightening=result["stage_tightening_history"],
        final_accepted_Sigma_eps=result["final_accepted_Sigma_eps"],
    )
    events_path.write_text(_strict_json_dumps(result["events"]), encoding="utf-8")
    event_sha256 = hashlib.sha256(events_path.read_bytes()).hexdigest()
    summary = {
        "experiment": "copyBL dual-gated online composite-innovation covariance SMPC",
        "steps": steps,
        "sample_time_seconds": BH.SAMPLE_TIME_SECONDS,
        "hard_bound_standardized": hard_bound,
        "source_hard_bound_standardized": float(source_hard_bound),
        "window_samples": WINDOW_SAMPLES,
        "fixed_objects": list(FIXED_OBJECTS),
        "online_updated_objects": list(ONLINE_UPDATED_OBJECTS),
        "backtracking_factors": list(BACKTRACKING_FACTORS),
        "value_evaluation_interval_steps": VALUE_EVALUATION_INTERVAL_STEPS,
        "value_replay_steps": VALUE_REPLAY_STEPS,
        "value_deadband_fraction": VALUE_DEADBAND,
        "shared_disturbance_sha256": result["disturbance_sha256"],
        "event_log_sha256": event_sha256,
        "copyBL": _summary(result),
        "control_value_semantics": (
            "Lin-inspired paired closed-loop realized-cost proxy; not Lin et al.'s "
            "strict prediction power and not an optimal-policy comparison"
        ),
        "offline_replay_scope": (
            "The value gate uses the fixed next 180 simulated samples in this offline "
            "feasibility experiment; it is not a deployable online oracle."
        ),
        "covariance_semantics": (
            "composite closed-loop latent innovation covariance; not separately "
            "identified process, measurement, or physical-wind covariance"
        ),
        "claim_boundary": (
            "Single-seed model-in-the-loop evidence only; no recursive feasibility, "
            "stability, probability calibration, or flight-validation claim."
        ),
        "artifacts": {
            "data": npz_path.name,
            "events": events_path.name,
            "summary": summary_path.name,
        },
    }
    summary_path.write_text(_strict_json_dumps(summary), encoding="utf-8")
    return summary


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--steps", type=int, default=DEFAULT_STEPS)
    args = parser.parse_args()
    print(_strict_json_dumps(run_experiment(args.output_dir, args.steps)))


if __name__ == "__main__":
    _main()