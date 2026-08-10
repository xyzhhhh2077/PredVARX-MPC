import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
BL_SCRIPT = (
    ROOT / "experiments" / "copyBL_pelican_dual_gated_covariance_smpc"
    / "run_dual_gated_covariance_smpc.py"
)


def load_module(name):
    spec = importlib.util.spec_from_file_location(name, BL_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_feasibility_gate_backtracks_until_candidate_is_accepted():
    bl = load_module("copybl_backtracking_contract")
    accepted = np.eye(2)
    raw = 5.0 * np.eye(2)

    result = bl.backtracking_feasibility_gate(
        accepted,
        raw,
        evaluator=lambda covariance: (
            float(covariance[0, 0]) <= 3.1,
            {"minimum_predicted_slack": 3.1 - float(covariance[0, 0])},
        ),
        factors=(1.0, 0.5, 0.25),
    )

    assert result.accepted
    assert result.factor == 0.5
    assert np.allclose(result.covariance, 3.0 * np.eye(2))
    assert [attempt["factor"] for attempt in result.attempts] == [1.0, 0.5]


def test_feasibility_gate_preserves_last_accepted_covariance_when_all_attempts_fail():
    bl = load_module("copybl_rejection_contract")
    accepted = np.diag([1.0, 2.0])

    result = bl.backtracking_feasibility_gate(
        accepted,
        8.0 * np.eye(2),
        evaluator=lambda covariance: (False, {"reason": "infeasible"}),
        factors=(1.0, 0.5, 0.25),
    )

    assert not result.accepted
    assert result.factor == 0.0
    assert np.array_equal(result.covariance, accepted)
    assert len(result.attempts) == 3


def test_control_value_gate_accepts_only_improvements_beyond_deadband():
    bl = load_module("copybl_value_acceptance_contract")

    accepted = bl.control_value_gate(
        baseline_cost=100.0,
        candidate_cost=98.0,
        deadband=0.01,
        candidate_completed=True,
    )
    rejected = bl.control_value_gate(
        baseline_cost=100.0,
        candidate_cost=99.5,
        deadband=0.01,
        candidate_completed=True,
    )

    assert accepted.accepted
    assert accepted.relative_improvement == 0.02
    assert accepted.reason == "improvement_exceeds_deadband"
    assert not rejected.accepted
    assert rejected.reason == "improvement_below_deadband"


def test_control_value_gate_rejects_failed_candidate_replay():
    bl = load_module("copybl_value_failure_contract")

    result = bl.control_value_gate(
        baseline_cost=100.0,
        candidate_cost=float("inf"),
        deadband=0.01,
        candidate_completed=False,
    )

    assert not result.accepted
    assert result.reason == "candidate_replay_failed"


def test_control_value_gate_rejects_incomplete_baseline_comparison():
    bl = load_module("copybl_baseline_failure_contract")

    result = bl.control_value_gate(
        baseline_cost=100.0,
        candidate_cost=80.0,
        deadband=0.01,
        candidate_completed=True,
        baseline_completed=False,
    )

    assert not result.accepted
    assert result.reason == "baseline_replay_failed"


def test_control_value_gate_deadband_is_strict_and_rejects_nonfinite_completed_costs():
    bl = load_module("copybl_value_deadband_contract")

    boundary = bl.control_value_gate(100.0, 99.0, 0.01, True)
    assert not boundary.accepted
    assert boundary.reason == "improvement_below_deadband"
    with pytest.raises(ValueError):
        bl.control_value_gate(100.0, float("nan"), 0.01, True)


def test_paired_replay_uses_global_future_reference_and_does_not_mutate_start(monkeypatch):
    bl = load_module("copybl_paired_reference_contract")
    seen_horizons = []

    class FakeController:
        N = 2

        def __init__(self, *args, **kwargs):
            pass

        def step(self, state, refs):
            seen_horizons.append(np.asarray(refs).copy())
            return np.array([0.0]), {
                "status": "optimal",
                "maximum_constraint_residual": 0.0,
                "minimum_predicted_slack": 1.0,
            }

    monkeypatch.setattr(bl.BJ, "AdaptiveBoundaryController", FakeController)
    monkeypatch.setattr(
        bl.BI,
        "realized_stage_cost",
        lambda task, reference, control, mean: np.sum((task - reference) ** 2, axis=1),
    )
    model = {
        "A": np.eye(2),
        "B": np.zeros((2, 1)),
        "P": np.zeros((3, 2)),
        "y_mean": np.zeros((3, 1)),
        "u_mean": np.zeros((1, 1)),
    }
    reference = np.arange(120, dtype=float).reshape(40, 3)
    residual = np.zeros((11, 1))
    wind = np.zeros((11, 1))
    state = np.array([[1.0], [2.0]])
    control = np.array([3.0])
    state_before, control_before = state.copy(), control.copy()

    result = bl._paired_replay(
        model, np.eye(3), 1e6, reference, 10, residual, wind,
        state, control, np.eye(1),
    )

    assert result["completed"]
    assert np.array_equal(state, state_before)
    assert np.array_equal(control, control_before)
    expected = [
        bl.BH.future_reference_horizon(reference, global_step, 2,
                                       control_interval=bl.BH.CONTROL_INTERVAL_STEPS)
        for global_step in (10, 20)
    ]
    assert len(seen_horizons) == 2
    assert all(np.array_equal(actual, wanted) for actual, wanted in zip(seen_horizons, expected))


def test_paired_replay_fails_closed_when_candidate_controller_cannot_be_built(monkeypatch):
    bl = load_module("copybl_paired_failure_contract")

    class BrokenController:
        def __init__(self, *args, **kwargs):
            raise ValueError("empty chance interval")

    monkeypatch.setattr(bl.BJ, "AdaptiveBoundaryController", BrokenController)
    model = {
        "A": np.eye(2), "B": np.zeros((2, 1)), "P": np.zeros((3, 2)),
        "y_mean": np.zeros((3, 1)), "u_mean": np.zeros((1, 1)),
    }
    result = bl._paired_replay(
        model, np.eye(3), 3.8, np.zeros((20, 3)), 0,
        np.zeros((10, 1)), np.zeros((10, 1)), np.zeros((2, 1)),
        np.zeros(1), np.eye(1),
    )
    assert not result["completed"]
    assert result["total_cost"] == float("inf")
    assert result["failure_status"].startswith("controller_setup_failed:")


def test_json_safe_replaces_every_nonfinite_number_and_strict_json_rejects_none_never():
    bl = load_module("copybl_json_contract")
    payload = {"a": float("nan"), "b": [float("inf"), np.float64(-np.inf)], "c": 2.0}
    safe = bl._json_safe(payload)
    assert safe == {"a": None, "b": [None, None], "c": 2.0}
    encoded = bl._strict_json_dumps(payload)
    assert "NaN" not in encoded and "Infinity" not in encoded
