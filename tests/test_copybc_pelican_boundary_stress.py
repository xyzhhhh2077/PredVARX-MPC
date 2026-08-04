import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "copyBC_pelican_boundary_stress" / "run_boundary_stress.py"


def load_module():
    spec = importlib.util.spec_from_file_location("copybc", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_critical_reference_is_inside_hard_bound_but_outside_tightened_bound():
    m = load_module()
    hard = 3.9124
    tightened = np.array([3.5750, 3.5892])
    ref = m.critical_reference(hard, margin=0.02, axis=0)

    assert np.all(np.abs(ref) < hard)
    assert ref[0] > tightened[0]
    assert ref[1] == 0.0


def test_common_disturbances_are_identical_for_both_controllers():
    m = load_module()
    covariance = np.diag([0.04, 0.01])
    disturbances = m.sample_disturbances(covariance, steps=20, seed=17)

    deterministic = disturbances.copy()
    stochastic = disturbances.copy()
    assert np.array_equal(deterministic, stochastic)
    assert disturbances.shape == (20, 2, 1)


def test_boundary_metrics_count_actual_violations_and_minimum_margin():
    m = load_module()
    trajectory = np.array([[3.8, 0.0], [3.93, 0.1], [3.7, -4.0]])
    metrics = m.boundary_metrics(trajectory, hard_bound=3.9124)

    assert metrics["violation_steps"] == 2
    assert metrics["minimum_hard_margin"] < 0.0
    assert metrics["maximum_abs_position"] == [3.93, 4.0]


def test_tightening_is_zero_only_for_deterministic_controller():
    m = load_module()
    stochastic = m.controller_tightening(alpha_joint=0.05, axes=2, horizon=18,
                                         standard_deviation=0.1)
    deterministic = m.controller_tightening(alpha_joint=None, axes=2, horizon=18,
                                            standard_deviation=0.1)

    assert stochastic > 0.3
    assert deterministic == 0.0


def test_initial_latent_state_reproduces_requested_task_position():
    m = load_module()
    model, E, _, hard_bound, _ = m.BB.load_model_and_data()
    target = m.critical_reference(hard_bound, margin=0.005, axis=0)
    z0 = m.initial_latent_state(model, E, target)
    realized = (E.T @ (model["y_mean"] + model["P"] @ z0)).ravel()

    assert np.allclose(realized, target, atol=1e-10)


def test_both_boundary_scenarios_target_one_axis_each():
    m = load_module()
    scenarios = m.boundary_scenarios(hard_bound=3.9124, margin=0.005)

    assert list(scenarios) == ["x_critical", "y_critical"]
    assert np.allclose(scenarios["x_critical"], [3.9074, 0.0])
    assert np.allclose(scenarios["y_critical"], [0.0, 3.9074])


def test_moving_boundary_reference_stays_near_x_bound_and_moves_in_y():
    m = load_module()
    reference = m.moving_boundary_reference(
        hard_bound=3.9124, steps=1200, margin=0.005,
        y_amplitude=1.0, period_steps=600,
    )

    assert reference.shape == (1200, 2)
    assert np.allclose(reference[:, 0], 3.9074)
    assert np.isclose(reference[:, 1].max(), 1.0)
    assert np.isclose(reference[:, 1].min(), -1.0)
    assert np.all(np.abs(reference) < 3.9124)
