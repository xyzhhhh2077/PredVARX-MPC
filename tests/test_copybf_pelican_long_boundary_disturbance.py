import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT / "experiments" / "copyBF_pelican_long_3d_boundary_disturbance"
    / "run_long_3d_boundary_disturbance.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location("copybf", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_reference_keeps_selected_axis_close_to_the_hard_bound():
    m = load_module()
    reference = m.long_boundary_reference(
        hard_bound=3.5, steps=12000, margin=0.02,
        inward_excursion=0.08, period_steps=1200, axis=2,
    )

    assert reference.shape == (12000, 3)
    assert np.all(reference < 3.5)
    assert np.all(reference >= 0.0)
    assert np.allclose(reference[:, :2], 0.0)
    assert reference[:, 2].max() > 3.47
    assert reference[:, 2].min() > 3.39
    assert np.sum(reference[:, 2] > 3.45) > 100
    assert np.allclose(reference[0], reference[-1], atol=1e-12)


def test_external_gusts_are_separate_and_outward_on_all_task_axes():
    m = load_module()
    task_map = np.eye(3)
    gusts, task_gusts = m.external_gust_sequence(
        task_map=task_map, steps=12000, period_steps=1200,
        gust_width=30, gust_amplitude=0.025,
    )

    assert gusts.shape == (12000, 3)
    assert task_gusts.shape == (12000, 3)
    assert np.allclose(gusts, task_gusts)
    assert np.all(task_gusts.max(axis=0) >= 0.025)
    assert np.all(task_gusts.min(axis=0) <= -0.025)
    assert np.all(np.count_nonzero(task_gusts, axis=0) > 0)
    for start in range(0, 12000, 1200):
        assert np.allclose(np.sum(task_gusts[start:start + 1200], axis=0), 0.0)


def test_external_gust_can_be_restricted_to_one_critical_axis():
    m = load_module()
    _, task_gusts = m.external_gust_sequence(
        task_map=np.eye(3), steps=1200, critical_axis=1,
    )

    assert np.allclose(task_gusts[:, 0], 0.0)
    assert np.any(task_gusts[:, 1] > 0.0)
    assert np.any(task_gusts[:, 1] < 0.0)
    assert np.allclose(task_gusts[:, 2], 0.0)


def test_cycle_metrics_report_long_term_error_and_violations():
    m = load_module()
    reference = np.zeros((12, 3))
    trajectory = reference.copy()
    trajectory[7, 2] = 1.1
    metrics = m.cycle_metrics(
        trajectory, reference, hard_bound=1.0, period_steps=4,
    )

    assert len(metrics) == 3
    assert metrics[0]["violation_steps"] == 0
    assert metrics[1]["violation_steps"] == 1
    assert metrics[1]["rmse"][2] == 0.55
    assert metrics[2]["violation_steps"] == 0


def test_equilibrium_initial_state_reproduces_requested_task():
    m = load_module()
    model, E, _, _ = m.BE.load_model()
    target = np.array([0.0, 0.0, 3.5])
    state, control = m.equilibrium_initial_condition(model, E, target)

    realized = (E.T @ (model["y_mean"] + model["P"] @ state)).ravel()
    assert np.allclose(realized, target, atol=1e-6)
    assert np.allclose(
        state,
        model["A"] @ state + model["B"] @ (
            control.reshape(-1, 1) - model["u_mean"]
        ),
        atol=1e-6,
    )
    assert np.all(np.abs(control) <= 6.0 + 1e-8)
