import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "experiments" / "copyBD_pelican_3d_fixed_anchor_smpc" / "run_3d_boundary_trajectory.py"


def load_module():
    spec = importlib.util.spec_from_file_location("copybd", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_fixed_anchor_is_standardized_xyz_position():
    mod = load_module()
    E = mod.fixed_position_anchor_3d(10)
    expected = np.zeros((10, 3))
    expected[:3, :3] = np.eye(3)
    assert np.array_equal(E, expected)
    assert np.allclose(E.T @ E, np.eye(3))


def test_3d_reference_moves_on_all_axes_and_stays_inside_bounds():
    mod = load_module()
    hard_bound = 3.9
    ref = mod.moving_reference_3d(steps=1200, hard_bound=hard_bound)
    assert ref.shape == (1200, 3)
    assert np.allclose(np.ptp(ref, axis=0), [0.5, 1.1, 0.64], atol=0.01)
    assert np.max(np.abs(ref)) < 0.6
    assert np.linalg.norm(ref[0] - ref[-1]) < 0.05
    centered = ref - np.mean(ref, axis=0)
    assert np.linalg.matrix_rank(centered, tol=1e-8) == 3


def test_risk_allocation_counts_three_axes():
    mod = load_module()
    got = mod.per_face_risk(alpha_joint=0.05, task_axes=3, horizon=18)
    assert np.isclose(got, 0.05 / (2 * 3 * 18))


def test_task_cost_only_penalizes_three_anchored_axes():
    mod = load_module()
    E = mod.fixed_position_anchor_3d(10)
    Q = mod.task_cost_matrix(E, 80.0)
    assert np.allclose(Q, 80.0 * E @ E.T)
    assert np.allclose(np.diag(Q)[:3], 80.0)
    assert np.allclose(Q[3:, :], 0.0)


def test_prediction_horizon_uses_future_decision_time_references():
    mod = load_module()
    reference = np.column_stack([
        np.arange(50),
        100 + np.arange(50),
        200 + np.arange(50),
    ])
    got = mod.future_reference_horizon(reference, k=5, d=4, horizon=3)
    expected = reference[[9, 13, 17]]
    assert np.array_equal(got, expected)


def test_prediction_horizon_holds_last_reference_past_trajectory_end():
    mod = load_module()
    reference = np.arange(18).reshape(6, 3)
    got = mod.future_reference_horizon(reference, k=4, d=2, horizon=3)
    assert np.array_equal(got, np.repeat(reference[[-1]], 3, axis=0))
