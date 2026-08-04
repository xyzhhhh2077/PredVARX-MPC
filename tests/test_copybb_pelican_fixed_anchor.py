import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "experiments" / "copyBB_pelican_fixed_anchor_smpc" / "run_fixed_anchor_closed_loop.py"


def load_module():
    spec = importlib.util.spec_from_file_location("copybb", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_task_cost_has_no_orthogonal_penalty():
    mod = load_module()
    E = np.eye(6)[:, :2]
    Q = mod.task_cost_matrix(E, 80.0)
    complement = np.eye(6) - E @ E.T
    assert np.allclose(Q, 80.0 * E @ E.T)
    assert np.linalg.norm(Q @ complement) < 1e-12


def test_lifted_noise_uses_physical_step_model():
    mod = load_module()
    A = np.diag([0.8, 0.6])
    Sigma = np.diag([0.2, 0.1])
    got = mod.lift_noise_covariance(A, Sigma, 3)
    expected = sum(
        np.linalg.matrix_power(A, i) @ Sigma @ np.linalg.matrix_power(A, i).T
        for i in range(3)
    )
    wrong = sum(
        np.linalg.matrix_power(np.linalg.matrix_power(A, 3), i)
        @ Sigma
        @ np.linalg.matrix_power(np.linalg.matrix_power(A, 3), i).T
        for i in range(3)
    )
    assert np.allclose(got, expected)
    assert not np.allclose(got, wrong)


def test_common_noise_drives_both_comparison_plants():
    mod = load_module()
    model = {
        "A": np.array([[0.9]]),
        "B": np.array([[1.0]]),
        "u_mean": np.array([[0.0]]),
    }
    eps = np.array([[0.25]])
    z0 = np.array([[0.4]])
    z_a = mod.plant_step_with_noise(model, z0, np.array([0.0]), eps)
    z_b = mod.plant_step_with_noise(model, z0, np.array([0.0]), eps)
    assert np.array_equal(z_a, z_b)


def test_fixed_anchor_is_standardized_xy_position():
    mod = load_module()
    E = mod.fixed_position_anchor(10)
    expected = np.zeros((10, 2))
    expected[0, 0] = 1.0
    expected[1, 1] = 1.0
    assert np.array_equal(E, expected)
    assert np.allclose(E.T @ E, np.eye(2))
