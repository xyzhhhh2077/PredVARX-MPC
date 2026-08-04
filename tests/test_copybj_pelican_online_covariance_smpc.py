import importlib.util
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BJ_SCRIPT = (
    ROOT / "experiments" / "copyBJ_pelican_online_covariance_smpc"
    / "run_online_covariance_smpc.py"
)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_online_covariance_window_discards_samples_older_than_40_steps():
    bj = load_module("copybj_window_contract", BJ_SCRIPT)
    estimator = bj.OnlineInnovationCovariance(
        dimension=2, window=40, initial_covariance=np.eye(2), floor=1e-8,
    )
    for value in range(45):
        estimator.update(np.array([value, -0.5 * value], dtype=float))

    retained = np.arange(5.0, 45.0)
    expected_samples = np.column_stack((retained, -0.5 * retained))
    expected = np.cov(expected_samples, rowvar=False, ddof=1) + 1e-8 * np.eye(2)

    assert estimator.count == 40
    assert np.allclose(estimator.covariance, expected)


def test_online_covariance_is_psd_and_keeps_offline_prior_until_five_samples():
    bj = load_module("copybj_psd_contract", BJ_SCRIPT)
    prior = np.diag([0.3, 0.2, 0.1])
    estimator = bj.OnlineInnovationCovariance(3, 40, prior, floor=1e-8)

    for value in range(4):
        estimator.update(np.array([value, 2.0 * value, -value], dtype=float))
    assert np.array_equal(estimator.covariance, prior)

    estimator.update(np.array([4.0, 8.0, -4.0]))
    assert np.min(np.linalg.eigvalsh(estimator.covariance)) >= 0.999e-8
    assert not np.allclose(estimator.covariance, prior)


def test_parameterized_controller_recomputes_tightening_without_rebuilding_qp():
    bj = load_module("copybj_controller_contract", BJ_SCRIPT)
    model, E, hard_bound, _, _ = bj.BH.load_identification_and_noise_objects()
    initial = np.asarray(model["Sigma_eps"], dtype=float)
    controller = bj.AdaptiveBoundaryController(
        model, E, hard_bound, initial,
        horizon=4, control_interval=2,
    )
    problem_id = id(controller.problem)
    tightening_before = controller.stage_tightening.copy()

    controller.update_process_covariance(1.8 * initial)

    assert id(controller.problem) == problem_id
    assert np.all(controller.stage_tightening > tightening_before)
    assert np.allclose(
        controller.stage_bound_parameter.value,
        controller.stage_mean_bound,
    )


def test_online_experiment_keeps_geometry_fixed_and_uses_same_disturbance():
    bj = load_module("copybj_fairness_contract", BJ_SCRIPT)

    assert bj.WINDOW_SAMPLES == 40
    assert bj.UPDATE_MIN_SAMPLES == 5
    assert bj.ONLINE_UPDATED_OBJECTS == ("Sigma_eps", "Sigma_obs_proxy")
    assert bj.FIXED_OBJECTS == ("E", "A", "B", "P", "R")


def test_shrunk_online_covariance_is_explicit_convex_combination():
    bj = load_module("copybj_shrinkage_contract", BJ_SCRIPT)
    prior = np.diag([1.0, 2.0])
    online = np.array([[4.0, 0.5], [0.5, 6.0]])

    effective = bj.shrink_online_covariance(prior, online, online_weight=0.8)

    assert np.allclose(effective, 0.2 * prior + 0.8 * online)
    with np.testing.assert_raises_regex(ValueError, "online_weight"):
        bj.shrink_online_covariance(prior, online, online_weight=1.01)
