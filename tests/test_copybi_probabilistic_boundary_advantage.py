import importlib.util
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BI_SCRIPT = (
    ROOT / "experiments" / "copyBI_pelican_probabilistic_boundary_advantage"
    / "run_probabilistic_boundary_advantage.py"
)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_slow_z_hard_edge_reference_has_two_long_dwells_and_smooth_transitions():
    bi = load_module("copybi_reference_contract", BI_SCRIPT)
    steps = 12_000
    edge = 3.88
    reference = bi.slow_z_hard_edge_reference(steps=steps, edge=edge)

    assert reference.shape == (steps, 3)
    assert np.allclose(reference[0], 0.0)
    assert np.allclose(reference[-1], 0.0)
    assert np.max(np.abs(reference[:, :2])) == 0.0
    assert np.isclose(np.max(reference[:, 2]), edge)
    assert np.isclose(np.min(reference[:, 2]), -edge)
    assert np.count_nonzero(reference[:, 2] >= edge - 1e-12) >= 0.19 * steps
    assert np.count_nonzero(reference[:, 2] <= -edge + 1e-12) >= 0.19 * steps
    assert np.max(np.abs(np.diff(reference[:, 2]))) < 0.01


def test_copybi_reuses_frozen_copybh_controller_and_noise_functions():
    bi = load_module("copybi_reuse_contract", BI_SCRIPT)

    assert bi.BH.INPUT_COMMAND_BOUND_STANDARDIZED == 7.0
    assert bi.CachedBoundaryController is bi.BH.CachedBoundaryController
    assert bi.run_controller_pair is bi.BH.run_controller_pair
    assert bi.innovation_sequence is bi.BH.innovation_sequence
    assert bi.physical_wind_velocity is bi.BH.physical_wind_velocity


def test_hard_edge_face_petal_reference_uses_frozen_copybh_geometry():
    bi = load_module("copybi_face_petal_contract", BI_SCRIPT)
    reference = bi.hard_edge_face_petal_reference(steps=18_000)
    expected = bi.BH.face_petal_stress_reference(
        steps=18_000,
        pressures=bi.FACE_PRESSURES,
        slide=bi.BH.FACE_PETAL_SLIDE,
    )

    assert np.array_equal(reference, expected)
    assert np.array_equal(bi.FACE_PRESSURES, np.array([3.46, 3.46, 3.46]))
    assert np.allclose(np.max(np.abs(reference), axis=0), [3.46, 3.46, 3.46])
    assert bi.EXPERIMENT_HARD_BOUND == 3.8
    assert np.all(bi.FACE_PRESSURES < bi.EXPERIMENT_HARD_BOUND)
    _, _, model_hard, _, _ = bi.BH.load_identification_and_noise_objects()
    assert np.all(bi.FACE_PRESSURES < model_hard)


def test_validate_sampled_advantage_requires_shared_noise_and_full_safe_smpc():
    bi = load_module("copybi_validation_contract", BI_SCRIPT)
    smpc = {
        "completed_steps": 100,
        "violation_steps": 0,
        "disturbance_sha256": "same",
        "fallback_count": 0,
        "qp_failure_step": None,
    }
    mpc = {
        "completed_steps": 100,
        "violation_steps": 3,
        "disturbance_sha256": "same",
        "fallback_count": 0,
        "qp_failure_step": None,
    }

    bi.validate_sampled_advantage(smpc, mpc, expected_steps=100)

    mpc["disturbance_sha256"] = "different"
    with np.testing.assert_raises_regex(RuntimeError, "identical disturbance"):
        bi.validate_sampled_advantage(smpc, mpc, expected_steps=100)


def test_realized_stage_cost_matches_tracking_and_control_penalties():
    bi = load_module("copybi_stage_cost_contract", BI_SCRIPT)
    task = np.array([[1.0, 2.0, 3.0], [2.0, 2.0, 2.0]])
    reference = np.array([[0.0, 2.0, 4.0], [1.0, 2.0, 3.0]])
    control = np.array([[2.0, 1.0], [1.0, 3.0]])
    control_mean = np.array([1.0, 1.0])

    cost = bi.realized_stage_cost(
        task, reference, control, control_mean, q_weight=2.0, ru=0.5,
    )

    assert np.allclose(cost, [4.5, 6.0])
