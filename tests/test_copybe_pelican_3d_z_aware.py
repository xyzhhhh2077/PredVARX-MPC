import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "experiments" / "copyBE_pelican_3d_z_aware_smpc" / "run_3d_z_aware_trajectory.py"


def load_module():
    spec = importlib.util.spec_from_file_location("copybe", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_collective_motor_direction_is_unit_common_mode():
    mod = load_module()
    collective = mod.collective_motor_direction(4)
    assert collective.shape == (4, 1)
    assert np.allclose(collective, np.ones((4, 1)) / 2.0)
    assert np.isclose((collective.T @ collective).item(), 1.0)


def test_second_order_companion_contract():
    mod = load_module()
    A1 = np.diag([0.8, 0.7])
    A2 = np.diag([0.1, 0.2])
    B = np.arange(8, dtype=float).reshape(2, 4)
    A_aug, B_aug = mod.second_order_companion(A1, A2, B)
    assert np.array_equal(A_aug[:2], np.hstack([A1, A2]))
    assert np.array_equal(A_aug[2:], np.hstack([np.eye(2), np.zeros((2, 2))]))
    assert np.array_equal(B_aug[:2], B)
    assert np.array_equal(B_aug[2:], np.zeros((2, 4)))


def test_augmented_output_uses_current_latent_state_only():
    mod = load_module()
    P = np.arange(15, dtype=float).reshape(3, 5)
    got = mod.augment_output_loading(P)
    assert got.shape == (3, 10)
    assert np.array_equal(got[:, :5], P)
    assert np.array_equal(got[:, 5:], np.zeros((3, 5)))


def test_xyz_anchor_remains_exact_in_augmented_output_map():
    mod = load_module()
    E = mod.base.fixed_position_anchor_3d(10)
    P = np.column_stack([E, np.zeros((10, 2))])
    P_aug = mod.augment_output_loading(P)
    R = P.copy()
    extractor = np.vstack([R.T, np.zeros((5, 10))])
    assert np.allclose(P_aug @ extractor @ E, E)
