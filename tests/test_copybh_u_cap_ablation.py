import ast
import importlib.util
import inspect
import os
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BG_SCRIPT = (
    ROOT / "experiments" / "copyBG_pelican_physical_wind_smpc"
    / "run_physical_wind_smpc.py"
)
BH_SCRIPT = (
    ROOT / "experiments" / "copyBH_pelican_u_cap_ablation"
    / "run_u_cap_ablation.py"
)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def definition_ast(module, name):
    source = inspect.getsource(getattr(module, name))
    return ast.dump(ast.parse(source), include_attributes=False)


def test_copybh_changes_configuration_not_controller_or_noise_algorithms(monkeypatch):
    monkeypatch.delenv("COPYBH_INPUT_CAP", raising=False)
    bg = load_module("copybg_for_copybh_contract", BG_SCRIPT)
    bh = load_module("copybh_contract", BH_SCRIPT)

    assert bg.INPUT_COMMAND_BOUND_STANDARDIZED == 6.0
    assert bh.INPUT_COMMAND_BOUND_STANDARDIZED == 7.0
    assert bh.GROWING_DC_U_INF_CAP == bh.INPUT_COMMAND_BOUND_STANDARDIZED

    fixed_configuration = (
        "SAMPLE_TIME_SECONDS",
        "CONTROL_INTERVAL_STEPS",
        "HORIZON_STEPS",
        "Q_WEIGHT",
        "RU",
        "ALPHA_JOINT",
        "SIGMA_WIND_MPS",
        "WIND_SEED",
        "INNOVATION_SEED",
    )
    for name in fixed_configuration:
        assert getattr(bh, name) == getattr(bg, name)

    unchanged_definitions = (
        "mstep_lift",
        "lift_noise_covariance",
        "build_process_covariances",
        "physical_wind_velocity",
        "wind_velocity_to_latent_increment",
        "innovation_sequence",
        "CachedBoundaryController",
        "run_controller_pair",
    )
    for name in unchanged_definitions:
        assert definition_ast(bh, name) == definition_ast(bg, name), name


def test_copybh_controller_uses_selected_cap_without_other_qp_changes(monkeypatch):
    monkeypatch.delenv("COPYBH_INPUT_CAP", raising=False)
    bh = load_module("copybh_cap_contract", BH_SCRIPT)
    model, E, hard_bound, scales, _noise = bh.load_identification_and_noise_objects()
    covariance = bh.build_process_covariances(
        model, E, scales["y_scale"][:3], sigma_wind_mps=bh.SIGMA_WIND_MPS,
    )
    controller = bh.CachedBoundaryController(
        model, E, hard_bound, covariance["Sigma_eps_aug"], horizon=3,
        alpha_joint=bh.ALPHA_JOINT,
    )

    assert np.array_equal(controller.u_min, np.full(4, -7.0))
    assert np.array_equal(controller.u_max, np.full(4, 7.0))
    assert controller.N == 3
    assert controller.d == bh.CONTROL_INTERVAL_STEPS
    assert controller.problem.is_dpp()


def test_copybh_cap_override_is_restricted_to_declared_sweep_values():
    old = os.environ.get("COPYBH_INPUT_CAP")
    try:
        for value in ("6.0", "6.5", "7.0"):
            os.environ["COPYBH_INPUT_CAP"] = value
            bh = load_module(f"copybh_cap_{value.replace('.', '_')}", BH_SCRIPT)
            assert bh.INPUT_COMMAND_BOUND_STANDARDIZED == float(value)
        os.environ["COPYBH_INPUT_CAP"] = "7.5"
        try:
            load_module("copybh_bad_cap", BH_SCRIPT)
        except ValueError as exc:
            assert "6.0, 6.5, or 7.0" in str(exc)
        else:
            raise AssertionError("undeclared input cap must be rejected")
    finally:
        if old is None:
            os.environ.pop("COPYBH_INPUT_CAP", None)
        else:
            os.environ["COPYBH_INPUT_CAP"] = old


def test_xy_face_dwell_reference_is_long_safe_and_center_initialized(monkeypatch):
    monkeypatch.delenv("COPYBH_INPUT_CAP", raising=False)
    bh = load_module("copybh_xy_face_dwell", BH_SCRIPT)
    reference = bh.xy_face_dwell_reference(
        steps=12_000,
        x_pressure=3.57,
        y_pressure=3.57,
    )
    near_hard = 3.9124117280376636 - 0.4

    assert "xy-face-dwell" in bh.REFERENCE_MODES
    assert reference.shape == (12_000, 3)
    assert np.allclose(reference[0], 0.0)
    assert np.allclose(reference[-1], 0.0)
    assert np.max(np.abs(reference[:, 2])) == 0.0
    assert np.max(np.abs(reference[:, :2])) <= 3.57 + 1e-12
    assert np.count_nonzero(reference[:, 0] >= near_hard) >= 2_300
    assert np.count_nonzero(reference[:, 1] >= near_hard) >= 2_300
    assert not np.any(
        (reference[:, 0] >= near_hard) & (reference[:, 1] >= near_hard)
    )
    assert np.max(np.linalg.norm(np.diff(reference, axis=0), axis=1)) < 0.01


def test_face_petal_stress_reference_covers_six_faces(monkeypatch):
    monkeypatch.delenv("COPYBH_INPUT_CAP", raising=False)
    bh = load_module("copybh_face_petal", BH_SCRIPT)
    steps = 1_800
    reference = bh.face_petal_stress_reference(
        steps=steps,
        pressures=bh.FACE_PETAL_PRESSURES,
        slide=bh.FACE_PETAL_SLIDE,
    )
    px, py, pz = bh.FACE_PETAL_PRESSURES
    hard = 3.9124117280376636
    mean_x, mean_y, mean_z = 2.935, 3.149, 3.583

    assert "face-petal-stress" in bh.REFERENCE_MODES
    assert reference.shape == (steps, 3)
    assert np.allclose(reference[0], 0.0, atol=1e-12)
    assert np.allclose(reference[-1], 0.0, atol=1e-9)
    assert np.max(np.abs(reference)) < hard
    # All six face directions are visited near the design pressures.
    assert np.max(reference[:, 0]) >= 0.95 * px
    assert np.min(reference[:, 0]) <= -0.95 * px
    assert np.max(reference[:, 1]) >= 0.95 * py
    assert np.min(reference[:, 1]) <= -0.95 * py
    assert np.max(reference[:, 2]) >= 0.95 * pz
    assert np.min(reference[:, 2]) <= -0.95 * pz
    # Stress intent: x/y/z peaks sit outside the terminal chance mean box.
    assert np.max(np.abs(reference[:, 0])) > mean_x
    assert np.max(np.abs(reference[:, 1])) > mean_y
    assert np.max(np.abs(reference[:, 2])) > mean_z - 0.05
    # On-face slides use the free axis (not pure axial dwell).
    on_pos_x = reference[:, 0] >= 0.9 * px
    assert np.max(np.abs(reference[on_pos_x, 1])) >= 0.5 * bh.FACE_PETAL_SLIDE
    # Smooth enough for the slow Pelican closed loop at the default density.
    # Short unit-test horizons are coarser; scale the step budget with 1/steps.
    max_step = float(np.max(np.linalg.norm(np.diff(reference, axis=0), axis=1)))
    assert max_step < max(0.08, 500.0 / steps)