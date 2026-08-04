import importlib.util
from pathlib import Path

import numpy as np
import scipy.io


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT / "experiments" / "copyBG_pelican_physical_wind_smpc"
    / "run_physical_wind_smpc.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location("copybg", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_copyau_time_contract_is_18_decisions_lifted_by_d10():
    m = load_module()
    assert m.SAMPLE_TIME_SECONDS == 0.01
    assert m.CONTROL_INTERVAL_STEPS == 10
    assert m.HORIZON_STEPS == 18
    assert (
        m.HORIZON_STEPS * m.CONTROL_INTERVAL_STEPS * m.SAMPLE_TIME_SECONDS
        == 1.8
    )


def test_mstep_lift_and_noise_match_single_step_power():
    m = load_module()
    rng = np.random.default_rng(0)
    A = rng.normal(size=(4, 4)) * 0.1 + np.eye(4) * 0.9
    B = rng.normal(size=(4, 2))
    Sigma = rng.normal(size=(4, 4))
    Sigma = Sigma @ Sigma.T
    d = 10
    Ad, Bd = m.mstep_lift(A, B, d)
    assert np.allclose(Ad, np.linalg.matrix_power(A, d))
    Ai = np.eye(4)
    Bd_ref = np.zeros_like(B)
    Sigma_d_ref = np.zeros_like(Sigma)
    for _ in range(d):
        Bd_ref = Bd_ref + Ai @ B
        Sigma_d_ref = Sigma_d_ref + Ai @ Sigma @ Ai.T
        Ai = Ai @ A
    assert np.allclose(Bd, Bd_ref)
    assert np.allclose(m.lift_noise_covariance(A, Sigma, d), Sigma_d_ref)


def test_future_reference_horizon_uses_control_interval_spacing():
    m = load_module()
    reference = np.arange(100).reshape(-1, 1) * np.ones((1, 3))
    refs = m.future_reference_horizon(reference, step=0, horizon=4, control_interval=10)
    assert refs.shape == (4, 3)
    assert np.allclose(refs[:, 0], [10, 20, 30, 40])


def test_physical_wind_velocity_is_reproducible_and_has_mps_scale():
    m = load_module()
    n = 200_000
    a = m.physical_wind_velocity(n, sigma_mps=0.5, seed=7007)
    b = m.physical_wind_velocity(n, sigma_mps=0.5, seed=7007)
    c = m.physical_wind_velocity(n, sigma_mps=0.5, seed=7008)
    assert a.shape == (n, 3)
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)
    # Colored/modulated wind → fewer independent samples; looser mean gate.
    assert np.all(np.abs(a.mean(axis=0)) < 0.05)
    # Long-run RMS (pooled) renormalized to sigma; per-axis std ≈ sigma.
    assert np.allclose(a.std(axis=0), 0.5, rtol=0.05)
    # Gusty: rolling intensity must breathe (not flat white).
    win = 800  # 8 s @ 100 Hz
    roll = np.mean(
        np.column_stack([m._rolling_rms(a[:, i], win) for i in range(3)]),
        axis=1,
    )
    # Drop warm-up; coefficient of variation of rolling std should be material.
    tail = roll[win:]
    cv = float(np.std(tail) / max(np.mean(tail), 1e-15))
    assert cv > 0.12
    # Default scale is larger than the old 0.02 white-noise demo.
    assert m.SIGMA_WIND_MPS >= 0.04
    assert m.WIND_MODULATION_DEPTH >= 0.5


def test_wind_intensity_envelope_breathes():
    m = load_module()
    env = m.wind_intensity_envelope(20_000, dt_seconds=0.01)
    assert env.shape == (20_000,)
    assert float(np.min(env)) >= 0.08
    assert float(np.max(env)) > float(np.min(env)) * 1.5
    assert float(np.std(env)) > 0.15


def test_wind_mapping_integrates_velocity_then_standardizes_position():
    m = load_module()
    task_map = np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
    ])
    y_scale_m = np.array([0.9, 1.2, 0.75])
    velocity_mps = np.array([[0.5, -0.25, 1.0], [-0.1, 0.4, 0.0]])
    latent, G_w = m.wind_velocity_to_latent_increment(
        velocity_mps, task_map, y_scale_m, dt_seconds=0.01,
    )
    expected = velocity_mps * 0.01 / y_scale_m
    reconstructed = (task_map @ latent.T).T
    assert latent.shape == (2, 5)
    assert G_w.shape == (5, 3)
    assert np.allclose(reconstructed, expected, atol=1e-14)
    assert np.allclose(latent[:, 3:], 0.0, atol=1e-14)


def test_offline_noise_objects_are_recomputed_from_training_flights():
    m = load_module()
    model, E, _, scales, noise = m.load_identification_and_noise_objects()
    assert noise["training_segments"] == 36
    assert noise["transition_count"] > 1_000_000
    assert noise["Sigma_eps"].shape == (5, 5)
    assert noise["Sigma_obs_proxy"].shape == (10, 10)
    assert np.allclose(noise["Sigma_eps"], model["Sigma_eps"][:5, :5], rtol=1e-9)
    assert np.allclose(noise["Sigma_obs_proxy"], noise["Sigma_obs_proxy"].T)
    assert np.linalg.eigvalsh(noise["Sigma_eps"]).min() >= -1e-12
    assert np.linalg.eigvalsh(noise["Sigma_obs_proxy"]).min() >= -1e-10
    assert noise["uses_true_Sigma_n"] == 0
    assert noise["obs_tracked_leakage_norm"] < 1e-9
    assert scales["position_unit"] == "m"
    assert np.allclose(E.T @ model["P"][:, :5], np.c_[np.eye(3), np.zeros((3, 2))], atol=1e-9)


def test_process_covariances_keep_identified_and_wind_sources_separate():
    m = load_module()
    model, E, _, scales, noise = m.load_identification_and_noise_objects()
    cov = m.build_process_covariances(
        model, E, scales["y_scale"][:3], sigma_wind_mps=0.5,
        dt_seconds=0.01,
    )
    assert cov["Sigma_eps_aug"].shape == model["A"].shape
    assert cov["Sigma_wind_aug"].shape == model["A"].shape
    assert np.allclose(
        cov["Sigma_total_aug"],
        cov["Sigma_eps_aug"] + cov["Sigma_wind_aug"],
    )
    assert np.allclose(cov["Sigma_eps_aug"][:5, :5], noise["Sigma_eps"])
    assert np.allclose(cov["Sigma_eps_aug"][5:, :], 0.0)
    assert np.allclose(cov["Sigma_wind_aug"][5:, :], 0.0)
    assert np.all(np.diag(cov["Sigma_wind_aug"])[:3] > 0.0)
    assert np.allclose(np.diag(cov["Sigma_wind_aug"])[3:], 0.0)
    assert "Sigma_obs_proxy" not in cov


def test_near_boundary_reference_is_continuous_closed_and_nonplanar():
    m = load_module()
    pressures = np.array([3.2, 3.2, 3.2])
    reference = m.continuous_near_boundary_reference(
        steps=1200, pressures=pressures,
    )

    assert reference.shape == (1200, 3)
    assert np.allclose(reference[0], [0.0, 0.0, 0.0])
    assert np.allclose(reference[-1], reference[0])
    assert np.isclose(reference[:, 0].max(), pressures[0], atol=2e-3)
    assert np.isclose(reference[:, 1].max(), pressures[1], atol=2e-3)
    assert np.isclose(reference[:, 2].max(), pressures[2], atol=2e-3)
    assert np.max(np.linalg.norm(np.diff(reference, axis=0), axis=1)) < 0.05
    assert np.linalg.matrix_rank(reference - reference.mean(axis=0)) == 3


def test_square_helix_reference_lies_on_l_infinity_boundary():
    m = load_module()
    pressures = np.array([3.2, 3.2, 3.2])
    reference = m.square_helix_reference(steps=12000, pressures=pressures)

    assert reference.shape == (12000, 3)
    rx, ry, z_hi = pressures
    # Every point sits on the square |x|<=rx, |y|<=ry edge (L-infinity sphere).
    radial = np.maximum(
        np.abs(reference[:, 0]) / rx, np.abs(reference[:, 1]) / ry,
    )
    assert np.allclose(radial, 1.0, atol=1e-9)
    assert np.all(np.abs(reference[:, 2]) <= z_hi + 1e-12)
    # Corners (|x|=rx and |y|=ry simultaneously) are actually visited;
    # sampled theta never hits pi/4 exactly, so use a tolerance well below
    # the 0.05-level step size.
    corner_mask = (
        np.isclose(np.abs(reference[:, 0]), rx, atol=1e-3)
        & np.isclose(np.abs(reference[:, 1]), ry, atol=1e-3)
    )
    assert corner_mask.sum() >= 4
    # Both dwells sit at the bottom/top z extremes.
    assert np.allclose(reference[:12000 // 5, 2], -z_hi)
    assert np.allclose(reference[-12000 // 5:, 2], z_hi)


def test_circle_helix_reference_is_smooth_cylinder_and_z_monotonic():
    m = load_module()
    pressures = np.array([2.0, 2.0, 3.0])
    reference = m.circle_helix_reference(steps=12000, pressures=pressures)

    assert reference.shape == (12000, 3)
    rx, ry, z_hi = pressures
    # Horizontal section is the ellipse x^2/rx^2 + y^2/ry^2 = 1 (circle if rx=ry).
    radial = (reference[:, 0] / rx) ** 2 + (reference[:, 1] / ry) ** 2
    assert np.allclose(radial, 1.0, atol=1e-9)
    assert np.all(np.abs(reference[:, 2]) <= z_hi + 1e-12)
    # No square corners: simultaneous axis extremes are not forced.
    corner_mask = (
        np.isclose(np.abs(reference[:, 0]), rx, atol=1e-3)
        & np.isclose(np.abs(reference[:, 1]), ry, atol=1e-3)
    )
    assert corner_mask.sum() == 0
    assert np.allclose(reference[:12000 // 5, 2], -z_hi)
    assert np.allclose(reference[-12000 // 5:, 2], z_hi)
    step_norms = np.linalg.norm(np.diff(reference, axis=0), axis=1)
    assert np.max(step_norms) < 0.05
    assert np.all(np.diff(reference[:, 2]) >= -1e-12)
    assert np.isclose(reference[0, 2], -z_hi)
    assert np.isclose(reference[-1, 2], z_hi)
    # One rising turn ends on +x.
    assert np.isclose(reference[-1, 0], rx, atol=1e-9)
    assert np.isclose(reference[-1, 1], 0.0, atol=1e-9)
    # Bottom dwell holds a fixed azimuth (no horizontal spinning while parked).
    assert np.allclose(reference[:12000 // 5, 0], rx)
    assert np.allclose(reference[:12000 // 5, 1], 0.0)


def test_default_reference_mode_is_growing_square_envelope():
    m = load_module()
    # Full narrative default: fly large arcs on the upper/lower z hard faces (v4.0).
    assert m.DEFAULT_REFERENCE_MODE == "growing-smooth-square-helix"
    assert "boundary-tour" in m.REFERENCE_MODES
    assert "smooth-square-helix" in m.REFERENCE_MODES
    assert "growing-smooth-square-helix" in m.REFERENCE_MODES
    assert "circle-helix" in m.REFERENCE_MODES
    assert "square-spiral" in m.REFERENCE_MODES
    assert "tour" in m.REFERENCE_MODES
    assert m.SMOOTH_SQUARE_HELIX_PRESSURES == (1.80, 1.80, 2.80)
    assert m.GROWING_SMOOTH_SQUARE_START == (0.90, 0.90, 1.20)
    assert m.GROWING_SMOOTH_SQUARE_END == (1.20, 1.20, 3.57)
    assert m.GROWING_SMOOTH_SQUARE_PEAK == m.GROWING_SMOOTH_SQUARE_END
    assert m.GROWING_DC_LABEL_R == 3.57
    assert m.GROWING_PATH_FAMILY == "z-face-large-arc-cruise"
    assert m.GROWING_NEAR_HARD_BAND == 0.40
    assert m.GROWING_NEAR_HARD_TARGET_FRACTION == 0.50
    assert m.GROWING_ARC_CENTER_OFFSET >= 8.0
    assert 0.25 <= m.GROWING_AXIS_TIP_SLOWDOWN <= 0.6
    assert m.GROWING_DC_U_INF_CAP == 6.0
    assert m.GROWING_SMOOTH_SQUARE_DEFAULT_STEPS == 48000
    assert len(m.GROWING_PHASE_FRACTIONS) == 5
    assert sum(m.GROWING_PHASE_FRACTIONS[1::2]) >= 0.60
    assert len(m.GROWING_Z_FACE_TURNS) == 5
    assert m.GROWING_Z_FACE_TURNS[1] >= 2.0
    assert m.GROWING_Z_FACE_TURNS[3] >= 2.0
    assert m.SMOOTH_SQUARE_CORNER_POWER == 4.0
    assert m.SMOOTH_SQUARE_DEFAULT_STEPS == 36000
    assert m.MAX_STEPS == 90000
    assert m.SMOOTH_SQUARE_RISING_TURNS == 2.0
    assert m.SMOOTH_SQUARE_DWELL_FRACTION == 0.10
    assert m.SMOOTH_SQUARE_DC_U_INF_CAP == 5.5
    assert m.BOUNDARY_TOUR_TARGETS == (3.60, 3.55, 3.75)
    assert m.CIRCLE_HELIX_PRESSURES == (0.55, 0.55, 0.80)
    assert m.NEAR_BOUND_CIRCLE_HELIX_PRESSURES == (2.0, 2.0, 3.0)
    assert m.LARGE_CIRCLE_HELIX_PRESSURES == (2.0, 2.0, 3.0)
    assert m.CIRCLE_HELIX_RISING_TURNS == 1.0
    assert m.SQUARE_SPIRAL_PRESSURES == (3.2, 3.2, 3.2)


def test_smooth_square_helix_is_rounded_box_and_continuous():
    m = load_module()
    pressures = np.array(m.SMOOTH_SQUARE_HELIX_PRESSURES)
    # Geometry unit-test uses single-turn + historical 0.2 dwell layout.
    reference = m.smooth_square_helix_reference(
        steps=12000, pressures=pressures, corner_power=4.0,
        rising_turns=1.0, dwell_fraction=0.2,
    )
    assert reference.shape == (12000, 3)
    # Peaks reach the box half-widths on axes.
    assert np.isclose(reference[:, 0].max(), pressures[0], atol=1e-9)
    assert np.isclose(reference[:, 1].max(), pressures[1], atol=1e-9)
    assert np.isclose(reference[:, 2].max(), pressures[2], atol=1e-9)
    assert np.isclose(reference[:, 2].min(), -pressures[2], atol=1e-9)
    # Superellipse stays inside the axis-aligned box.
    assert np.all(np.abs(reference[:, 0]) <= pressures[0] + 1e-12)
    assert np.all(np.abs(reference[:, 1]) <= pressures[1] + 1e-12)
    # Corners are filleted: n=4 L1 at 45° is smaller than sharp square (2.0).
    sharp_corner_l1 = 2.0
    dwell = int(round(12000 * 0.2))
    rise = reference[dwell: 12000 - dwell]
    l1 = np.abs(rise[:, 0]) / pressures[0] + np.abs(rise[:, 1]) / pressures[1]
    assert np.max(l1) < sharp_corner_l1 - 0.05
    step_norms = np.linalg.norm(np.diff(reference, axis=0), axis=1)
    # Arc-length sampling keeps horizontal speed nearly constant (~perimeter/rise).
    assert np.max(step_norms) < 0.01
    rise_steps = step_norms[dwell: 12000 - dwell]
    assert np.max(rise_steps) < 0.006
    assert np.all(np.diff(reference[:, 2]) >= -1e-12)
    # Smoother than sharp θ-parameterized square at same pressures.
    sharp = m.square_helix_reference(steps=12000, pressures=pressures)
    sharp_peak = np.max(np.linalg.norm(np.diff(sharp, axis=0), axis=1))
    assert np.max(step_norms) < sharp_peak
    # Multi-turn packing on the DC-feasible orbit stays slow enough.
    dense = m.smooth_square_helix_reference(steps=36000, pressures=pressures)
    dense_peak = np.max(np.linalg.norm(np.diff(dense, axis=0), axis=1))
    assert dense_peak < 0.01
    # Default rising turns: multi-rev but not denser-than-bandwidth packing.
    assert m.SMOOTH_SQUARE_RISING_TURNS >= 2.0


def test_smooth_square_dc_feasibility_gate():
    m = load_module()
    model, E, hard_bound, scales, noise = m.load_identification_and_noise_objects()
    ok = m.smooth_square_dc_feasibility(
        model, E, pressures=m.SMOOTH_SQUARE_HELIX_PRESSURES,
    )
    assert ok["feasible"] is True
    assert ok["worst_u_inf"] <= m.SMOOTH_SQUARE_DC_U_INF_CAP
    bad = m.smooth_square_dc_feasibility(
        model, E, pressures=(3.5, 3.5, 3.5),
    )
    assert bad["feasible"] is False
    assert bad["worst_u_inf"] > 6.0
    _ = hard_bound, scales, noise


def test_steady_state_input_uses_motor_nullspace_to_minimize_peak_command():
    m = load_module()
    model, E, hard_bound, _scales, _noise = m.load_identification_and_noise_objects()

    # Four motor commands control three tracked position coordinates.  The
    # minimum-2-norm pseudoinverse solution exceeds ±6 on the y hard face, but
    # its one-dimensional nullspace contains an exact equilibrium below ±6.
    target = np.array([0.0, hard_bound, 0.0])
    info = m.steady_state_task_input(model, E, target)

    assert np.allclose(info["y_hat"], target, atol=1e-9)
    assert info["u_inf"] < 6.0


def test_z_face_large_arc_cruise_is_continuous_and_dc_feasible():
    m = load_module()
    start = np.array(m.GROWING_SMOOTH_SQUARE_START)
    peak = np.array(m.GROWING_SMOOTH_SQUARE_PEAK)
    steps = 48000
    reference = m.growing_smooth_square_helix_reference(
        steps=steps,
        pressures_start=start,
        pressures_end=peak,
    )
    assert reference.shape == (steps, 3)
    phase_ids, counts = m.growing_envelope_phase_index(steps)
    assert int(np.sum(counts)) == steps
    assert counts.shape == (5,)
    top = reference[phase_ids == 1]
    bottom = reference[phase_ids == 3]
    assert np.allclose(top[:, 2], peak[2], atol=1e-12)
    assert np.allclose(bottom[:, 2], -peak[2], atol=1e-12)
    # It is flying on both hard faces, not hovering at one point.
    assert float(np.ptp(top[:, 0])) > 1.5
    assert float(np.ptp(top[:, 1])) > 1.5
    assert float(np.ptp(bottom[:, 0])) > 1.5
    assert float(np.ptp(bottom[:, 1])) > 1.5
    # Horizontal large arcs still cut the forbidden square corners.
    joint = np.minimum(np.abs(reference[:, 0]) / peak[0], np.abs(reference[:, 1]) / peak[1])
    assert float(np.max(joint)) < 0.70
    # Continuous speed with no teleport between the five phases.
    sn = np.linalg.norm(np.diff(reference, axis=0), axis=1)
    assert float(np.max(sn)) < 0.01
    phases = m.growing_smooth_square_phase_masks(
        reference, hard_bound=3.91, phase_ids=phase_ids, peak_pressures=peak,
    )
    assert phases["fraction_near_hard"] >= 0.60
    assert float(phases["r_linf_max"]) <= peak[2] + 5e-2
    assert float(np.max(peak)) < float(3.91)
    # Both z faces plus the horizontal orbit are sustainable under ±6 input hard limits.
    model, E, hard_bound, _scales, _noise = m.load_identification_and_noise_objects()
    dc = m.z_face_large_arc_cruise_dc_feasibility(
        model, E, xy_radius=peak[0], z_face=peak[2], u_inf_cap=m.GROWING_DC_U_INF_CAP,
    )
    assert dc["feasible"] is True
    assert float(dc["worst_u_inf"]) <= m.GROWING_DC_U_INF_CAP
    assert float(dc["max_joint_xy_unit"]) < 0.55
    assert float(hard_bound) > peak[2]


def test_large_arc_diamond_unit_cuts_square_corners():
    m = load_module()
    table = m._unit_large_arc_diamond_table(
        center_offset=m.GROWING_ARC_CENTER_OFFSET, n_dense=4096,
    )
    assert table["family"] == "large-arc-diamond"
    assert table["max_joint_xy"] < 0.55
    s = np.linspace(0.0, table["perimeter"], 720, endpoint=False)
    x, y = m._unit_xy_at_arc(s, table)
    assert np.isclose(np.max(np.abs(x)), 1.0, atol=1e-6)
    assert np.isclose(np.max(np.abs(y)), 1.0, atol=1e-6)
    # No sample near the forbidden square corner (1,1).
    assert not np.any((np.abs(x) > 0.85) & (np.abs(y) > 0.85))


def test_default_growing_reference_has_majority_near_hard_flight():
    m = load_module()
    hard = 3.9124117280376636
    reference = m.growing_smooth_square_helix_reference(steps=48_000)
    occupancy = m.near_hard_occupancy(
        reference, hard_bound=hard, band=m.GROWING_NEAR_HARD_BAND,
    )
    assert occupancy["threshold"] == hard - m.GROWING_NEAR_HARD_BAND
    assert occupancy["fraction"] >= 0.60
    assert occupancy["violation_steps"] == 0


def test_noise_actual_vs_estimated_summary_and_figure(tmp_path):
    m = load_module()
    rng = np.random.default_rng(0)
    ell = 5
    Sigma_eps = np.diag(np.array([0.04, 0.03, 0.02, 0.01, 0.005]) ** 2)
    Sigma_wind = np.diag(np.array([0.001, 0.001, 0.0005, 0.0, 0.0]) ** 2)
    n = 8000
    eps = rng.multivariate_normal(np.zeros(ell), Sigma_eps, size=n)
    wind_z = rng.multivariate_normal(np.zeros(ell), Sigma_wind, size=n)
    sigma_w = 0.02
    wind_v = sigma_w * rng.standard_normal((n, 3))
    Sigma_eps_aug = np.zeros((2 * ell, 2 * ell))
    Sigma_eps_aug[:ell, :ell] = Sigma_eps
    Sigma_wind_aug = np.zeros_like(Sigma_eps_aug)
    Sigma_wind_aug[:ell, :ell] = Sigma_wind
    Sigma_total_aug = Sigma_eps_aug + Sigma_wind_aug
    summary = m.summarize_noise_actual_vs_estimated(
        eps, wind_v, wind_z, Sigma_eps, Sigma_wind_aug, Sigma_total_aug, sigma_w,
    )
    assert summary["n_steps"] == n
    assert np.allclose(
        summary["wind_velocity_std_ratio_actual_over_est"], 1.0, atol=0.05,
    )
    assert np.allclose(
        summary["eps_std_ratio_actual_over_est"], 1.0, atol=0.08,
    )
    assert summary["fro_rel_emp_eps_vs_Sigma_eps"] < 0.15
    path = tmp_path / "noise.png"
    meta = m.plot_noise_actual_vs_estimated(
        path, eps, wind_v, wind_z, Sigma_eps, Sigma_wind_aug, Sigma_total_aug,
        sigma_w,
    )
    assert path.is_file() and path.stat().st_size > 5000
    assert meta["figure"] == path.name


def test_forward_view_gif_writes_file(tmp_path):
    m = load_module()
    steps = 400
    reference = m.smooth_square_helix_reference(
        steps=steps, pressures=(1.2, 1.2, 1.0),
        rising_turns=2.0, dwell_fraction=0.1,
    )
    # Slight lag so the craft is not identical to the reference ghost.
    trajectory = reference.copy()
    trajectory[1:] = 0.92 * reference[:-1] + 0.08 * reference[1:]
    path = tmp_path / "forward.gif"
    meta = m.render_forward_view_gif(
        path, trajectory=trajectory, reference=reference,
        hard_bound=3.9, stage_mean_bound=np.tile([2.5, 2.5, 2.5], (4, 1)),
        max_frames=24, fps=8,
    )
    assert path.is_file()
    assert path.stat().st_size > 1000
    assert meta["frames"] == 24
    assert meta["completed_steps"] == steps


def test_continuous_boundary_tour_matches_copybf_geometry():
    m = load_module()
    targets = np.array(m.BOUNDARY_TOUR_TARGETS)
    reference = m.continuous_boundary_tour(steps=12000, targets=targets)
    assert reference.shape == (12000, 3)
    assert np.allclose(reference[0], [targets[0], 0.0, 0.0], atol=1e-12)
    assert np.allclose(reference[-1], [targets[0], 0.0, 0.0], atol=1e-9)
    # Midpoints of the three legs approach the axis targets.
    assert np.allclose(reference[4000], [0.0, targets[1], 0.0], atol=0.05)
    assert np.allclose(reference[8000], [0.0, 0.0, targets[2]], atol=0.05)
    step_norms = np.linalg.norm(np.diff(reference, axis=0), axis=1)
    assert np.max(step_norms) < 0.01
    assert np.all(np.max(np.abs(reference), axis=0) <= targets + 1e-12)


def test_square_helix_reference_is_continuous_and_z_monotonic():
    m = load_module()
    pressures = np.array([2.0, 2.5, 3.0])
    reference = m.square_helix_reference(steps=5000, pressures=pressures)

    step_norms = np.linalg.norm(np.diff(reference, axis=0), axis=1)
    assert np.max(step_norms) < 0.1
    # z rises monotonically across the whole path (bottom dwell, rise, top dwell).
    assert np.all(np.diff(reference[:, 2]) >= -1e-12)
    assert np.isclose(reference[0, 2], -3.0)
    assert np.isclose(reference[-1, 2], 3.0)
    assert np.isclose(reference[-1, 0], 2.0)  # theta=10*pi ends back on +x edge


def test_equilibrium_initial_condition_matches_task_and_dynamics():
    m = load_module()
    model, E, _, _, _ = m.load_identification_and_noise_objects()
    target = np.array([2.8, 0.0, 0.0])
    state, control = m.equilibrium_initial_condition(model, E, target)

    task = (E.T @ (model["y_mean"] + model["P"] @ state)).ravel()
    successor = model["A"] @ state + model["B"] @ (
        control.reshape(-1, 1) - model["u_mean"]
    )
    assert np.allclose(task, target, atol=1e-7)
    assert np.allclose(successor, state, atol=1e-7)
    assert np.all(np.abs(control) <= 6.0 + 1e-9)


def test_cached_controller_reuses_one_dpp_qp_and_bounds_full_plan():
    m = load_module()
    model, E, hard_bound, scales, _ = m.load_identification_and_noise_objects()
    cov = m.build_process_covariances(model, E, scales["y_scale"][:3])
    target = np.array([2.8, 0.0, 0.0])
    state, _ = m.equilibrium_initial_condition(model, E, target)
    controller = m.CachedBoundaryController(
        model, E, hard_bound, cov["Sigma_total_aug"],
        horizon=4, alpha_joint=0.05,
    )
    problem_identity = id(controller.problem)
    control, info = controller.step(state, np.tile(target, (4, 1)))

    assert id(controller.problem) == problem_identity
    assert controller.problem.is_dpp()
    assert controller.U.shape == (4, model["B"].shape[1])
    assert controller.state_parameter.shape == (model["A"].shape[0],)
    assert controller.reference_parameter.shape == (4, 3)
    assert control is not None
    assert info["planned_controls"].shape == (4, model["B"].shape[1])
    assert np.max(np.abs(info["planned_controls"])) <= 6.0 + 1e-5
    assert controller.stage_tightening.shape == (4, 3)
    assert np.all(controller.stage_tightening > 0.0)


def test_short_pair_uses_identical_noise_and_stays_on_qp_path():
    m = load_module()
    model, E, hard_bound, scales, _ = m.load_identification_and_noise_objects()
    cov = m.build_process_covariances(model, E, scales["y_scale"][:3])
    # 30 plant steps => 3 QP solves under d=10 hold.
    steps = 30
    reference = np.tile(np.array([0.4, 0.0, 0.0]), (steps, 1))
    innovation = np.zeros((steps, 5))
    wind_latent = np.zeros((steps, 5))
    smpc, mpc = m.run_controller_pair(
        model, E, hard_bound, reference, innovation, wind_latent,
        cov["Sigma_total_aug"], horizon=3,
    )

    assert smpc["S"].shape == (steps, 3)
    assert mpc["S"].shape == (steps, 3)
    assert smpc["qp_count"] == mpc["qp_count"] == 3
    assert smpc["fallback_count"] == mpc["fallback_count"] == 0
    assert smpc["disturbance_sha256"] == mpc["disturbance_sha256"]
    assert smpc["control_interval_steps"] == 10
    assert np.allclose(smpc["S"][0], reference[0], atol=1e-7)
    assert np.allclose(mpc["S"][0], reference[0], atol=1e-7)


def test_controller_state_is_reconstructed_from_current_and_previous_outputs():
    m = load_module()
    model, _, _, _, _ = m.load_identification_and_noise_objects()
    ell = model["A"].shape[0] // 2
    current = np.linspace(-0.3, 0.4, ell).reshape(-1, 1)
    previous = np.linspace(0.2, -0.1, ell).reshape(-1, 1)
    output = model["y_mean"] + model["P"][:, :ell] @ current

    observed = m.observed_augmented_state(model, output, previous)

    assert observed.shape == model["A"].shape[:1] + (1,)
    assert np.allclose(observed[:ell], current, atol=1e-12)
    assert np.allclose(observed[ell:], previous, atol=1e-12)


def test_qp_failure_stops_simulation_without_fallback():
    m = load_module()
    model, E, _, scales, _ = m.load_identification_and_noise_objects()
    cov = m.build_process_covariances(
        model, E, scales["y_scale"][:3], sigma_wind_mps=0.0,
    )
    steps = 5
    # Target far outside a near-zero hard box => first QP is infeasible.
    reference = np.tile(np.array([2.8, 0.0, 0.0]), (steps, 1))
    zeros = np.zeros((steps, 5))

    result = m.simulate_controller(
        model, E, 1e-6, reference, zeros, zeros,
        cov["Sigma_total_aug"], alpha_joint=None, horizon=3,
    )

    assert result["completed_steps"] == 0
    assert result["qp_failure_count"] == 1
    assert result["qp_failure_step"] == 0
    assert result["fallback_count"] == 0
    assert result["S"].shape == (0, 3)
    assert result["U"].shape == (0, 4)


def test_mat_payload_preserves_three_noise_channels_and_meter_units(tmp_path):
    m = load_module()
    model, E, hard_bound, scales, noise = m.load_identification_and_noise_objects()
    cov = m.build_process_covariances(model, E, scales["y_scale"][:3])
    n = 5
    reference = np.zeros((n, 3))
    result = {
        "S": np.zeros((n, 3)), "U": np.zeros((n, 4)),
        "violation_steps": 0, "violation_rate": 0.0,
        "fallback_count": 0, "active_qp_steps": 2,
        "positive_dual_qp_steps": 2,
        "minimum_hard_margin": 1.0,
        "minimum_hard_margin_per_axis": np.ones(3),
        "rmse": np.zeros(3), "input_saturation_steps": 0,
    }
    innovation = np.zeros((n, 5))
    wind_velocity = np.full((n, 3), 0.5)
    wind_latent, _ = m.wind_velocity_to_latent_increment(
        wind_velocity, E.T @ model["P"][:, :5], scales["y_scale"][:3],
    )
    tightening = np.ones((m.HORIZON_STEPS, 3)) * 0.2
    payload = m.build_mat_payload(
        reference, result, result, innovation, wind_velocity, wind_latent,
        cov, noise, scales, hard_bound, tightening,
    )
    required = {
        "reference_standardized", "reference_position_m",
        "smpc_trajectory_standardized", "smpc_position_m",
        "deterministic_mpc_trajectory_standardized",
        "deterministic_mpc_position_m", "smpc_control",
        "deterministic_mpc_control", "residual_innovation_latent",
        "wind_velocity_mps", "wind_position_delta_m",
        "wind_position_delta_standardized", "wind_latent_increment",
        "Sigma_eps", "Sigma_obs_proxy", "Sigma_wind_aug",
        "Sigma_total_aug", "G_w", "dt_seconds", "sample_rate_hz",
        "horizon_steps", "horizon_seconds", "sigma_wind_mps",
        "smpc_hard_violation_rate",
        "deterministic_mpc_hard_violation_rate",
        "smpc_fallback_count", "deterministic_mpc_fallback_count",
        "stage_tightening_standardized", "stage_mean_bound_standardized",
        "uses_true_Sigma_n", "unit_note",
    }
    assert required <= payload.keys()
    assert np.allclose(payload["wind_position_delta_m"], 0.005)
    assert float(payload["dt_seconds"][0]) == 0.01
    assert "m/s" in payload["unit_note"]

    path = tmp_path / "copybg_contract.mat"
    scipy.io.savemat(path, payload, do_compression=True)
    saved = scipy.io.loadmat(path)
    assert required <= saved.keys()
    assert saved["wind_velocity_mps"].shape == (n, 3)
    assert saved["Sigma_eps"].shape == (5, 5)
    assert saved["Sigma_obs_proxy"].shape == (10, 10)


def test_json_ready_preserves_boolean_types():
    m = load_module()
    assert m._json_ready(True) is True
    assert m._json_ready(False) is False
    assert m._json_ready(np.bool_(True)) is True
    assert m._json_ready({"met": np.bool_(False)}) == {"met": False}


def test_short_smoke_writes_all_artifacts_without_fallback(tmp_path):
    m = load_module()
    summary = m.run_experiment(
        steps=40, sigma_wind_mps=0.01, hard_bound=3.4,
        pressures=(0.2, 0.2, 0.2), reference_mode="circle-helix",
        output_dir=tmp_path, make_forward_gif=False,
    )

    assert summary["requested_steps"] == 40
    assert summary["shared_noise_exact"] == 1
    assert summary["smpc"]["fallback_count"] == 0
    assert summary["deterministic_mpc"]["fallback_count"] == 0
    assert summary["smpc"]["maximum_qp_constraint_residual"] <= 1e-6
    assert summary["deterministic_mpc"]["maximum_qp_constraint_residual"] <= 1e-6
    for artifact in summary["artifacts"].values():
        path = tmp_path / artifact
        assert path.is_file()
        assert path.stat().st_size > 0

    arrays = np.load(tmp_path / summary["artifacts"]["npz"])
    assert arrays["reference_standardized"].shape == (40, 3)
    assert arrays["wind_velocity_mps"].shape == (40, 3)
    assert np.array_equal(arrays["shared_noise_exact"], [1])
