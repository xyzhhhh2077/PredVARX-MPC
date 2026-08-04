"""Pelican frozen-model SMPC with separately modeled physical Gaussian wind.

The identified model and both data-derived covariance objects come only from
Waterloo Pelican training flights 1-36. Closed-loop wind is an independent
velocity disturbance in m/s, integrated for one 100 Hz sample and mapped into
standardized position coordinates before entering the latent state.
"""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import cvxpy as cp
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
from scipy.linalg import null_space
from scipy.io import loadmat
from scipy.optimize import linprog
from scipy.stats import norm


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BE_DIR = HERE.parent / "copyBE_pelican_3d_z_aware_smpc"
BE_SCRIPT = BE_DIR / "run_3d_z_aware_trajectory.py"
DATA_MAT = (
    HERE.parent / "copyAY_pelican_soft_preference" / "data"
    / "copyAY_pelican_dataset.mat"
)
MODEL_MAT = BE_DIR / "results" / "copyBE_pelican_3d_z_aware_data.mat"

_spec = importlib.util.spec_from_file_location("copybg_copybe", BE_SCRIPT)
BE = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(BE)

SAMPLE_TIME_SECONDS = 0.01
SAMPLE_RATE_HZ = 100.0
# Multirate lift (copyBA/BE/BF): N=18 at every 100 Hz sample only sees 0.18 s
# and the QP backs off on slow Pelican modes. Decision every d=10 samples
# (0.1 s) makes N=18 cover 1.8 s. Q/Ru/N/dt/input bounds stay frozen.
CONTROL_INTERVAL_STEPS = 10
HORIZON_STEPS = 18
Q_WEIGHT = 80.0
RU = 0.18
ALPHA_JOINT = 0.05
# Physical wind velocity scale (m/s). v3.8: 2× old 0.02 + gusty intensity envelope.
# 0.06+ finishes 48k but RMSE/sat degrade hard; 0.08 can QP-fail mid-run.
# Chance QP still uses Σε only; wind hits the plant only.
SIGMA_WIND_MPS = 0.04
# Slow intensity envelope depth in [0, 0.95]: rolling std breathes between
# roughly (1-depth)…(1+depth) × long-run σ (after RMS renormalization).
WIND_MODULATION_DEPTH = 0.75
WIND_MODULATION_PERIOD_S = 40.0  # primary gust-intensity cycle
# Optional AR(1) color at 100 Hz. Default 0 = white × envelope only.
# Strong color (e.g. 0.97) integrates into large drift and can make the QP
# infeasible; keep ≤0.85 if enabling.
WIND_COLOR_RHO = 0.0
# Helix-mode default box (stricter than source). Boundary-tour uses the
# identification source hard bound (~3.91) unless --hard-bound overrides.
HARD_BOUND_STANDARDIZED = 3.4
REFERENCE_PRESSURES_STANDARDIZED = (1.5, 1.5, 3.2)
# copyBF continuous 3D boundary tour: large axial targets that still track
# (RMSE ~0.10) with d=10, plant-state feedback, and source hard bound.
BOUNDARY_TOUR_TARGETS = (3.60, 3.55, 3.75)
# Rounded-square rising helix (smooth corners). Pressures are a *closed-loop
# feasibility calibration*, not a theory object and not CRTE/copyAU.
#
# DC audit on the frozen Pelican model (do not retune identification):
#   simultaneous superellipse corner (n=4) at r_xy=3.5 needs u_inf≈8.79 > 6
#   → structurally unreachable; closed loop saturates ~85% and shrinks radius.
#   r_xy=1.8 keeps corner steady u_inf≈4.61 (headroom under ±6); z is cheap.
# Large single-axis near-bound demos stay on boundary-tour (3.60/3.55/3.75).
SMOOTH_SQUARE_HELIX_PRESSURES = (1.80, 1.80, 2.80)
SMOOTH_SQUARE_CORNER_POWER = 4.0  # 2=ellipse, →∞=sharp square; 4–6 = rounded box
# Steady-state |u|_inf headroom gate used only as experiment diagnostic.
SMOOTH_SQUARE_DC_U_INF_CAP = 5.5
DEFAULT_STEPS = 12000
# Dense multi-turn rounded box on the DC-feasible orbit.
SMOOTH_SQUARE_DEFAULT_STEPS = 36000
MAX_STEPS = 90000
# Rising circular helix defaults. Smaller radii stay trackable as an optional
# cylindrical mode; not the primary large-range demonstration path.
CIRCLE_HELIX_PRESSURES = (0.55, 0.55, 0.80)
NEAR_BOUND_CIRCLE_HELIX_PRESSURES = (2.0, 2.0, 3.0)
SQUARE_SPIRAL_PRESSURES = (3.2, 3.2, 3.2)
LARGE_CIRCLE_HELIX_PRESSURES = (2.0, 2.0, 3.0)
# One rising turn (was three): keeps mean |Δref| near the BE orbit (~0.001–0.002).
CIRCLE_HELIX_RISING_TURNS = 1.0
# Smooth-square default: 2 rising turns (readable orbit; still multi-rev).
SMOOTH_SQUARE_RISING_TURNS = 2.0
# Fraction of total samples spent at each of bottom/top dwells (rise gets the rest).
SMOOTH_SQUARE_DWELL_FRACTION = 0.10
# Input coordinates are standardized commanded motor speeds:
#   u_std = (Motors_CMD - u_offset) / u_scale.
# The ±6 cap is an inherited copyBA/copyBE experiment-design boundary, roughly
# a six-training-standard-deviation envelope. It is not a Pelican hardware
# rating. Dataset documentation gives a dimensionless Motors_CMD domain [0,218],
# but the 54 prepared flight records occupy a substantially narrower range.
INPUT_COMMAND_BOUND_STANDARDIZED = 6.0
# Growing envelope (v4.0): fly horizontal large arcs on z=±3.57 hard faces.
# Minimax DC audit: x near-hard exceeds |u_std|<=6 (u_inf=6.152), while y is
# input-feasible (5.151) but its horizon-18 chance mean bound is only 3.148.
# The full z-face orbit stays below u_inf=2.281 and z=3.57 stays inside its
# horizon-18 chance mean bound 3.583, so z is the sustainable default face.
GROWING_SMOOTH_SQUARE_START = (0.90, 0.90, 1.20)
GROWING_SMOOTH_SQUARE_END = (1.20, 1.20, 3.57)
GROWING_SMOOTH_SQUARE_PEAK = GROWING_SMOOTH_SQUARE_END
GROWING_SMOOTH_SQUARE_DEFAULT_STEPS = 48000
GROWING_SMOOTH_SQUARE_RISING_TURNS = 1.50
GROWING_SMOOTH_SQUARE_PEAK_TURNS = 3.0
# Fractions: approach top | top cruise | transfer | bottom cruise | return.
GROWING_PHASE_FRACTIONS = (0.15, 0.30, 0.10, 0.30, 0.15)
GROWING_Z_FACE_TURNS = (1.50, 3.0, 1.0, 3.0, 1.50)
GROWING_SMOOTH_SQUARE_DWELL_FRACTION = 0.0
# Path family for growing envelope horizontal section.
# "z-face-large-arc-cruise" = two sustained z-face moving cruises (v4.0).
# "large-arc-diamond" = legacy unit loop via 4 large circular arcs (v3.9).
# "superellipse" = legacy n-power box (use corner_power>=2).
GROWING_PATH_FAMILY = "z-face-large-arc-cruise"
# Near-hard means inside the task hard box and no farther than this band.
GROWING_NEAR_HARD_BAND = 0.40
GROWING_NEAR_HARD_TARGET_FRACTION = 0.50
# Center offset (unit path) for large arcs: larger → flatter sides, closer to
# pure diamond, weaker multi-axis mid-arc demand. ≥8 strongly cuts square corners.
GROWING_ARC_CENTER_OFFSET = 12.0
# Mild tip slowdown: strong tip-only slow made mid-arc too fast at large R.
GROWING_AXIS_TIP_SLOWDOWN = 0.35
# DC r label for phase masks / notes (peak axial half-width).
GROWING_DC_LABEL_R = 3.57
# Sustained face cruise must fit the implemented standardized-command boundary.
GROWING_DC_U_INF_CAP = INPUT_COMMAND_BOUND_STANDARDIZED
# Playback GIF (MATLAB-style fixed camera — not chase/FPV).
PLAYBACK_GIF_MAX_FRAMES = 200  # denser temporal sample → less GIF "teleport" shake
PLAYBACK_GIF_FRAME_MS = 120  # ~8.3 fps; 200 frames ≈ 24 s
# Legacy chase-cam knobs (kept only if someone calls render_forward_view_gif).
FORWARD_GIF_MAX_FRAMES = 120
FORWARD_GIF_FPS = 4

REFERENCE_MODES = (
    "boundary-tour",
    "smooth-square-helix",
    "growing-smooth-square-helix",
    "circle-helix",
    "square-spiral",
    "tour",
)
# Default: large-arc axial near-hard envelope (v3.9) — chance-QP friendly.
DEFAULT_REFERENCE_MODE = "growing-smooth-square-helix"
WIND_SEED = 7007
INNOVATION_SEED = 7
OUTPUT_STEM = "copyBG_pelican_physical_wind_smpc"


def _model_h5_array(group, key):
    return np.asarray(group[key], dtype=float).T


def _symmetrize_psd(matrix, floor=0.0):
    matrix = (np.asarray(matrix, dtype=float) + np.asarray(matrix, dtype=float).T) / 2.0
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * np.maximum(values, float(floor))) @ vectors.T


def load_identification_and_noise_objects():
    """Load frozen copyBE geometry and recompute CRTE noise objects from train data."""
    model, E, hard_bound, scales = BE.load_model()
    raw = loadmat(DATA_MAT, squeeze_me=True, struct_as_record=False)
    y_raw = np.asarray(raw["y"], dtype=float)
    u_raw = np.asarray(raw["u"], dtype=float)
    segment_id = np.asarray(raw["segment_id"], dtype=int).ravel()
    train = np.isin(segment_id, np.arange(1, 37))
    y = (y_raw[:, train] - scales["y_offset"][:, None]) / scales["y_scale"][:, None]

    with h5py.File(MODEL_MAT, "r") as f:
        m = f["model"]
        A1 = _model_h5_array(m, "A1")
        A2 = _model_h5_array(m, "A2")
        B_base = _model_h5_array(m, "B_base")
        R = _model_h5_array(m, "R")
        P_base = _model_h5_array(m, "P_base")
        y_mean = _model_h5_array(m, "y_mean")
        u_mean = _model_h5_array(m, "u_mean")
        u_offset = np.asarray(f["u_offset"], dtype=float).ravel()
        u_scale = np.asarray(f["u_scale"], dtype=float).ravel()

    u = (u_raw[:, train] - u_offset[:, None]) / u_scale[:, None]
    train_segments = segment_id[train]
    z = R.T @ (y - y_mean)
    uc = u - u_mean
    valid = np.flatnonzero(
        (train_segments[:-2] == train_segments[1:-1])
        & (train_segments[1:-1] == train_segments[2:])
    ) + 1
    innovations = (
        z[:, valid + 1] - A1 @ z[:, valid] - A2 @ z[:, valid - 1]
        - B_base @ uc[:, valid]
    )
    Sigma_eps = innovations @ innovations.T / innovations.shape[1]
    Sigma_eps = _symmetrize_psd(Sigma_eps)

    projection_residual = (np.eye(y.shape[0]) - P_base @ R.T) @ (y - y_mean)
    projection_residual -= projection_residual.mean(axis=1, keepdims=True)
    Sigma_obs_proxy = (
        projection_residual @ projection_residual.T
        / max(projection_residual.shape[1] - 1, 1)
    )
    Sigma_obs_proxy = _symmetrize_psd(Sigma_obs_proxy)

    ell = R.shape[1]
    Sigma_eps_aug = np.zeros_like(model["A"])
    Sigma_eps_aug[:ell, :ell] = Sigma_eps
    model["Sigma_eps"] = Sigma_eps_aug
    scales = dict(scales)
    scales.update({
        "u_offset": u_offset,
        "u_scale": u_scale,
        "position_unit": "m",
    })
    tracked_residual = E.T @ projection_residual
    noise = {
        "Sigma_eps": Sigma_eps,
        "Sigma_obs_proxy": Sigma_obs_proxy,
        "training_segments": int(np.unique(train_segments).size),
        "training_samples": int(y.shape[1]),
        "transition_count": int(valid.size),
        "uses_true_Sigma_n": 0,
        "obs_tracked_leakage_norm": float(
            np.linalg.norm(tracked_residual, ord="fro")
            / max(np.linalg.norm(projection_residual, ord="fro"), 1e-15)
        ),
    }
    return model, E, hard_bound, scales, noise


def physical_wind_velocity(
        steps,
        sigma_mps=SIGMA_WIND_MPS,
        seed=WIND_SEED,
        dt_seconds=SAMPLE_TIME_SECONDS,
        modulation_depth=WIND_MODULATION_DEPTH,
        modulation_period_s=WIND_MODULATION_PERIOD_S,
        color_rho=WIND_COLOR_RHO):
    """Physical wind velocity [m/s], larger and intensity-modulated.

    - Default: white driving noise × slow multiplicative intensity envelope
      → rolling std visibly breathes (panel 2), amplitude larger than 0.02
    - Optional mild AR(1) color via ``color_rho`` (keep low; strong color drifts)
    - Long-run RMS (all axes pooled) renormalized to ``sigma_mps`` so the
      offline stationary Σ_wind (built from the same σ) still matches global std

    Chance QP still ignores wind (Σε only); this sequence only enters the plant
    via ``wind_velocity_to_latent_increment``. Not a CRTE/copyAU change.
    """
    steps = int(steps)
    if steps <= 0:
        raise ValueError("steps must be positive")
    rng = np.random.default_rng(int(seed))
    sigma = float(sigma_mps)
    dt = float(dt_seconds)
    depth = float(np.clip(modulation_depth, 0.0, 0.95))
    period = max(float(modulation_period_s), 1.0)
    rho = float(np.clip(color_rho, 0.0, 0.999))

    white = rng.standard_normal((steps, 3))
    if rho > 0.0:
        innov_scale = float(np.sqrt(max(1.0 - rho * rho, 1e-12)))
        colored = np.empty((steps, 3), dtype=float)
        colored[0] = white[0]
        for t in range(1, steps):
            colored[t] = rho * colored[t - 1] + innov_scale * white[t]
        drive = colored
    else:
        drive = white

    t = np.arange(steps, dtype=float) * dt
    # Two incommensurate periods → irregular gust breathing, not a pure sine.
    env = (
        1.0
        + depth * np.sin(2.0 * np.pi * t / period)
        + 0.45 * depth * np.sin(
            2.0 * np.pi * t / (0.37 * period + 5.0) + 0.7
        )
    )
    env = np.maximum(env, 0.08)
    raw = drive * env[:, None]

    rms = float(np.sqrt(np.mean(raw * raw)))
    if rms < 1e-15 or sigma <= 0.0:
        return np.zeros((steps, 3), dtype=float)
    return (sigma / rms) * raw


def wind_intensity_envelope(
        steps,
        dt_seconds=SAMPLE_TIME_SECONDS,
        modulation_depth=WIND_MODULATION_DEPTH,
        modulation_period_s=WIND_MODULATION_PERIOD_S):
    """Prescribed relative intensity envelope (mean-normalized shape, pre-RMS)."""
    steps = int(steps)
    dt = float(dt_seconds)
    depth = float(np.clip(modulation_depth, 0.0, 0.95))
    period = max(float(modulation_period_s), 1.0)
    t = np.arange(steps, dtype=float) * dt
    env = (
        1.0
        + depth * np.sin(2.0 * np.pi * t / period)
        + 0.45 * depth * np.sin(
            2.0 * np.pi * t / (0.37 * period + 5.0) + 0.7
        )
    )
    return np.maximum(env, 0.08)


def wind_velocity_to_latent_increment(velocity_mps, task_map, y_scale_m,
                                      dt_seconds=SAMPLE_TIME_SECONDS):
    """Integrate wind velocity and map its standardized xyz displacement to z."""
    velocity_mps = np.asarray(velocity_mps, dtype=float)
    task_map = np.asarray(task_map, dtype=float)
    y_scale_m = np.asarray(y_scale_m, dtype=float).reshape(3)
    standardized_position_map = np.diag(float(dt_seconds) / y_scale_m)
    G_w = np.linalg.pinv(task_map) @ standardized_position_map
    latent_increment = (G_w @ velocity_mps.T).T
    reconstruction_error = np.max(np.abs(
        (task_map @ latent_increment.T).T
        - velocity_mps * float(dt_seconds) / y_scale_m
    ))
    if reconstruction_error > 1e-10:
        raise ValueError(f"wind-to-latent mapping error {reconstruction_error:.3e}")
    return latent_increment, G_w


def build_process_covariances(model, E, y_scale_m,
                               sigma_wind_mps=SIGMA_WIND_MPS,
                               dt_seconds=SAMPLE_TIME_SECONDS):
    ell = model["A"].shape[0] // 2
    task_map = E.T @ model["P"][:, :ell]
    _, G_w = wind_velocity_to_latent_increment(
        np.zeros((1, 3)), task_map, y_scale_m, dt_seconds,
    )
    Sigma_wind = G_w @ (float(sigma_wind_mps) ** 2 * np.eye(3)) @ G_w.T
    Sigma_eps_aug = np.array(model["Sigma_eps"], dtype=float, copy=True)
    Sigma_wind_aug = np.zeros_like(Sigma_eps_aug)
    Sigma_wind_aug[:ell, :ell] = _symmetrize_psd(Sigma_wind)
    return {
        "G_w": G_w,
        "Sigma_eps_aug": Sigma_eps_aug,
        "Sigma_wind_aug": Sigma_wind_aug,
        "Sigma_total_aug": _symmetrize_psd(Sigma_eps_aug + Sigma_wind_aug),
        "sigma_wind_mps": float(sigma_wind_mps),
        "dt_seconds": float(dt_seconds),
        "y_scale_m": np.asarray(y_scale_m, dtype=float).reshape(3),
    }


def innovation_sequence(model, steps, seed=INNOVATION_SEED):
    ell = model["A"].shape[0] // 2
    rng = np.random.default_rng(int(seed))
    return rng.multivariate_normal(
        np.zeros(ell), model["Sigma_eps"][:ell, :ell], size=int(steps),
    )


def _smoothstep(values):
    values = np.asarray(values, dtype=float)
    return values * values * (3.0 - 2.0 * values)


def continuous_near_boundary_reference(
        steps=1200, pressures=REFERENCE_PRESSURES_STANDARDIZED):
    """Create a C1 center-to-x-to-y-to-z-to-center pressure tour."""
    steps = int(steps)
    pressures = np.asarray(pressures, dtype=float).reshape(3)
    if steps < 40:
        raise ValueError("At least 40 samples are required for the continuous tour")
    if np.any(pressures <= 0.0):
        raise ValueError("All pressure magnitudes must be positive")

    px, py, pz = pressures
    anchors = np.array([
        [0.0, 0.0, 0.0],
        [px, 0.0, 0.0],
        [0.0, py, 0.0],
        [0.0, 0.0, pz],
        [0.0, 0.0, 0.0],
    ])
    transition_count = len(anchors) - 1
    intervals = np.full(
        transition_count, (steps - 1) // transition_count, dtype=int,
    )
    intervals[:(steps - 1) - int(intervals.sum())] += 1

    points = [anchors[0]]
    for segment, count in enumerate(intervals):
        tau = np.arange(1, count + 1, dtype=float) / float(count)
        blend = _smoothstep(tau)[:, None]
        values = anchors[segment] + blend * (anchors[segment + 1] - anchors[segment])
        points.extend(values)
    reference = np.asarray(points, dtype=float)
    if reference.shape != (steps, 3):
        raise RuntimeError(f"Reference construction returned {reference.shape}, expected {(steps, 3)}")
    return reference


def continuous_boundary_tour(steps=12000, targets=BOUNDARY_TOUR_TARGETS):
    """Closed x→y→z→x tour matching copyBF's large-range demonstration path.

    Three slow C1-smooth axial legs (no simultaneous multi-axis spinning).
    Default targets sit near the source hard box so chance constraints activate
    while remaining trackable under d=10 + plant-state feedback.
    """
    steps = int(steps)
    targets = np.asarray(targets, dtype=float).reshape(3)
    if steps < 40:
        raise ValueError("At least 40 samples are required for the boundary tour")
    if np.any(targets <= 0.0):
        raise ValueError("All boundary-tour targets must be positive")

    anchors = np.array([
        [targets[0], 0.0, 0.0],
        [0.0, targets[1], 0.0],
        [0.0, 0.0, targets[2]],
        [targets[0], 0.0, 0.0],
    ])
    t = np.linspace(0.0, 3.0, steps, endpoint=True)
    segment = np.minimum(np.floor(t).astype(int), 2)
    local = t - segment
    blend = local * local * (3.0 - 2.0 * local)
    reference = (
        (1.0 - blend[:, None]) * anchors[segment]
        + blend[:, None] * anchors[segment + 1]
    )
    if reference.shape != (steps, 3):
        raise RuntimeError(
            f"Boundary tour construction returned {reference.shape}, "
            f"expected {(steps, 3)}"
        )
    return reference


def _helix_theta_and_z(steps, z_hi, rising_turns=CIRCLE_HELIX_RISING_TURNS,
                       dwell_fraction=0.20):
    """Shared bottom-dwell / rising / top-dwell schedule for helix references.

    Rising phase sweeps ``rising_turns`` revolutions while z goes -z_hi→+z_hi.
    ``dwell_fraction`` is the share of samples at each of bottom and top
    (default 0.20 matches the historical steps//5 layout).
    """
    steps = int(steps)
    z_hi = float(z_hi)
    z_lo = -z_hi
    rising_turns = float(rising_turns)
    dwell_fraction = float(dwell_fraction)
    if rising_turns <= 0.0:
        raise ValueError("rising_turns must be positive")
    if not 0.0 <= dwell_fraction < 0.5:
        raise ValueError("dwell_fraction must satisfy 0 <= f < 0.5")
    bottom_steps = int(round(steps * dwell_fraction))
    top_steps = int(round(steps * dwell_fraction))
    if bottom_steps + top_steps >= steps:
        raise ValueError("dwell fractions leave no rising samples")
    rise_steps = steps - bottom_steps - top_steps
    # Bottom dwell: hold one azimuth while parked at z_lo.
    theta_bottom = np.zeros(bottom_steps)
    theta_rise = np.linspace(
        0.0, 2.0 * np.pi * rising_turns, rise_steps, endpoint=False,
    )
    # Top dwell: hold the final azimuth.
    theta_top = np.full(top_steps, 2.0 * np.pi * rising_turns)
    theta = np.concatenate((theta_bottom, theta_rise, theta_top))
    z = np.concatenate((
        np.full(bottom_steps, z_lo),
        np.linspace(z_lo, z_hi, rise_steps, endpoint=False),
        np.full(top_steps, z_hi),
    ))
    return theta, z, bottom_steps, top_steps, rise_steps


def circle_helix_reference(
        steps=12000, pressures=CIRCLE_HELIX_PRESSURES,
        rising_turns=CIRCLE_HELIX_RISING_TURNS):
    """Bottom dwell, rising circular/elliptical helix, top dwell.

    Horizontal cross-section is a smooth ellipse (circle when rx=ry). Segment
    layout: bottom dwell steps//5 at z=-z_hi, ``rising_turns`` revolutions with
    z linear -z_hi→+z_hi, top dwell steps//5. Default one rising turn keeps
    the path inside the identified plant's trackable band under Q/Ru/N/d=10.
    """
    steps = int(steps)
    pressures = np.asarray(pressures, dtype=float).reshape(3)
    if steps < 40:
        raise ValueError("At least 40 samples are required for the circle helix")
    if np.any(pressures <= 0.0):
        raise ValueError("All pressure magnitudes must be positive")

    rx, ry, z_hi = pressures
    theta, z, _, _, _ = _helix_theta_and_z(steps, z_hi, rising_turns=rising_turns)
    reference = np.column_stack((rx * np.cos(theta), ry * np.sin(theta), z))
    if reference.shape != (steps, 3):
        raise RuntimeError(
            f"Circle helix construction returned {reference.shape}, "
            f"expected {(steps, 3)}"
        )
    return reference


def square_helix_reference(
        steps=12000, pressures=SQUARE_SPIRAL_PRESSURES,
        rising_turns=CIRCLE_HELIX_RISING_TURNS):
    """Bottom dwell, rising square helix, top dwell, without resets.

    Optional sustained-pressure mode: every point satisfies
    |x| <= rx, |y| <= ry with at least one coordinate on the square edge.
    Prefer smooth_square_helix_reference for large-range tracking with
    filleted corners; keep this sharp L∞ square for corner-stress tests.
    """
    steps = int(steps)
    pressures = np.asarray(pressures, dtype=float).reshape(3)
    if steps < 40:
        raise ValueError("At least 40 samples are required for the square helix")
    if np.any(pressures <= 0.0):
        raise ValueError("All pressure magnitudes must be positive")

    rx, ry, z_hi = pressures
    theta, z, _, _, _ = _helix_theta_and_z(steps, z_hi, rising_turns=rising_turns)
    denom = np.maximum(np.abs(np.cos(theta)), np.abs(np.sin(theta)))
    x = rx * np.cos(theta) / denom
    y = ry * np.sin(theta) / denom
    reference = np.column_stack((x, y, z))
    if reference.shape != (steps, 3):
        raise RuntimeError(
            f"Square helix construction returned {reference.shape}, "
            f"expected {(steps, 3)}"
        )
    return reference


def steady_state_task_input(model, E, task):
    """Minimum-norm steady input for a constant standardized task target.

    Experiment diagnostic only. Uses the frozen (A,B,P,E) plant; does not
    change CRTE identification, Sigma_eps, or copyAU preference geometry.
    """
    A = np.asarray(model["A"], dtype=float)
    B = np.asarray(model["B"], dtype=float)
    P = np.asarray(model["P"], dtype=float)
    E = np.asarray(E, dtype=float)
    task = np.asarray(task, dtype=float).reshape(-1)
    dc = E.T @ P @ np.linalg.solve(np.eye(A.shape[0]) - A, B)
    # Minimize the implemented actuator metric ||u||_inf directly.  The
    # pseudoinverse minimizes ||u||_2; with four motors and three tracked
    # coordinates it can report a peak above the box limit even though the
    # motor nullspace contains an exact feasible redistribution.
    n_u = dc.shape[1]
    objective = np.r_[np.zeros(n_u), 1.0]
    peak_constraints = np.block([
        [np.eye(n_u), -np.ones((n_u, 1))],
        [-np.eye(n_u), -np.ones((n_u, 1))],
    ])
    result = linprog(
        objective,
        A_ub=peak_constraints,
        b_ub=np.zeros(2 * n_u),
        A_eq=np.c_[dc, np.zeros(dc.shape[0])],
        b_eq=task,
        bounds=[(None, None)] * n_u + [(0.0, None)],
        method="highs",
    )
    if not result.success:
        raise RuntimeError(f"steady-state minimax input solve failed: {result.message}")
    u = result.x[:n_u]
    y_hat = dc @ u
    return {
        "task": task,
        "u": u,
        "u_inf": float(np.max(np.abs(u))),
        "u_rms": float(np.sqrt(np.mean(u * u))),
        "y_hat": y_hat,
        "task_residual": task - y_hat,
        "dc_gain": dc,
    }


def smooth_square_dc_feasibility(
        model, E, pressures=SMOOTH_SQUARE_HELIX_PRESSURES,
        corner_power=SMOOTH_SQUARE_CORNER_POWER,
        u_inf_cap=SMOOTH_SQUARE_DC_U_INF_CAP,
        n_theta=720):
    """Worst-case steady |u|_inf on the superellipse orbit (incl. corners).

    Returns a dict; ``feasible`` is True when corner/edge steady demand stays
    under ``u_inf_cap`` (default 5.5 < actuator 6, leaving dynamic headroom).
    """
    pressures = np.asarray(pressures, dtype=float).reshape(3)
    rx, ry, z_hi = pressures
    n = float(corner_power)
    exp = 2.0 / n
    theta = np.linspace(0.0, 2.0 * np.pi, int(n_theta), endpoint=False)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    x = rx * np.sign(cos_t) * np.maximum(np.abs(cos_t), 1e-30) ** exp
    y = ry * np.sign(sin_t) * np.maximum(np.abs(sin_t), 1e-30) ** exp
    worst = None
    for xi, yi in zip(x, y):
        for z in (0.0, float(z_hi)):
            info = steady_state_task_input(model, E, np.array([xi, yi, z]))
            if worst is None or info["u_inf"] > worst["u_inf"]:
                worst = info
    if worst is None:
        raise RuntimeError("smooth-square DC scan produced no samples")
    return {
        "pressures": pressures,
        "corner_power": n,
        "u_inf_cap": float(u_inf_cap),
        "worst_u_inf": float(worst["u_inf"]),
        "worst_task": np.asarray(worst["task"], dtype=float),
        "worst_u": np.asarray(worst["u"], dtype=float),
        "feasible": bool(worst["u_inf"] <= float(u_inf_cap)),
        "note": (
            "Steady-state standardized-command demand on the rounded-square orbit. "
            "Infeasible means the reference is not DC-reachable under the implemented ±6 cap; "
            "closed-loop RMSE then reflects saturation/shrink, not CRTE/copyAU."
        ),
    }


def smooth_square_helix_reference(
        steps=12000, pressures=SMOOTH_SQUARE_HELIX_PRESSURES,
        rising_turns=None,
        corner_power=SMOOTH_SQUARE_CORNER_POWER,
        dwell_fraction=None):
    """DC-feasible rounded-square helix (uniform arc-length superellipse).

    Horizontal section is a superellipse
        |x/rx|^n + |y/ry|^n = 1
    sampled at **uniform arc length** (not uniform θ). Uniform θ makes
    |d r / dθ| blow up near the axes for n>2 and overloads the tracker;
    arc-length reparameterization keeps a nearly constant horizontal speed
    while preserving filleted corners (n=2 ellipse → n→∞ sharp square;
    default n=4 = rounded box).

    Default pressures are calibrated so the *corner* steady input stays under
    the actuator bound with headroom (see ``smooth_square_dc_feasibility``).
    Near-hard large-range demos use ``boundary-tour`` (single-axis legs).

    Vertical schedule: bottom dwell, slow rise over ``rising_turns``
    revolutions, top dwell.
    """
    steps = int(steps)
    pressures = np.asarray(pressures, dtype=float).reshape(3)
    n = float(corner_power)
    turns = float(
        SMOOTH_SQUARE_RISING_TURNS if rising_turns is None else rising_turns
    )
    dwell = float(
        SMOOTH_SQUARE_DWELL_FRACTION if dwell_fraction is None else dwell_fraction
    )
    if steps < 40:
        raise ValueError("At least 40 samples are required for the smooth square helix")
    if np.any(pressures <= 0.0):
        raise ValueError("All pressure magnitudes must be positive")
    if n < 2.0:
        raise ValueError("corner_power must be >= 2 (2=ellipse, larger→squarer)")
    if turns <= 0.0:
        raise ValueError("rising_turns must be positive")

    rx, ry, z_hi = pressures
    _, z, bottom_steps, top_steps, rise_steps = _helix_theta_and_z(
        steps, z_hi, rising_turns=turns, dwell_fraction=dwell,
    )

    # Dense superellipse samples over one turn, then arc-length table.
    n_dense = max(4096, 512 * int(np.ceil(turns)))
    theta_dense = np.linspace(0.0, 2.0 * np.pi, n_dense, endpoint=False)
    exp = 2.0 / n
    cos_t = np.cos(theta_dense)
    sin_t = np.sin(theta_dense)
    x_one = rx * np.sign(cos_t) * np.maximum(np.abs(cos_t), 1e-30) ** exp
    y_one = ry * np.sign(sin_t) * np.maximum(np.abs(sin_t), 1e-30) ** exp
    # Close the loop for arc length of one revolution.
    x_closed = np.concatenate((x_one, x_one[:1]))
    y_closed = np.concatenate((y_one, y_one[:1]))
    seg = np.hypot(np.diff(x_closed), np.diff(y_closed))
    s_tab = np.concatenate(([0.0], np.cumsum(seg)))
    perimeter = float(s_tab[-1])
    if perimeter <= 0.0:
        raise RuntimeError("smooth square helix perimeter collapsed")

    def xy_at_s(arc_length):
        s = np.mod(np.asarray(arc_length, dtype=float), perimeter)
        return np.interp(s, s_tab, x_closed), np.interp(s, s_tab, y_closed)

    # Dwell holds arc position; rise advances turns * perimeter at const speed.
    s_bottom = np.zeros(bottom_steps)
    s_rise = np.linspace(0.0, turns * perimeter, rise_steps, endpoint=False)
    s_top = np.full(top_steps, turns * perimeter)
    s_all = np.concatenate((s_bottom, s_rise, s_top))
    if s_all.shape[0] != steps:
        raise RuntimeError(
            f"arc-length schedule length {s_all.shape[0]} != steps {steps}"
        )
    x, y = xy_at_s(s_all)
    reference = np.column_stack((x, y, z))
    if reference.shape != (steps, 3):
        raise RuntimeError(
            f"Smooth square helix construction returned {reference.shape}, "
            f"expected {(steps, 3)}"
        )
    return reference


def _smoothstep01(t):
    """C1-ish ramp on [0,1]: 0 at 0, 1 at 1, zero endpoint derivatives."""
    t = np.clip(np.asarray(t, dtype=float), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _eased_unit_grid(n, ease_frac=0.06):
    """Monotone samples on [0,1) with cosine ease-in/out of along-track speed.

    Removes hard start/stop Δ² spikes when a constant-speed rise is glued to
    zero-speed dwells (dominant reference-side jerk at dwell junctions).
    """
    n = int(n)
    if n < 2:
        return np.zeros(max(n, 0))
    ease_frac = float(np.clip(ease_frac, 0.0, 0.45))
    t = np.linspace(0.0, 1.0, n, endpoint=False)
    v = np.ones(n, dtype=float)
    if ease_frac > 0.0:
        left = t < ease_frac
        right = t > (1.0 - ease_frac)
        if np.any(left):
            v[left] = _smoothstep01(t[left] / ease_frac)
        if np.any(right):
            v[right] = _smoothstep01((1.0 - t[right]) / ease_frac)
        v = np.maximum(v, 1e-6)
    c = np.cumsum(v)
    c = (c - c[0]) / (c[-1] - c[0])
    # endpoint=False style: stay strictly below 1
    return np.minimum(c, 1.0 - 1e-15)


def _unit_superellipse_arc_table(corner_power=SMOOTH_SQUARE_CORNER_POWER, n_dense=8192):
    """Closed unit superellipse |x|^n + |y|^n = 1 as arc-length tables."""
    n = float(corner_power)
    if n < 1.0:
        raise ValueError("corner_power must be >= 1 (1=diamond, 2=ellipse, larger→squarer)")
    n_dense = int(n_dense)
    exp = 2.0 / n
    theta = np.linspace(0.0, 2.0 * np.pi, n_dense, endpoint=False)
    c = np.cos(theta)
    s = np.sin(theta)
    x = np.sign(c) * np.maximum(np.abs(c), 1e-30) ** exp
    y = np.sign(s) * np.maximum(np.abs(s), 1e-30) ** exp
    x_closed = np.concatenate((x, x[:1]))
    y_closed = np.concatenate((y, y[:1]))
    seg = np.hypot(np.diff(x_closed), np.diff(y_closed))
    s_tab = np.concatenate(([0.0], np.cumsum(seg)))
    perimeter = float(s_tab[-1])
    if perimeter <= 0.0:
        raise RuntimeError("unit superellipse perimeter collapsed")
    return {
        "s_tab": s_tab,
        "x_closed": x_closed,
        "y_closed": y_closed,
        "perimeter": perimeter,
        "corner_power": n,
        "family": "superellipse",
    }


def _unit_large_arc_diamond_table(
        center_offset=GROWING_ARC_CENTER_OFFSET, n_dense=8192):
    """Unit loop through axis tips (±1,0)/(0,±1) via four large circular arcs.

    Each arc is centered on the *opposite* diagonal so the path bows inward and
    never visits the multi-axis square corner (1,1). Larger ``center_offset``
    → flatter sides (closer to a diamond) and weaker simultaneous |x|∧|y|.
    """
    a = float(center_offset)
    if a < 0.05:
        raise ValueError("center_offset must be > 0 (larger = stronger corner cut)")
    n_dense = int(n_dense)
    # Four arcs CCW: +x→+y, +y→−x, −x→−y, −y→+x.
    # Centers sit in the opposite quadrant relative to each arc's square corner.
    arcs = (
        # from, to, center
        ((1.0, 0.0), (0.0, 1.0), (-a, -a)),
        ((0.0, 1.0), (-1.0, 0.0), (a, -a)),
        ((-1.0, 0.0), (0.0, -1.0), (a, a)),
        ((0.0, -1.0), (1.0, 0.0), (-a, a)),
    )
    xs = []
    ys = []
    n_per = max(int(n_dense // 4), 64)
    for (x0, y0), (x1, y1), (cx, cy) in arcs:
        v0 = np.array([x0 - cx, y0 - cy], dtype=float)
        v1 = np.array([x1 - cx, y1 - cy], dtype=float)
        r0 = float(np.hypot(v0[0], v0[1]))
        r1 = float(np.hypot(v1[0], v1[1]))
        if abs(r0 - r1) > 1e-9 * max(r0, 1.0):
            raise RuntimeError("large-arc diamond radii mismatch")
        th0 = float(np.arctan2(v0[1], v0[0]))
        th1 = float(np.arctan2(v1[1], v1[0]))
        # CCW delta in (-pi, pi] then force positive CCW sweep.
        dth = th1 - th0
        while dth <= 0.0:
            dth += 2.0 * np.pi
        while dth > 2.0 * np.pi:
            dth -= 2.0 * np.pi
        th = th0 + np.linspace(0.0, dth, n_per, endpoint=False)
        xs.append(cx + r0 * np.cos(th))
        ys.append(cy + r0 * np.sin(th))
    x = np.concatenate(xs)
    y = np.concatenate(ys)
    # Renormalize so axis tips land exactly at unit axes (numerical guard).
    ax = float(np.max(np.abs(x)))
    ay = float(np.max(np.abs(y)))
    if ax < 1e-12 or ay < 1e-12:
        raise RuntimeError("large-arc diamond collapsed")
    x = x / ax
    y = y / ay
    x_closed = np.concatenate((x, x[:1]))
    y_closed = np.concatenate((y, y[:1]))
    seg = np.hypot(np.diff(x_closed), np.diff(y_closed))
    s_tab = np.concatenate(([0.0], np.cumsum(seg)))
    perimeter = float(s_tab[-1])
    if perimeter <= 0.0:
        raise RuntimeError("large-arc diamond perimeter collapsed")
    # Simultaneous multi-axis peak on the path (should be << 1).
    max_joint = float(np.max(np.minimum(np.abs(x), np.abs(y))))
    return {
        "s_tab": s_tab,
        "x_closed": x_closed,
        "y_closed": y_closed,
        "perimeter": perimeter,
        "corner_power": None,
        "family": "large-arc-diamond",
        "center_offset": a,
        "max_joint_xy": max_joint,
    }


def _unit_path_table(family=None, corner_power=SMOOTH_SQUARE_CORNER_POWER,
                     center_offset=None, n_dense=8192):
    """Dispatch unit horizontal path table for growing / smooth envelopes."""
    fam = GROWING_PATH_FAMILY if family is None else str(family)
    if fam in ("large-arc-diamond", "large_arc_diamond", "axial-arc", "diamond-arc"):
        off = GROWING_ARC_CENTER_OFFSET if center_offset is None else center_offset
        return _unit_large_arc_diamond_table(
            center_offset=off, n_dense=n_dense,
        )
    if fam in ("superellipse", "smooth-square", "box"):
        return _unit_superellipse_arc_table(
            corner_power=corner_power, n_dense=n_dense,
        )
    raise ValueError(f"unknown path family: {fam}")


def _unit_xy_at_arc(arc_length, table):
    s = np.mod(np.asarray(arc_length, dtype=float), table["perimeter"])
    x = np.interp(s, table["s_tab"], table["x_closed"])
    y = np.interp(s, table["s_tab"], table["y_closed"])
    return x, y


def _axis_tip_speed_weights(x, y, slowdown=GROWING_AXIS_TIP_SLOWDOWN):
    """Lower along-track speed near axis tips (high |x| or |y| peaks).

    ``slowdown`` in [0, 0.95]: weight = 1 - slowdown * tipness, tipness→1 at
    pure axis extremes of the current sample set.
    """
    slow = float(np.clip(slowdown, 0.0, 0.95))
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ax = np.max(np.abs(x)) + 1e-15
    ay = np.max(np.abs(y)) + 1e-15
    # tipness high when one axis is near its max and the other is near 0.
    nx = np.abs(x) / ax
    ny = np.abs(y) / ay
    tip = np.maximum(nx * (1.0 - ny), ny * (1.0 - nx))
    tip = np.clip(tip, 0.0, 1.0)
    return np.maximum(1.0 - slow * tip, 1e-3)


def large_arc_diamond_dc_feasibility(
        model, E, pressures=None,
        center_offset=None,
        u_inf_cap=None,
        n_samples=720):
    """Worst-case steady |u|_inf on the large-arc axial diamond orbit."""
    pressures = np.asarray(
        GROWING_SMOOTH_SQUARE_PEAK if pressures is None else pressures,
        dtype=float,
    ).reshape(3)
    off = GROWING_ARC_CENTER_OFFSET if center_offset is None else float(center_offset)
    cap = GROWING_DC_U_INF_CAP if u_inf_cap is None else float(u_inf_cap)
    table = _unit_large_arc_diamond_table(center_offset=off, n_dense=max(n_samples * 4, 2048))
    # Uniform arc samples on one rev.
    s = np.linspace(0.0, table["perimeter"], int(n_samples), endpoint=False)
    ux, uy = _unit_xy_at_arc(s, table)
    rx, ry, z_hi = pressures
    worst = None
    for xi, yi in zip(rx * ux, ry * uy):
        for z in (0.0, float(z_hi), -float(z_hi)):
            info = steady_state_task_input(model, E, np.array([xi, yi, z]))
            if worst is None or info["u_inf"] > worst["u_inf"]:
                worst = info
    if worst is None:
        raise RuntimeError("large-arc DC scan produced no samples")
    return {
        "pressures": pressures,
        "family": "large-arc-diamond",
        "center_offset": off,
        "max_joint_xy_unit": float(table["max_joint_xy"]),
        "u_inf_cap": float(cap),
        "worst_u_inf": float(worst["u_inf"]),
        "worst_task": np.asarray(worst["task"], dtype=float),
        "worst_u": np.asarray(worst["u"], dtype=float),
        "feasible": bool(worst["u_inf"] <= float(cap)),
        "note": (
            "Steady actuator demand on large-arc axial diamond (no square "
            "corners). Cap is diagnostic (boundary-tour band ~7); ±6 may be "
            "briefly saturated in closed loop while chance QP still binds."
        ),
    }


def z_face_large_arc_cruise_dc_feasibility(
        model, E, xy_radius=None, z_face=None, u_inf_cap=None, n_samples=720):
    """Worst steady |u|_inf on both moving z-face large-arc cruises."""
    peak = np.asarray(GROWING_SMOOTH_SQUARE_PEAK, dtype=float)
    radius = float(peak[0] if xy_radius is None else xy_radius)
    z_abs = float(peak[2] if z_face is None else z_face)
    cap = GROWING_DC_U_INF_CAP if u_inf_cap is None else float(u_inf_cap)
    table = _unit_large_arc_diamond_table(
        center_offset=GROWING_ARC_CENTER_OFFSET,
        n_dense=max(int(n_samples) * 4, 2048),
    )
    s = np.linspace(0.0, table["perimeter"], int(n_samples), endpoint=False)
    ux, uy = _unit_xy_at_arc(s, table)
    worst = None
    for z in (z_abs, -z_abs):
        for xi, yi in zip(radius * ux, radius * uy):
            info = steady_state_task_input(model, E, np.array([xi, yi, z]))
            if worst is None or info["u_inf"] > worst["u_inf"]:
                worst = info
    if worst is None:
        raise RuntimeError("z-face cruise DC scan produced no samples")
    return {
        "pressures": np.array([radius, radius, z_abs]),
        "family": "z-face-large-arc-cruise",
        "center_offset": float(GROWING_ARC_CENTER_OFFSET),
        "max_joint_xy_unit": float(table["max_joint_xy"]),
        "u_inf_cap": float(cap),
        "worst_u_inf": float(worst["u_inf"]),
        "worst_task": np.asarray(worst["task"], dtype=float),
        "worst_u": np.asarray(worst["u"], dtype=float),
        "feasible": bool(worst["u_inf"] <= cap),
        "note": (
            "Steady standardized-command demand on both z=±face cruises, including the "
            "moving horizontal large-arc orbit; checked against the inherited ±6 "
            "experiment-design cap (not a Pelican hardware rating)."
        ),
    }


def _scaled_superellipse_helix_segment(
        steps,
        pressures_a,
        pressures_b,
        turns,
        corner_power=SMOOTH_SQUARE_CORNER_POWER,
        z_from=None,
        z_to=None,
        arc_offset=0.0,
        ease_frac=0.05):
    """Arc-length paced superellipse helix from box A → box B (one leg)."""
    steps = int(steps)
    a = np.asarray(pressures_a, dtype=float).reshape(3)
    b = np.asarray(pressures_b, dtype=float).reshape(3)
    turns = float(turns)
    if steps < 2:
        raise ValueError("helix segment needs at least 2 samples")
    if turns < 0.0:
        raise ValueError("turns must be non-negative")
    table = _unit_superellipse_arc_table(corner_power=corner_power)
    perimeter = table["perimeter"]
    n_dense = max(4096, int(1024 * max(turns, 1.0)))
    tau = np.linspace(0.0, 1.0, n_dense)
    if np.allclose(a, b):
        alpha = np.ones_like(tau)
        scale = np.repeat(a[None, :], n_dense, axis=0)
    else:
        alpha = _smoothstep01(tau)
        scale = a[None, :] + alpha[:, None] * (b - a)[None, :]
    unit_arc = float(arc_offset) + tau * turns * perimeter
    ux, uy = _unit_xy_at_arc(unit_arc, table)
    z0 = float(a[2] if z_from is None else z_from)
    z1 = float(b[2] if z_to is None else z_to)
    z = z0 + _smoothstep01(tau) * (z1 - z0) if not np.isclose(z0, z1) else np.full(n_dense, z0)
    # When scale is constant, still allow linear z via smoothstep
    if np.allclose(a, b) and not np.isclose(z0, z1):
        z = z0 + _smoothstep01(tau) * (z1 - z0)
    xyz = np.column_stack((scale[:, 0] * ux, scale[:, 1] * uy, z))
    seg = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    s_tab = np.concatenate(([0.0], np.cumsum(seg)))
    total = float(s_tab[-1])
    if total <= 1e-12:
        # Degenerate (no motion): hold point.
        return np.repeat(xyz[:1], steps, axis=0)
    s_query = _eased_unit_grid(steps, ease_frac=float(ease_frac)) * total
    out = np.column_stack([
        np.interp(s_query, s_tab, xyz[:, axis]) for axis in range(3)
    ])
    return out


def _allocate_phase_steps(total, fractions):
    fractions = np.asarray(fractions, dtype=float).reshape(-1)
    if fractions.size == 0:
        raise ValueError("fractions empty")
    if np.any(fractions < 0.0):
        raise ValueError("fractions must be non-negative")
    s = float(np.sum(fractions))
    if s <= 0.0:
        raise ValueError("fractions sum must be positive")
    fractions = fractions / s
    raw = fractions * int(total)
    steps = np.floor(raw).astype(int)
    for i, f in enumerate(fractions):
        if f > 0.0 and steps[i] < 2 and int(total) >= 2 * len(fractions):
            steps[i] = 2
    deficit = int(total) - int(np.sum(steps))
    order = np.argsort(-(raw - steps))
    k = 0
    while deficit > 0:
        steps[order[k % len(steps)]] += 1
        deficit -= 1
        k += 1
    while deficit < 0:
        j = int(np.argmax(steps))
        if steps[j] <= 2:
            break
        steps[j] -= 1
        deficit += 1
    if int(np.sum(steps)) != int(total):
        steps[-1] += int(total) - int(np.sum(steps))
    return steps.astype(int)


def z_face_large_arc_cruise_reference(
        steps,
        pressures_start=GROWING_SMOOTH_SQUARE_START,
        pressures_peak=GROWING_SMOOTH_SQUARE_PEAK,
        phase_fractions=GROWING_PHASE_FRACTIONS,
        phase_turns=GROWING_Z_FACE_TURNS,
        center_offset=GROWING_ARC_CENTER_OFFSET):
    """Five-phase moving cruise on the upper/lower z task-hard faces."""
    steps = int(steps)
    start = np.asarray(pressures_start, dtype=float).reshape(3)
    peak = np.asarray(pressures_peak, dtype=float).reshape(3)
    turns = np.asarray(phase_turns, dtype=float).reshape(-1)
    counts = _allocate_phase_steps(steps, phase_fractions)
    if counts.size != 5 or turns.size != 5 or np.any(turns <= 0.0):
        raise ValueError("z-face cruise needs five positive phases and turn counts")
    if peak[2] <= start[2] or np.any(start <= 0.0) or np.any(peak <= 0.0):
        raise ValueError("z face must exceed a positive start condition")

    table = _unit_large_arc_diamond_table(
        center_offset=float(center_offset), n_dense=8192,
    )
    perimeter = float(table["perimeter"])
    reference = np.zeros((steps, 3), dtype=float)
    phase_ids = np.repeat(np.arange(5, dtype=int), counts)
    start_turns = np.concatenate(([0.0], np.cumsum(turns[:-1])))
    for phase_id, count in enumerate(counts):
        idx = np.flatnonzero(phase_ids == phase_id)
        local = np.arange(int(count), dtype=float) / float(max(1, int(count)))
        smooth = _smoothstep01(local)
        unit_arc = (start_turns[phase_id] + turns[phase_id] * local) * perimeter
        ux, uy = _unit_xy_at_arc(unit_arc, table)

        if phase_id == 0:  # approach upper face
            rx = start[0] + (peak[0] - start[0]) * smooth
            ry = start[1] + (peak[1] - start[1]) * smooth
            z = start[2] + (peak[2] - start[2]) * smooth
        elif phase_id == 1:  # upper-face cruise
            rx = np.full(count, peak[0])
            ry = np.full(count, peak[1])
            z = np.full(count, peak[2])
        elif phase_id == 2:  # upper→lower transfer
            rx = np.full(count, peak[0])
            ry = np.full(count, peak[1])
            z = peak[2] - 2.0 * peak[2] * smooth
        elif phase_id == 3:  # lower-face cruise
            rx = np.full(count, peak[0])
            ry = np.full(count, peak[1])
            z = np.full(count, -peak[2])
        else:  # return to the initial z/radius condition
            rx = peak[0] + (start[0] - peak[0]) * smooth
            ry = peak[1] + (start[1] - peak[1]) * smooth
            z = -peak[2] + (start[2] + peak[2]) * smooth

        reference[idx, 0] = rx * ux
        reference[idx, 1] = ry * uy
        reference[idx, 2] = z
    return reference


def growing_smooth_square_helix_reference(
        steps=None,
        pressures_end=None,
        pressures_start=None,
        rising_turns=None,
        peak_turns=None,
        corner_power=SMOOTH_SQUARE_CORNER_POWER,
        dwell_fraction=None,
        phase_fractions=None,
        path_family=None,
        center_offset=None,
        axis_tip_slowdown=None):
    """Generate the selected continuous near-hard reference family.

    v4.0 defaults to moving large-arc cruises on z=±3.57. Legacy v3.9
    large-arc-diamond and superellipse families remain available explicitly.
    """
    steps = int(GROWING_SMOOTH_SQUARE_DEFAULT_STEPS if steps is None else steps)
    peak = np.asarray(
        GROWING_SMOOTH_SQUARE_PEAK if pressures_end is None else pressures_end,
        dtype=float,
    ).reshape(3)
    start = np.asarray(
        GROWING_SMOOTH_SQUARE_START if pressures_start is None else pressures_start,
        dtype=float,
    ).reshape(3)
    turns_g = float(
        GROWING_SMOOTH_SQUARE_RISING_TURNS if rising_turns is None else rising_turns
    )
    turns_p = float(
        GROWING_SMOOTH_SQUARE_PEAK_TURNS if peak_turns is None else peak_turns
    )
    fracs = GROWING_PHASE_FRACTIONS if phase_fractions is None else phase_fractions
    _ = dwell_fraction  # legacy kw ignored (continuous envelope)
    fam = GROWING_PATH_FAMILY if path_family is None else str(path_family)
    off = GROWING_ARC_CENTER_OFFSET if center_offset is None else float(center_offset)
    tip_slow = (
        GROWING_AXIS_TIP_SLOWDOWN if axis_tip_slowdown is None
        else float(axis_tip_slowdown)
    )
    n = float(corner_power)
    if steps < 200:
        raise ValueError("growing envelope needs at least 200 samples")
    if np.any(start <= 0.0) or np.any(peak <= 0.0):
        raise ValueError("start/peak pressures must be positive")
    if np.any(peak < start):
        raise ValueError("peak pressures must be componentwise >= start")
    if turns_g <= 0.0 or turns_p <= 0.0:
        raise ValueError("turns must be positive")

    if fam in ("z-face-large-arc-cruise", "z_face_large_arc_cruise", "z-face"):
        z_turns = (
            turns_g,
            turns_p,
            float(GROWING_Z_FACE_TURNS[2]),
            turns_p,
            turns_g,
        )
        return z_face_large_arc_cruise_reference(
            steps=steps,
            pressures_start=start,
            pressures_peak=peak,
            phase_fractions=fracs,
            phase_turns=z_turns,
            center_offset=off,
        )

    n_grow, n_peak, n_shrink = _allocate_phase_steps(steps, fracs)
    table = _unit_path_table(
        family=fam, corner_power=n, center_offset=off,
    )
    perimeter = float(table["perimeter"])

    # Dense continuous geometry (one path) then global arc-length resample so
    # phase junctions don't create speed spikes.
    dens_counts = []
    for n_seg, turns_seg in (
        (max(n_grow * 4, 2048), turns_g),
        (max(n_peak * 4, 2048), turns_p),
        (max(n_shrink * 4, 2048), turns_g),
    ):
        dens_counts.append(int(n_seg))
    n_dg, n_dp, n_ds = dens_counts
    # scale alpha: 0→1 grow, 1 hold peak, 1→0 shrink
    a_g = _smoothstep01(np.linspace(0.0, 1.0, n_dg, endpoint=False))
    a_p = np.ones(n_dp)
    a_s = 1.0 - _smoothstep01(np.linspace(0.0, 1.0, n_ds, endpoint=False))
    alpha = np.concatenate((a_g, a_p, a_s))
    scale = start[None, :] + alpha[:, None] * (peak - start)[None, :]
    # continuous unit-arc angle
    turns_total = turns_g + turns_p + turns_g
    unit_arc = np.linspace(0.0, turns_total * perimeter, alpha.size, endpoint=False)
    ux, uy = _unit_xy_at_arc(unit_arc, table)
    # z: climb during grow, swing +peak→−peak during peak loiter (3D pressure),
    # climb down during shrink back to −start_z
    z_g = -start[2] + a_g * (peak[2] + start[2])
    z_p = peak[2] + (-peak[2] - peak[2]) * _smoothstep01(
        np.linspace(0.0, 1.0, n_dp, endpoint=False)
    )
    z_s = z_p[-1] + (-start[2] - z_p[-1]) * _smoothstep01(
        np.linspace(0.0, 1.0, n_ds, endpoint=False)
    )
    z = np.concatenate((z_g, z_p, z_s))
    xyz = np.column_stack((scale[:, 0] * ux, scale[:, 1] * uy, z))
    # Tip slowdown: densify samples near axis extremes via weighted arc length
    # so the final time-resample spends more dwell near axial peaks.
    raw_seg = np.linalg.norm(np.diff(xyz, axis=0), axis=1)
    w_pt = _axis_tip_speed_weights(xyz[:, 0], xyz[:, 1], slowdown=tip_slow)
    w_seg = 0.5 * (w_pt[:-1] + w_pt[1:])
    # Lower speed at tips → stretch effective arc budget there (weight = 1/speed).
    seg = raw_seg / np.maximum(w_seg, 1e-3)
    s_tab = np.concatenate(([0.0], np.cumsum(seg)))
    total = float(s_tab[-1])
    if total <= 1e-12:
        raise RuntimeError("growing large-arc envelope path length collapsed")
    # Map final samples: use cumulative arc budget proportional to phase
    # fractions so grow/peak/shrink keep intended time shares after resample.
    dens_arc = s_tab
    # phase boundaries in dense index
    i0, i1, i2 = 0, n_dg, n_dg + n_dp
    s0 = 0.0
    s1 = float(dens_arc[i1])
    s2 = float(dens_arc[i2])
    s3 = total
    s_query = np.concatenate((
        s0 + _eased_unit_grid(n_grow, ease_frac=0.04) * (s1 - s0),
        s1 + _eased_unit_grid(n_peak, ease_frac=0.02) * (s2 - s1),
        s2 + _eased_unit_grid(n_shrink, ease_frac=0.04) * (s3 - s2),
    ))
    # Keep strictly increasing for interp stability
    s_query = np.maximum.accumulate(s_query)
    s_query[-1] = min(s_query[-1], total)
    reference = np.column_stack([
        np.interp(s_query, s_tab, xyz[:, axis]) for axis in range(3)
    ])
    if reference.shape != (steps, 3):
        raise RuntimeError(
            f"growing large-arc envelope shape {reference.shape}, expected {(steps, 3)}"
        )
    return reference


def growing_envelope_phase_index(steps, phase_fractions=None):
    """Per-sample phase id: 0 grow, 1 peak-loiter, 2 shrink."""
    steps = int(steps)
    fracs = GROWING_PHASE_FRACTIONS if phase_fractions is None else phase_fractions
    counts = _allocate_phase_steps(steps, fracs)
    ids = np.concatenate([
        np.full(int(c), i, dtype=int) for i, c in enumerate(counts)
    ])
    if ids.size < steps:
        ids = np.concatenate((ids, np.full(steps - ids.size, len(counts) - 1, dtype=int)))
    return ids[:steps], counts


def growing_smooth_square_phase_masks(
        reference, r_dc=None, hard_bound=None, phase_ids=None,
        peak_pressures=None):
    """Label samples by envelope chapter + DC/hard proximity."""
    reference = np.asarray(reference, dtype=float)
    r_dc = float(GROWING_DC_LABEL_R if r_dc is None else r_dc)
    r_box = np.maximum(np.abs(reference[:, 0]), np.abs(reference[:, 1]))
    r_xy = np.hypot(reference[:, 0], reference[:, 1])
    r_linf = np.max(np.abs(reference), axis=1)
    dc = r_box <= r_dc + 1e-6
    peak = np.asarray(
        GROWING_SMOOTH_SQUARE_PEAK if peak_pressures is None else peak_pressures,
        dtype=float,
    ).reshape(3)
    out = {
        "r_xy": r_xy,
        "r_box": r_box,
        "r_linf": r_linf,
        "dc_feasible_mask": dc,
        "beyond_dc_mask": ~dc,
        "r_xy_max": float(np.max(r_xy)),
        "r_box_max": float(np.max(r_box)),
        "r_linf_max": float(np.max(r_linf)),
        "r_dc": r_dc,
        "peak_pressures": peak,
    }
    if hard_bound is not None:
        hard = float(hard_bound)
        out["near_hard_mask"] = (
            (r_linf <= hard)
            & (r_linf >= hard - GROWING_NEAR_HARD_BAND)
        )
        out["fraction_near_hard"] = float(np.mean(out["near_hard_mask"]))
    if phase_ids is not None:
        phase_ids = np.asarray(phase_ids, dtype=int).reshape(-1)
        if int(np.max(phase_ids, initial=-1)) == 4:
            names = (
                "approach_top", "top_cruise", "face_transfer",
                "bottom_cruise", "return",
            )
        else:
            names = ("grow", "peak_square", "shrink")
        out["phase_ids"] = phase_ids
        out["phase_names"] = names
        for i, name in enumerate(names):
            m = phase_ids == i
            out[f"fraction_{name}"] = float(np.mean(m)) if phase_ids.size else 0.0
        out["fraction_smpc_peak"] = (
            out.get("fraction_peak_square", 0.0)
            + out.get("fraction_top_cruise", 0.0)
            + out.get("fraction_bottom_cruise", 0.0)
        )
    return out


def near_hard_occupancy(trajectory, hard_bound, band=GROWING_NEAR_HARD_BAND):
    """Measure time inside, and within ``band`` of, any task hard face."""
    values = np.asarray(trajectory, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError("trajectory must have shape (N, 3)")
    hard = float(hard_bound)
    width = float(band)
    if hard <= 0.0 or not 0.0 < width < hard:
        raise ValueError("hard_bound and band must satisfy 0 < band < hard_bound")
    linf = np.max(np.abs(values), axis=1)
    inside = linf <= hard
    near = inside & (linf >= hard - width)
    violations = linf > hard
    return {
        "hard_bound": hard,
        "band": width,
        "threshold": hard - width,
        "fraction": float(np.mean(near)),
        "steps": int(np.sum(near)),
        "violation_steps": int(np.sum(violations)),
        "min_hard_margin": float(np.min(hard - linf)),
        "max_abs": float(np.max(linf)),
    }


def future_reference_horizon(reference, step, horizon, control_interval=CONTROL_INTERVAL_STEPS):
    """Reference samples aligned with the N lifted (d-step) predicted states."""
    reference = np.asarray(reference, dtype=float)
    d = int(control_interval)
    indices = np.minimum(
        int(step) + d * np.arange(1, int(horizon) + 1), len(reference) - 1,
    )
    return reference[indices]


def mstep_lift(A, B, d):
    """Lift single-step (A,B) to a d-step decision model with held input."""
    d = int(d)
    Ad = np.linalg.matrix_power(A, d)
    Bd = np.zeros_like(B)
    Ai = np.eye(A.shape[0])
    for _ in range(d):
        Bd = Bd + Ai @ B
        Ai = Ai @ A
    return Ad, Bd


def lift_noise_covariance(A, Sigma, d):
    """Accumulate single-step process noise over a d-step control interval."""
    d = int(d)
    Sigma_d = np.zeros_like(Sigma)
    Ai = np.eye(A.shape[0])
    for _ in range(d):
        Sigma_d = Sigma_d + Ai @ Sigma @ Ai.T
        Ai = Ai @ A
    return _symmetrize_psd(Sigma_d)


def equilibrium_initial_condition(model, E, target):
    """Solve the joint second-order state/input equilibrium at a task target."""
    n = model["A"].shape[0]
    nu = model["B"].shape[1]
    target = np.asarray(target, dtype=float).reshape(3, 1)
    equations = np.block([
        [np.eye(n) - model["A"], -model["B"]],
        [E.T @ model["P"], np.zeros((3, nu))],
    ])
    rhs = np.vstack([
        -model["B"] @ model["u_mean"],
        target - E.T @ model["y_mean"],
    ])
    solution = np.linalg.lstsq(equations, rhs, rcond=None)[0]
    basis = null_space(equations)
    if basis.shape[1] == 1:
        direction = basis[n:, 0]
        lower, upper = -np.inf, np.inf
        for value, slope in zip(solution[n:, 0], direction):
            if abs(slope) < 1e-14:
                if abs(value) > 6.0 + 1e-9:
                    raise ValueError("No bounded-input equilibrium")
                continue
            interval = sorted(((-6.0 - value) / slope, (6.0 - value) / slope))
            lower, upper = max(lower, interval[0]), min(upper, interval[1])
        if lower > upper:
            raise ValueError(f"No bounded-input equilibrium for task target {target.ravel()}")
        solution += basis[:, :1] * np.clip(0.0, lower, upper)
    elif basis.shape[1] > 1:
        raise ValueError("Equilibrium null space has unsupported dimension greater than one")

    state = solution[:n]
    control = solution[n:].ravel()
    dynamic_residual = np.max(np.abs(
        (np.eye(n) - model["A"]) @ state
        - model["B"] @ (control.reshape(-1, 1) - model["u_mean"])
    ))
    task_residual = np.max(np.abs(
        E.T @ (model["y_mean"] + model["P"] @ state) - target
    ))
    if dynamic_residual > 1e-7 or task_residual > 1e-7:
        raise ValueError(
            f"Equilibrium residual too large: dynamics={dynamic_residual:.3e}, "
            f"task={task_residual:.3e}"
        )
    if np.any(np.abs(control) > 6.0 + 1e-8):
        raise ValueError(f"Equilibrium input exceeds bounds for task target {target.ravel()}")
    return state, control


def initialization_diagnostics(model, E, target):
    state, control = equilibrium_initial_condition(model, E, target)
    successor = model["A"] @ state + model["B"] @ (
        control.reshape(-1, 1) - model["u_mean"]
    )
    task = E.T @ (model["y_mean"] + model["P"] @ state)
    target = np.asarray(target, dtype=float).reshape(3, 1)
    return {
        "state": state,
        "control": control,
        "dynamic_residual": float(np.max(np.abs(successor - state))),
        "task_residual": float(np.max(np.abs(task - target))),
        "maximum_absolute_input": float(np.max(np.abs(control))),
    }


class CachedBoundaryController:
    """Full-horizon parameterized QP on the d-step lifted plant (d=10)."""

    def __init__(self, model, E, task_bound, process_covariance,
                 horizon=HORIZON_STEPS, q_weight=Q_WEIGHT, ru=RU,
                 alpha_joint=ALPHA_JOINT,
                 control_interval=CONTROL_INTERVAL_STEPS):
        self.model = model
        self.E = np.asarray(E, dtype=float)
        self.N = int(horizon)
        self.d = int(control_interval)
        self.nu = model["B"].shape[1]
        self.n = model["A"].shape[0]
        self.task_axes = self.E.shape[1]
        self.task_bound = float(task_bound)
        self.alpha_joint = alpha_joint
        self.u_min = np.full(self.nu, -INPUT_COMMAND_BOUND_STANDARDIZED)
        self.u_max = np.full(self.nu, INPUT_COMMAND_BOUND_STANDARDIZED)
        self.u_mean = model["u_mean"].ravel()
        self.y_mean = model["y_mean"].ravel()
        self.Ad, self.Bd = mstep_lift(model["A"], model["B"], self.d)
        # Chance-constraint noise must be d-step accumulated (see
        # predvarx-mpc/references/real-model-in-the-loop.md §5).
        process_covariance = lift_noise_covariance(
            model["A"], _symmetrize_psd(process_covariance), self.d,
        )
        covariance = np.zeros((self.n, self.n))
        stage_tightening = []
        quantile = 0.0
        if alpha_joint is not None:
            per_face = float(alpha_joint) / (2.0 * self.task_axes * self.N)
            quantile = float(norm.ppf(1.0 - per_face))
        for _ in range(self.N):
            covariance = (
                self.Ad @ covariance @ self.Ad.T + process_covariance
            )
            covariance = _symmetrize_psd(covariance)
            task_covariance = (
                self.E.T @ model["P"] @ covariance @ model["P"].T @ self.E
            )
            variances = np.maximum(np.diag(task_covariance), 1e-15)
            stage_tightening.append(quantile * np.sqrt(variances))
        self.stage_tightening = np.asarray(stage_tightening)
        self.stage_mean_bound = self.task_bound - self.stage_tightening
        if np.any(self.stage_mean_bound <= 0.0):
            raise ValueError("Chance tightening leaves an empty task interval")

        self.U = cp.Variable((self.N, self.nu), name="control_plan")
        self.state_parameter = cp.Parameter(self.n, name="current_latent_state")
        self.reference_parameter = cp.Parameter(
            (self.N, self.task_axes), name="future_task_reference",
        )
        constraints = [self.U >= self.u_min, self.U <= self.u_max]
        objective = float(ru) * cp.sum_squares(
            self.U - np.tile(self.u_mean, (self.N, 1))
        )
        state = self.state_parameter
        self.task_expressions = []
        self.task_constraints = []
        for stage in range(self.N):
            state = (
                self.Ad @ state
                + self.Bd @ (self.U[stage] - self.u_mean)
            )
            task = self.E.T @ (self.y_mean + model["P"] @ state)
            self.task_expressions.append(task)
            objective += float(q_weight) * cp.sum_squares(
                task - self.reference_parameter[stage]
            )
            bound = self.stage_mean_bound[stage]
            upper = task <= bound
            lower = task >= -bound
            constraints.extend([upper, lower])
            self.task_constraints.append((upper, lower))
        self.problem = cp.Problem(cp.Minimize(objective), constraints)
        if not self.problem.is_dpp():
            raise RuntimeError("Cached control problem must satisfy DPP")

    def step(self, state, reference_horizon):
        state = np.asarray(state, dtype=float).reshape(self.n)
        reference_horizon = np.asarray(reference_horizon, dtype=float)
        expected = (self.N, self.task_axes)
        if reference_horizon.shape != expected:
            raise ValueError(f"Reference horizon must have shape {expected}")
        self.state_parameter.value = state
        self.reference_parameter.value = reference_horizon
        try:
            self.problem.solve(
                solver=cp.OSQP, warm_start=True, verbose=False,
                eps_abs=1e-8, eps_rel=1e-8, max_iter=20000, polishing=True,
            )
        except cp.error.SolverError as exc:
            return None, {"status": f"solver_error: {exc}"}
        if self.problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
            return None, {"status": self.problem.status}

        planned = np.asarray(self.U.value, dtype=float)
        predicted_task = np.vstack([
            np.asarray(expression.value, dtype=float).reshape(1, -1)
            for expression in self.task_expressions
        ])
        slacks = self.stage_mean_bound - np.abs(predicted_task)
        dual_values = []
        for upper, lower in self.task_constraints:
            dual_values.extend(np.asarray(upper.dual_value, dtype=float).ravel())
            dual_values.extend(np.asarray(lower.dual_value, dtype=float).ravel())
        input_residual = max(
            float(np.max(self.u_min - planned)),
            float(np.max(planned - self.u_max)),
        )
        task_residual = float(np.max(
            np.abs(predicted_task) - self.stage_mean_bound
        ))
        return np.clip(planned[0], self.u_min, self.u_max), {
            "status": self.problem.status,
            "planned_controls": planned,
            "predicted_task": predicted_task,
            "minimum_predicted_slack": float(np.min(slacks)),
            "maximum_constraint_dual": float(np.max(dual_values)),
            "maximum_constraint_residual": max(input_residual, task_residual),
            "solve_time_seconds": float(self.problem.solver_stats.solve_time or 0.0),
        }


def _disturbance_sha256(residual_innovation, wind_latent_increment):
    digest = hashlib.sha256()
    for value in (residual_innovation, wind_latent_increment):
        array = np.ascontiguousarray(np.asarray(value, dtype=np.float64))
        digest.update(array.tobytes())
    return digest.hexdigest()


def observed_augmented_state(model, output, previous_latent):
    """Construct [z_k; z_(k-1)] using outputs available to the controller."""
    output = np.asarray(output, dtype=float).reshape(-1, 1)
    previous_latent = np.asarray(previous_latent, dtype=float).reshape(-1, 1)
    current_latent = model["R"].T @ (output - model["y_mean"])
    if current_latent.shape != previous_latent.shape:
        raise ValueError("Current and previous latent blocks must have equal shapes")
    return np.vstack((current_latent, previous_latent))


def simulate_controller(model, E, hard_bound, reference, residual_innovation,
                        wind_latent_increment, process_covariance, alpha_joint,
                        horizon=HORIZON_STEPS, q_weight=Q_WEIGHT, ru=RU,
                        control_interval=CONTROL_INTERVAL_STEPS,
                        use_plant_state=True, initial_task=None):
    """Simulate one controller on the frozen plant.

    Default ``use_plant_state=True`` matches copyBE/BF model-in-the-loop:
    the QP sees the true augmented plant state (not R.T output reconstruction).
    Reconstruction remains available for output-only diagnostics.

    ``initial_task`` defaults to ``reference[0]``. For large boundary-tour
    targets the input-bounded equilibrium may not exist at the first anchor;
    callers should pass a feasible start (copyBF uses terminal mean bound − ε).
    """
    reference = np.asarray(reference, dtype=float)
    residual_innovation = np.asarray(residual_innovation, dtype=float)
    wind_latent_increment = np.asarray(wind_latent_increment, dtype=float)
    steps = len(reference)
    d = int(control_interval)
    ell = model["A"].shape[0] // 2
    if reference.shape != (steps, 3):
        raise ValueError("reference must have shape (steps, 3)")
    if residual_innovation.shape != (steps, ell):
        raise ValueError("residual innovation shape does not match the latent model")
    if wind_latent_increment.shape != (steps, ell):
        raise ValueError("wind latent increment shape does not match the latent model")

    controller = CachedBoundaryController(
        model, E, hard_bound, process_covariance, horizon=horizon,
        q_weight=q_weight, ru=ru, alpha_joint=alpha_joint,
        control_interval=d,
    )
    if initial_task is None:
        initial_task = reference[0]
    initial_task = np.asarray(initial_task, dtype=float).reshape(3)
    state, control = equilibrium_initial_condition(model, E, initial_task)
    previous_latent = state[ell:].copy()
    task_history, control_history = [], []
    slacks, duals, residuals, solve_times = [], [], [], []
    status_counts = {}
    qp_count = 0
    qp_failure_step = None
    qp_failure_status = ""
    for step in range(steps):
        output = model["y_mean"] + model["P"] @ state
        task = (E.T @ output).ravel()
        if use_plant_state:
            controller_state = state.copy()
        else:
            controller_state = observed_augmented_state(
                model, output, previous_latent,
            )
        if step % d == 0:
            refs = future_reference_horizon(
                reference, step, horizon, control_interval=d,
            )
            candidate, info = controller.step(controller_state, refs)
            qp_count += 1
            status = info["status"]
            status_counts[status] = status_counts.get(status, 0) + 1
            if candidate is None:
                qp_failure_step = int(step)
                qp_failure_status = str(status)
                break
            control = candidate
            slacks.append(info["minimum_predicted_slack"])
            duals.append(info["maximum_constraint_dual"])
            residuals.append(info["maximum_constraint_residual"])
            solve_times.append(info["solve_time_seconds"])
        task_history.append(task)
        control_history.append(control.copy())
        previous_latent = controller_state[:ell].copy()
        state = (
            model["A"] @ state
            + model["B"] @ (control.reshape(-1, 1) - model["u_mean"])
        )
        state[:ell, 0] += residual_innovation[step]
        state[:ell, 0] += wind_latent_increment[step]

    task_history = np.asarray(task_history, dtype=float).reshape(-1, 3)
    control_history = np.asarray(control_history, dtype=float).reshape(
        -1, model["B"].shape[1],
    )
    completed_steps = int(len(task_history))
    axis_violations = np.abs(task_history) > float(hard_bound)
    step_violations = (
        np.any(axis_violations, axis=1)
        if completed_steps else np.zeros(0, dtype=bool)
    )
    hard_margin = float(hard_bound) - np.abs(task_history)
    saturation = (
        np.any(np.abs(control_history) >= 6.0 - 1e-4, axis=1)
        if completed_steps else np.zeros(0, dtype=bool)
    )
    slack_array = np.asarray(slacks, dtype=float)
    dual_array = np.asarray(duals, dtype=float)
    active_qp = (
        (slack_array <= 5e-4) | (dual_array > 1e-6)
        if slack_array.size else np.zeros(0, dtype=bool)
    )
    # Expand per-decision activity to plant-step length for plotting/metrics.
    active_steps = np.zeros(completed_steps, dtype=bool)
    if completed_steps and active_qp.size:
        for qi, is_active in enumerate(active_qp):
            start = qi * d
            if start >= completed_steps:
                break
            stop = min(start + d, completed_steps)
            if is_active:
                active_steps[start:stop] = True
    empty_axis_float = np.full(3, np.nan)
    empty_axis_int = np.zeros(3, dtype=int)
    return {
        "S": task_history,
        "U": control_history,
        "completed_steps": completed_steps,
        "qp_count": int(qp_count),
        "qp_failure_count": int(qp_failure_step is not None),
        "qp_failure_step": qp_failure_step,
        "qp_failure_status": qp_failure_status,
        "fallback_count": 0,
        "qp_status_counts": status_counts,
        "active_qp_steps": int(np.sum(active_qp)),
        "positive_dual_qp_steps": int(np.sum(dual_array > 1e-6)),
        "active_history": active_steps.astype(np.uint8),
        "minimum_predicted_slack_history": slack_array,
        "maximum_constraint_dual_history": dual_array,
        "maximum_qp_constraint_residual_history": np.asarray(
            residuals, dtype=float,
        ),
        "qp_solve_time_seconds": np.asarray(solve_times, dtype=float),
        "minimum_predicted_slack": (
            float(np.min(slack_array)) if slack_array.size else float("nan")
        ),
        "maximum_constraint_dual": (
            float(np.max(dual_array)) if dual_array.size else float("nan")
        ),
        "mean_qp_solve_time_seconds": (
            float(np.mean(solve_times)) if solve_times else float("nan")
        ),
        "maximum_qp_constraint_residual": (
            float(np.max(residuals)) if residuals else float("nan")
        ),
        "violation_steps": int(np.sum(step_violations)),
        "violation_rate": (
            float(np.mean(step_violations)) if completed_steps else float("nan")
        ),
        "violation_steps_per_axis": (
            np.sum(axis_violations, axis=0).astype(int)
            if completed_steps else empty_axis_int
        ),
        "violation_rate_per_axis": (
            np.mean(axis_violations, axis=0)
            if completed_steps else empty_axis_float.copy()
        ),
        "control_interval_steps": d,
        "horizon_seconds": float(horizon * d * SAMPLE_TIME_SECONDS),
        "minimum_hard_margin": (
            float(np.min(hard_margin)) if completed_steps else float("nan")
        ),
        "minimum_hard_margin_per_axis": (
            np.min(hard_margin, axis=0)
            if completed_steps else empty_axis_float.copy()
        ),
        "rmse": (
            np.sqrt(np.mean(
                (task_history - reference[:completed_steps]) ** 2, axis=0,
            ))
            if completed_steps else empty_axis_float.copy()
        ),
        "input_saturation_steps": int(np.sum(saturation)),
        "disturbance_sha256": _disturbance_sha256(
            residual_innovation, wind_latent_increment,
        ),
        "state_resets": 0,
        "controller_state_source": (
            "true plant augmented state (model-in-the-loop)"
            if use_plant_state
            else "R.T output reconstruction plus previous output"
        ),
        "stage_tightening": controller.stage_tightening,
        "stage_mean_bound": controller.stage_mean_bound,
    }


def run_controller_pair(model, E, hard_bound, reference, residual_innovation,
                        wind_latent_increment, chance_covariance,
                        horizon=HORIZON_STEPS, initial_task=None):
    """Run SMPC and deterministic MPC on identical plant noise/wind.

    ``chance_covariance`` is the process noise used only inside the chance-
    constrained SMPC prediction (identified Σ_eps_aug). Physical wind still
    enters the plant through ``wind_latent_increment`` for both controllers.
    """
    smpc = simulate_controller(
        model, E, hard_bound, reference, residual_innovation,
        wind_latent_increment, chance_covariance, ALPHA_JOINT,
        horizon=horizon, initial_task=initial_task,
    )
    deterministic_mpc = simulate_controller(
        model, E, hard_bound, reference, residual_innovation,
        wind_latent_increment, chance_covariance, None,
        horizon=horizon, initial_task=initial_task,
    )
    if smpc["disturbance_sha256"] != deterministic_mpc["disturbance_sha256"]:
        raise RuntimeError("Controller comparison did not use identical disturbances")
    return smpc, deterministic_mpc


def _metric_payload(result, prefix):
    failure_step = result.get("qp_failure_step")
    return {
        f"{prefix}_completed_steps": np.asarray([
            result.get("completed_steps", len(result["S"]))
        ]),
        f"{prefix}_qp_failure_count": np.asarray([
            result.get("qp_failure_count", 0)
        ]),
        f"{prefix}_qp_failure_step": np.asarray([
            -1 if failure_step is None else failure_step
        ]),
        f"{prefix}_qp_failure_status": str(
            result.get("qp_failure_status", "")
        ),
        f"{prefix}_hard_violation_steps": np.asarray([result["violation_steps"]]),
        f"{prefix}_hard_violation_rate": np.asarray([result["violation_rate"]]),
        f"{prefix}_hard_violation_steps_per_axis": np.asarray(
            result.get("violation_steps_per_axis", np.zeros(3)), dtype=float,
        ),
        f"{prefix}_minimum_hard_margin_standardized": np.asarray(
            result["minimum_hard_margin_per_axis"], dtype=float,
        ),
        f"{prefix}_fallback_count": np.asarray([result["fallback_count"]]),
        f"{prefix}_maximum_qp_constraint_residual": np.asarray([
            result.get("maximum_qp_constraint_residual", 0.0)
        ]),
        f"{prefix}_minimum_hard_margin": np.asarray([
            result.get("minimum_hard_margin", np.min(
                result["minimum_hard_margin_per_axis"]
            ))
        ]),
        f"{prefix}_active_qp_steps": np.asarray([result["active_qp_steps"]]),
        f"{prefix}_positive_dual_qp_steps": np.asarray([
            result["positive_dual_qp_steps"]
        ]),
        f"{prefix}_input_saturation_steps": np.asarray([
            result["input_saturation_steps"]
        ]),
        f"{prefix}_rmse_standardized": np.asarray(result["rmse"], dtype=float),
    }


def build_mat_payload(reference, smpc, deterministic_mpc, residual_innovation,
                      wind_velocity_mps, wind_latent_increment, covariance,
                      noise, scales, hard_bound, stage_tightening,
                      initialization=None):
    reference = np.asarray(reference, dtype=float)
    wind_velocity_mps = np.asarray(wind_velocity_mps, dtype=float)
    position_scale = np.asarray(scales["y_scale"][:3], dtype=float)
    position_offset = np.asarray(scales["y_offset"][:3], dtype=float)

    def physical_position(standardized):
        return np.asarray(standardized, dtype=float) * position_scale + position_offset

    dt_seconds = float(covariance.get("dt_seconds", SAMPLE_TIME_SECONDS))
    stage_tightening = np.asarray(stage_tightening, dtype=float)
    payload = {
        "reference_standardized": reference,
        "reference_position_m": physical_position(reference),
        "smpc_trajectory_standardized": np.asarray(smpc["S"], dtype=float),
        "smpc_position_m": physical_position(smpc["S"]),
        "deterministic_mpc_trajectory_standardized": np.asarray(
            deterministic_mpc["S"], dtype=float,
        ),
        "deterministic_mpc_position_m": physical_position(deterministic_mpc["S"]),
        "smpc_control": np.asarray(smpc["U"], dtype=float),
        "deterministic_mpc_control": np.asarray(deterministic_mpc["U"], dtype=float),
        "residual_innovation_latent": np.asarray(residual_innovation, dtype=float),
        "innovation_latent": np.asarray(residual_innovation, dtype=float),
        "wind_velocity_mps": wind_velocity_mps,
        "wind_position_delta_m": wind_velocity_mps * dt_seconds,
        "wind_position_delta_xyz_m": wind_velocity_mps * dt_seconds,
        "wind_position_delta_standardized": (
            wind_velocity_mps * dt_seconds / position_scale
        ),
        "wind_latent_increment": np.asarray(wind_latent_increment, dtype=float),
        "Sigma_eps": np.asarray(noise["Sigma_eps"], dtype=float),
        "Sigma_obs_proxy": np.asarray(noise["Sigma_obs_proxy"], dtype=float),
        "Sigma_wind_aug": np.asarray(covariance["Sigma_wind_aug"], dtype=float),
        "Sigma_total_aug": np.asarray(covariance["Sigma_total_aug"], dtype=float),
        "G_w": np.asarray(covariance["G_w"], dtype=float),
        "position_offset_m": position_offset,
        "position_scale_m": position_scale,
        "hard_bound_standardized": np.asarray([hard_bound], dtype=float),
        "hard_position_lower_m": position_offset - float(hard_bound) * position_scale,
        "hard_position_upper_m": position_offset + float(hard_bound) * position_scale,
        "stage_tightening_standardized": stage_tightening,
        "stage_mean_bound_standardized": float(hard_bound) - stage_tightening,
        "dt_seconds": np.asarray([dt_seconds], dtype=float),
        "dt_used_seconds": np.asarray([dt_seconds], dtype=float),
        "sample_rate_hz": np.asarray([1.0 / dt_seconds], dtype=float),
        "horizon_steps": np.asarray([stage_tightening.shape[0]], dtype=float),
        "horizon_seconds": np.asarray([
            stage_tightening.shape[0] * CONTROL_INTERVAL_STEPS * dt_seconds
        ], dtype=float),
        "control_interval_steps": np.asarray(
            [CONTROL_INTERVAL_STEPS], dtype=float,
        ),
        "sigma_wind_mps": np.asarray([
            covariance.get("sigma_wind_mps", SIGMA_WIND_MPS)
        ], dtype=float),
        "uses_true_Sigma_n": np.asarray([noise["uses_true_Sigma_n"]], dtype=float),
        "training_segments": np.asarray([noise["training_segments"]], dtype=float),
        "transition_count": np.asarray([noise["transition_count"]], dtype=float),
        "shared_noise_exact": np.asarray([1], dtype=float),
        "disturbance_sha256": str(smpc.get("disturbance_sha256", "")),
        "position_unit": "m",
        "unit_note": (
            "wind_velocity_mps is m/s; wind_position_delta_m = "
            "wind_velocity_mps * 0.01 s; physical positions are meters"
        ),
        **_metric_payload(smpc, "smpc"),
        **_metric_payload(deterministic_mpc, "deterministic_mpc"),
    }
    if initialization is not None:
        payload.update({
            "initial_state": np.asarray(initialization["state"], dtype=float),
            "initial_control": np.asarray(initialization["control"], dtype=float),
            "initial_dynamic_residual": np.asarray([
                initialization["dynamic_residual"]
            ]),
            "initial_task_residual": np.asarray([
                initialization["task_residual"]
            ]),
            "initial_maximum_absolute_input": np.asarray([
                initialization["maximum_absolute_input"]
            ]),
        })
    return payload


def _json_ready(value):
    if isinstance(value, dict):
        return {key: _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_ready(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    return value


def _controller_summary(result):
    trajectory = np.asarray(result["S"], dtype=float)
    peak = (
        np.max(np.abs(trajectory), axis=0)
        if trajectory.size else np.full(3, np.nan)
    )
    return _json_ready({
        "completed_steps": result["completed_steps"],
        "qp_count": result["qp_count"],
        "qp_failure_count": result["qp_failure_count"],
        "qp_failure_step": result["qp_failure_step"],
        "qp_failure_status": result["qp_failure_status"],
        "fallback_count": result["fallback_count"],
        "qp_status_counts": result["qp_status_counts"],
        "maximum_qp_constraint_residual": result[
            "maximum_qp_constraint_residual"
        ],
        "mean_qp_solve_time_seconds": result["mean_qp_solve_time_seconds"],
        "active_qp_steps": result["active_qp_steps"],
        "positive_dual_qp_steps": result["positive_dual_qp_steps"],
        "minimum_predicted_slack": result["minimum_predicted_slack"],
        "maximum_constraint_dual": result["maximum_constraint_dual"],
        "hard_violation_steps": result["violation_steps"],
        "hard_violation_rate": result["violation_rate"],
        "hard_violation_steps_per_axis": result["violation_steps_per_axis"],
        "hard_violation_rate_per_axis": result["violation_rate_per_axis"],
        "minimum_hard_margin_standardized": result["minimum_hard_margin"],
        "minimum_hard_margin_per_axis_standardized": result[
            "minimum_hard_margin_per_axis"
        ],
        "rmse_standardized": result["rmse"],
        "peak_absolute_position_standardized": peak,
        "input_saturation_steps": result["input_saturation_steps"],
        "state_resets": result["state_resets"],
        "controller_state_source": result["controller_state_source"],
        "disturbance_sha256": result["disturbance_sha256"],
    })


def _draw_bound_box(axis, bound, color, linestyle, label):
    values = (-float(bound), float(bound))
    first = True
    for x in values:
        for y in values:
            axis.plot(
                [x, x], [y, y], [values[0], values[1]],
                color=color, linestyle=linestyle, linewidth=0.8,
                label=label if first else None,
            )
            first = False
    for x in values:
        for z in values:
            axis.plot(
                [x, x], [values[0], values[1]], [z, z],
                color=color, linestyle=linestyle, linewidth=0.8,
            )
    for y in values:
        for z in values:
            axis.plot(
                [values[0], values[1]], [y, y], [z, z],
                color=color, linestyle=linestyle, linewidth=0.8,
            )


def plot_controller_diagnostics(path, controller_name, result, reference,
                                hard_bound, smpc_stage_mean_bound,
                                sigma_wind_mps):
    reference = np.asarray(reference, dtype=float)
    trajectory = np.asarray(result["S"], dtype=float)
    controls = np.asarray(result["U"], dtype=float)
    completed = result["completed_steps"]
    dt = SAMPLE_TIME_SECONDS
    time_full = np.arange(len(reference)) * dt
    time = np.arange(completed) * dt
    colors = ("#006D77", "#D1495B", "#6A4C93")
    labels = ("x", "y", "z")
    tightened = np.asarray(smpc_stage_mean_bound, dtype=float)[-1]

    figure = plt.figure(figsize=(16, 13), constrained_layout=True)
    grid = figure.add_gridspec(4, 3, height_ratios=(1.0, 1.15, 1.0, 0.85))
    position_axes = [figure.add_subplot(grid[0, axis]) for axis in range(3)]
    trajectory_axis = figure.add_subplot(grid[1:3, 0], projection="3d")
    error_axis = figure.add_subplot(grid[1, 1:])
    input_axis = figure.add_subplot(grid[2, 1:])
    margin_axis = figure.add_subplot(grid[3, :])

    for axis_index, axis in enumerate(position_axes):
        axis.plot(
            time_full, reference[:, axis_index], color="#555555",
            linestyle="--", linewidth=1.2, label="reference",
        )
        if completed:
            axis.plot(
                time, trajectory[:, axis_index], color=colors[axis_index],
                linewidth=1.4, label=controller_name,
            )
        axis.axhline(hard_bound, color="#202020", linewidth=1.0,
                    label="hard bound" if axis_index == 0 else None)
        axis.axhline(-hard_bound, color="#202020", linewidth=1.0)
        axis.axhline(tightened[axis_index], color="#E09F3E", linewidth=1.0,
                    linestyle=":", label="SMPC terminal mean bound")
        axis.axhline(-tightened[axis_index], color="#E09F3E",
                    linewidth=1.0, linestyle=":")
        axis.set_title(f"{labels[axis_index]} position")
        axis.set_xlabel("time [s]")
        axis.set_ylabel("standardized position")
        axis.grid(alpha=0.2)
        if axis_index == 0:
            axis.legend(loc="best", fontsize=8)

    trajectory_axis.plot(
        reference[:, 0], reference[:, 1], reference[:, 2],
        color="#777777", linestyle="--", linewidth=1.0, label="reference",
    )
    if completed:
        trajectory_axis.plot(
            trajectory[:, 0], trajectory[:, 1], trajectory[:, 2],
            color="#006D77", linewidth=1.4, label=controller_name,
        )
    _draw_bound_box(trajectory_axis, hard_bound, "#333333", "-", "hard box")
    _draw_bound_box(
        trajectory_axis, float(np.min(tightened)), "#E09F3E", ":",
        "SMPC terminal mean box",
    )
    trajectory_axis.set_xlabel("x standardized")
    trajectory_axis.set_ylabel("y standardized")
    trajectory_axis.set_zlabel("z standardized")
    trajectory_axis.set_title("3D trajectory and constraint boxes")
    trajectory_axis.legend(loc="upper left", fontsize=8)
    trajectory_axis.view_init(elev=24, azim=-55)

    if completed:
        error = trajectory - reference[:completed]
        for axis_index in range(3):
            error_axis.plot(
                time, error[:, axis_index], color=colors[axis_index],
                linewidth=1.0, label=f"{labels[axis_index]} error",
            )
        for channel in range(controls.shape[1]):
            input_axis.plot(
                time, controls[:, channel], linewidth=1.0,
                label=f"u{channel + 1}",
            )
        margin = hard_bound - np.max(np.abs(trajectory), axis=1)
        margin_axis.plot(time, margin, color="#006D77", linewidth=1.2,
                         label="minimum hard margin")
        active = np.asarray(result["active_history"], dtype=bool)
        if np.any(active):
            margin_axis.scatter(
                time[active], margin[active], color="#D1495B", s=16,
                zorder=3, label="active predicted bound",
            )
    error_axis.axhline(0.0, color="#333333", linewidth=0.8)
    error_axis.set_title("Tracking errors")
    error_axis.set_xlabel("time [s]")
    error_axis.set_ylabel("actual - reference")
    error_axis.grid(alpha=0.2)
    error_axis.legend(loc="best", ncol=3, fontsize=8)
    input_axis.axhline(6.0, color="#333333", linewidth=0.8)
    input_axis.axhline(-6.0, color="#333333", linewidth=0.8)
    input_axis.set_title("Normalized control inputs")
    input_axis.set_xlabel("time [s]")
    input_axis.set_ylabel("input")
    input_axis.grid(alpha=0.2)
    input_axis.legend(loc="best", ncol=4, fontsize=8)
    margin_axis.axhline(0.0, color="#202020", linewidth=1.0,
                        label="hard-bound contact")
    margin_axis.set_title("Realized hard margin and predicted-bound activation")
    margin_axis.set_xlabel("time [s]")
    margin_axis.set_ylabel("standardized margin")
    margin_axis.grid(alpha=0.2)
    margin_axis.legend(loc="best", ncol=3, fontsize=8)

    failure = result.get("qp_failure_step")
    if failure is not None:
        failure_time = failure * dt
        for axis in position_axes + [error_axis, input_axis, margin_axis]:
            axis.axvline(failure_time, color="#D1495B", linestyle="--",
                         linewidth=1.1)
        margin_axis.text(
            failure_time, 0.98, "QP infeasible: simulation stopped",
            transform=margin_axis.get_xaxis_transform(), ha="right",
            va="top", color="#A4243B", fontsize=9,
        )
    status = (
        "completed" if failure is None
        else f"stopped at step {failure} ({result['qp_failure_status']})"
    )
    figure.suptitle(
        f"copyBG {controller_name}: Pelican model-in-the-loop, "
        f"physical Gaussian wind sigma={sigma_wind_mps:.3f} m/s; {status}"
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)


def summarize_noise_actual_vs_estimated(
        residual_innovation,
        wind_velocity_mps,
        wind_latent_increment,
        Sigma_eps,
        Sigma_wind_aug,
        Sigma_total_aug,
        sigma_wind_mps,
        dt_seconds=SAMPLE_TIME_SECONDS):
    """Compare realized plant noise to the covariances used / estimated offline.

    Layers (must stay distinct):
      - physical wind v ~ N(0, sigma_wind^2 I) in m/s  (prescribed, known)
      - latent innovation eps ~ N(0, Sigma_eps) from train residuals (sampled)
      - plant step disturbance = eps + G_w @ (v * dt mapped) = eps + wind_z
      - SMPC chance covariance uses Sigma_eps only (not wind, not Sigma_n)
      - Sigma_obs_proxy is diagnostic only and is not compared here
    """
    eps = np.asarray(residual_innovation, dtype=float)
    wind_v = np.asarray(wind_velocity_mps, dtype=float)
    wind_z = np.asarray(wind_latent_increment, dtype=float)
    Sigma_eps = np.asarray(Sigma_eps, dtype=float)
    Sigma_wind_aug = np.asarray(Sigma_wind_aug, dtype=float)
    Sigma_total_aug = np.asarray(Sigma_total_aug, dtype=float)
    if eps.ndim != 2 or wind_z.ndim != 2 or wind_v.ndim != 2:
        raise ValueError("noise arrays must be 2-D")
    if eps.shape != wind_z.shape:
        raise ValueError("innovation and wind_latent shapes must match")
    if wind_v.shape[0] != eps.shape[0] or wind_v.shape[1] != 3:
        raise ValueError("wind_velocity must be (T,3)")
    n_steps, ell = eps.shape
    total = eps + wind_z

    def _emp_cov(x):
        return (x.T @ x) / max(x.shape[0], 1)

    emp_eps = _emp_cov(eps)
    emp_wind_z = _emp_cov(wind_z)
    emp_total = _emp_cov(total)
    Sigma_wind = Sigma_wind_aug[:ell, :ell]
    Sigma_total = Sigma_total_aug[:ell, :ell]
    Sigma_eps_ell = Sigma_eps if Sigma_eps.shape[0] == ell else Sigma_eps[:ell, :ell]

    wind_std_actual = wind_v.std(axis=0, ddof=0)
    wind_std_estimated = np.full(3, float(sigma_wind_mps))
    eps_std_actual = np.sqrt(np.maximum(np.diag(emp_eps), 0.0))
    eps_std_estimated = np.sqrt(np.maximum(np.diag(Sigma_eps_ell), 0.0))
    wind_z_std_actual = np.sqrt(np.maximum(np.diag(emp_wind_z), 0.0))
    wind_z_std_estimated = np.sqrt(np.maximum(np.diag(Sigma_wind), 0.0))
    total_std_actual = np.sqrt(np.maximum(np.diag(emp_total), 0.0))
    total_std_est_eps_only = eps_std_estimated  # what chance uses
    total_std_est_total = np.sqrt(np.maximum(np.diag(Sigma_total), 0.0))

    def _fro_rel(a, b):
        nb = float(np.linalg.norm(b, ord="fro"))
        denom = max(nb, 1e-15)
        return float(np.linalg.norm(a - b, ord="fro") / denom)

    return {
        "n_steps": int(n_steps),
        "ell": int(ell),
        "dt_seconds": float(dt_seconds),
        "sigma_wind_mps_prescribed": float(sigma_wind_mps),
        "wind_velocity_std_actual_mps": wind_std_actual,
        "wind_velocity_std_estimated_mps": wind_std_estimated,
        "wind_velocity_std_ratio_actual_over_est": (
            wind_std_actual / np.maximum(wind_std_estimated, 1e-15)
        ),
        "eps_std_actual": eps_std_actual,
        "eps_std_estimated": eps_std_estimated,
        "eps_std_ratio_actual_over_est": (
            eps_std_actual / np.maximum(eps_std_estimated, 1e-15)
        ),
        "wind_latent_std_actual": wind_z_std_actual,
        "wind_latent_std_estimated": wind_z_std_estimated,
        "total_plant_std_actual": total_std_actual,
        "total_plant_std_est_chance_Sigma_eps": total_std_est_eps_only,
        "total_plant_std_est_Sigma_total": total_std_est_total,
        "emp_cov_eps": emp_eps,
        "emp_cov_wind_latent": emp_wind_z,
        "emp_cov_total_plant": emp_total,
        "Sigma_eps": Sigma_eps_ell,
        "Sigma_wind": Sigma_wind,
        "Sigma_total": Sigma_total,
        "fro_rel_emp_eps_vs_Sigma_eps": _fro_rel(emp_eps, Sigma_eps_ell),
        "fro_rel_emp_wind_z_vs_Sigma_wind": _fro_rel(emp_wind_z, Sigma_wind),
        "fro_rel_emp_total_vs_Sigma_eps_chance": _fro_rel(
            emp_total, Sigma_eps_ell,
        ),
        "fro_rel_emp_total_vs_Sigma_total": _fro_rel(emp_total, Sigma_total),
        "chance_uses": "Sigma_eps only (not physical wind, not Sigma_n)",
        "plant_injects": "eps ~ N(0,Sigma_eps) + wind_latent (physical v)",
        "note": (
            "Actual = realized closed-loop sequences. Estimated = offline "
            "Sigma_eps from train innovations + analytic Sigma_wind from "
            "sigma_wind and G_w. SMPC chance tightening uses Sigma_eps only; "
            "Sigma_total is a diagnostic upper envelope including wind."
        ),
    }


def _rolling_rms(x, window):
    """Causal rolling RMS for 1-D series (expanding until window fills)."""
    x = np.asarray(x, dtype=float).ravel()
    n = int(x.shape[0])
    if n == 0:
        return x.copy()
    window = int(max(1, min(int(window), n)))
    c2 = np.concatenate(([0.0], np.cumsum(x * x, dtype=float)))
    idx = np.arange(1, n + 1)
    left = np.maximum(0, idx - window)
    count = (idx - left).astype(float)
    s2 = c2[idx] - c2[left]
    return np.sqrt(np.maximum(s2, 0.0) / np.maximum(count, 1.0))


def _downsample_idx(n, max_points):
    stride = max(1, int(np.ceil(float(n) / float(max(max_points, 1)))))
    return slice(None, None, stride)


def plot_noise_actual_vs_estimated(
        path,
        residual_innovation,
        wind_velocity_mps,
        wind_latent_increment,
        Sigma_eps,
        Sigma_wind_aug,
        Sigma_total_aug,
        sigma_wind_mps,
        dt_seconds=SAMPLE_TIME_SECONDS,
        max_trace_points=2500,
        roll_window_seconds=8.0):
    """Time-series actual-vs-estimated noise figure.

    Read top→bottom:
      1) actual wind waveform inside estimated ±σ band
      2) rolling wind std hugs estimated σ line
      3) dominant latent innovation inside estimated ±σ band
      4) rolling latent std hugs estimated σ lines
      5) rolling match ratio actual/est → 1 over time
    """
    metrics = summarize_noise_actual_vs_estimated(
        residual_innovation, wind_velocity_mps, wind_latent_increment,
        Sigma_eps, Sigma_wind_aug, Sigma_total_aug, sigma_wind_mps,
        dt_seconds=dt_seconds,
    )
    eps = np.asarray(residual_innovation, dtype=float)
    wind_v = np.asarray(wind_velocity_mps, dtype=float)
    n_steps = eps.shape[0]
    ell = int(metrics["ell"])
    sigma = float(sigma_wind_mps)
    dt = float(dt_seconds)
    time = np.arange(n_steps, dtype=float) * dt
    sl = _downsample_idx(n_steps, max_trace_points)
    win = max(20, int(round(float(roll_window_seconds) / max(dt, 1e-9))))

    wind_ratio = np.asarray(metrics["wind_velocity_std_ratio_actual_over_est"])
    eps_ratio = np.asarray(metrics["eps_std_ratio_actual_over_est"])
    wind_err_pct = float(np.max(100.0 * np.abs(wind_ratio - 1.0)))
    eps_err_pct = float(np.max(100.0 * np.abs(eps_ratio - 1.0)))
    fro_eps_pct = 100.0 * float(metrics["fro_rel_emp_eps_vs_Sigma_eps"])
    fro_tot_pct = 100.0 * float(metrics["fro_rel_emp_total_vs_Sigma_total"])
    ok = (wind_err_pct < 5.0) and (fro_eps_pct < 5.0) and (eps_err_pct < 10.0)
    verdict = "MATCH" if ok else "CHECK"
    color_v = "#1B7F4E" if ok else "#A4243B"

    eps_std_est = np.asarray(metrics["eps_std_estimated"], dtype=float)
    # Show the two largest latent modes (readable scale); keep others in ratio row.
    mode_order = np.argsort(eps_std_est)[::-1]
    mode_a = int(mode_order[0])
    mode_b = int(mode_order[1]) if ell > 1 else mode_a
    sa = float(max(eps_std_est[mode_a], 1e-15))
    sb = float(max(eps_std_est[mode_b], 1e-15))

    # Rolling statistics (full rate), then downsample for draw.
    roll_wind = np.mean(
        np.column_stack([_rolling_rms(wind_v[:, i], win) for i in range(3)]),
        axis=1,
    )
    roll_eps_a = _rolling_rms(eps[:, mode_a], win)
    roll_eps_b = _rolling_rms(eps[:, mode_b], win)
    roll_wind_ratio = roll_wind / max(sigma, 1e-15)
    roll_eps_a_ratio = roll_eps_a / sa
    roll_eps_b_ratio = roll_eps_b / sb

    figure = plt.figure(figsize=(15, 12), constrained_layout=True)
    axes = figure.subplots(5, 1, sharex=True)

    # --- 1) wind time series + estimated ±σ band ---
    ax = axes[0]
    ax.fill_between(
        time[sl], -sigma, sigma, color="#E09F3E", alpha=0.22,
        label=r"estimated $\pm\sigma$ band",
    )
    ax.plot(time[sl], wind_v[sl, 0], color="#006D77", linewidth=0.55,
            alpha=0.9, label="actual wind x")
    ax.plot(time[sl], wind_v[sl, 1], color="#D1495B", linewidth=0.45,
            alpha=0.55, label="actual wind y")
    ax.plot(time[sl], wind_v[sl, 2], color="#6A4C93", linewidth=0.45,
            alpha=0.55, label="actual wind z")
    ax.axhline(sigma, color="#8B5E00", linestyle="--", linewidth=1.1)
    ax.axhline(-sigma, color="#8B5E00", linestyle="--", linewidth=1.1)
    ax.axhline(0.0, color="#888888", linewidth=0.5)
    ax.set_ylabel("wind [m/s]")
    ax.set_title(
        "1) Actual wind over time  vs  long-run estimated ±σ band  "
        f"(σ={sigma:.3f} m/s, gusty).  High-gust stretches exceed the band.",
        loc="left", fontsize=11, fontweight="bold",
    )
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", ncol=4, fontsize=8)

    # --- 2) rolling wind std vs estimated σ (gusty: actual breathes) ---
    ax = axes[1]
    ax.plot(time[sl], roll_wind[sl], color="#006D77", linewidth=1.8,
            label=f"actual rolling std (window={roll_window_seconds:.0f}s)")
    ax.axhline(sigma, color="#E09F3E", linestyle="-", linewidth=2.0,
               label=r"long-run estimated $\sigma$ (stationary)")
    # Prescribed intensity shape × long-run σ (guide only; plant uses realized wind).
    try:
        env = wind_intensity_envelope(n_steps, dt_seconds=dt)
        env_scale = float(np.sqrt(np.mean(env * env)))
        guide = sigma * (env / max(env_scale, 1e-15))
        ax.plot(time[sl], guide[sl], color="#9B2226", linestyle="--",
                linewidth=1.3, alpha=0.85,
                label="prescribed gust envelope × σ")
    except Exception:
        guide = None
    ax.fill_between(
        time[sl], 0.95 * sigma, 1.05 * sigma,
        color="#1B7F4E", alpha=0.10, label="±5% of long-run σ",
    )
    ax.set_ylabel("std [m/s]")
    ax.set_title(
        "2) Rolling actual wind std  breathes over time  "
        f"(gusty plant wind; long-run actual/est max err {wind_err_pct:.2f}%)",
        loc="left", fontsize=11, fontweight="bold",
    )
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", ncol=2, fontsize=8)
    y_candidates = [float(np.nanmax(roll_wind)), sigma]
    if guide is not None:
        y_candidates.append(float(np.nanmax(guide)))
    ymax = max(y_candidates) * 1.35
    ax.set_ylim(0.0, max(ymax, 1e-6))

    # --- 3) dominant latent innovation + ±σ band ---
    ax = axes[2]
    ax.fill_between(
        time[sl], -sa, sa, color="#E09F3E", alpha=0.20,
        label=fr"est $\pm\sigma$  z{mode_a+1}",
    )
    ax.plot(time[sl], eps[sl, mode_a], color="#006D77", linewidth=0.55,
            alpha=0.9, label=f"actual eps z{mode_a+1}")
    if mode_b != mode_a:
        ax.plot(time[sl], eps[sl, mode_b], color="#D1495B", linewidth=0.45,
                alpha=0.55, label=f"actual eps z{mode_b+1}")
        ax.axhline(sb, color="#D1495B", linestyle=":", linewidth=1.0, alpha=0.8)
        ax.axhline(-sb, color="#D1495B", linestyle=":", linewidth=1.0, alpha=0.8)
    ax.axhline(sa, color="#8B5E00", linestyle="--", linewidth=1.1)
    ax.axhline(-sa, color="#8B5E00", linestyle="--", linewidth=1.1)
    ax.axhline(0.0, color="#888888", linewidth=0.5)
    ax.set_ylabel("latent eps")
    ax.set_title(
        f"3) Actual latent innovation over time  vs  estimated ±σ "
        f"(largest modes z{mode_a+1}, z{mode_b+1}).  "
        "Chance noise = Σε only.",
        loc="left", fontsize=11, fontweight="bold",
    )
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", ncol=4, fontsize=8)

    # --- 4) rolling latent std vs estimated ---
    ax = axes[3]
    ax.plot(time[sl], roll_eps_a[sl], color="#006D77", linewidth=1.6,
            label=f"actual roll std z{mode_a+1}")
    ax.axhline(sa, color="#E09F3E", linewidth=2.0,
               label=fr"est $\sigma$ z{mode_a+1}={sa:.4g}")
    if mode_b != mode_a:
        ax.plot(time[sl], roll_eps_b[sl], color="#D1495B", linewidth=1.3,
                alpha=0.85, label=f"actual roll std z{mode_b+1}")
        ax.axhline(sb, color="#D1495B", linestyle="--", linewidth=1.5,
                   label=fr"est $\sigma$ z{mode_b+1}={sb:.4g}")
    ax.set_ylabel("std")
    ax.set_title(
        "4) Rolling actual latent std  hugs  estimated σ lines  "
        f"(Σε Fro err {fro_eps_pct:.2f}%)",
        loc="left", fontsize=11, fontweight="bold",
    )
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", ncol=2, fontsize=8)
    ymax = max(float(np.nanmax(roll_eps_a)), float(np.nanmax(roll_eps_b)), sa, sb)
    ax.set_ylim(0.0, max(ymax * 1.35, 1e-9))

    # --- 5) rolling match ratio → 1 ---
    ax = axes[4]
    ax.plot(time[sl], roll_wind_ratio[sl], color="#006D77", linewidth=1.5,
            label="wind actual/est")
    ax.plot(time[sl], roll_eps_a_ratio[sl], color="#E09F3E", linewidth=1.5,
            label=f"eps z{mode_a+1} actual/est")
    if mode_b != mode_a:
        ax.plot(time[sl], roll_eps_b_ratio[sl], color="#D1495B", linewidth=1.2,
                alpha=0.85, label=f"eps z{mode_b+1} actual/est")
    ax.axhline(1.0, color="#111111", linewidth=2.0, label="perfect match (=1)")
    ax.axhspan(0.95, 1.05, color="#1B7F4E", alpha=0.12, label="±5%")
    ax.set_ylabel("actual / estimated")
    ax.set_xlabel("time [s]")
    ax.set_ylim(0.5, 1.5)
    ax.set_title(
        "5) Match ratio over time  (→ 1 means estimate tracks actual).  "
        f"VERDICT {verdict}: wind max err {wind_err_pct:.2f}%, "
        f"eps max err {eps_err_pct:.2f}%, Σε Fro {fro_eps_pct:.2f}%",
        loc="left", fontsize=11, fontweight="bold", color=color_v,
    )
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", ncol=3, fontsize=8)

    figure.suptitle(
        "copyBG actual noise vs estimated noise  |  time series  |  "
        f"wind long-run σ={sigma:.3f} m/s (gusty)  |  "
        f"roll window={roll_window_seconds:.0f}s  |  "
        "chance uses Σε only (uses_true_Σn=0)",
        fontsize=12, fontweight="bold",
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=170)
    plt.close(figure)

    summary = {
        key: (value.tolist() if isinstance(value, np.ndarray) else value)
        for key, value in metrics.items()
        if key not in (
            "emp_cov_eps", "emp_cov_wind_latent", "emp_cov_total_plant",
            "Sigma_eps", "Sigma_wind", "Sigma_total",
        )
    }
    summary.update({
        "figure": str(path.name),
        "verdict": verdict,
        "wind_max_abs_error_pct": wind_err_pct,
        "eps_max_abs_std_error_pct": eps_err_pct,
        "eps_fro_error_pct": fro_eps_pct,
        "total_fro_error_pct": fro_tot_pct,
        "roll_window_seconds": float(roll_window_seconds),
        "shown_latent_modes": [mode_a + 1, mode_b + 1],
        "read_guide": (
            "Top panels: actual waveform inside estimated ±σ band. "
            "Middle: rolling actual std should hug the estimated σ line. "
            "Bottom: actual/est ratio should stay near 1 over time."
        ),
    })
    return summary


def _smooth_unit_tangents(points, window=21):
    """Unit tangents along a 3D polyline; dwells fall back to last nonzero."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or points.shape[0] < 2:
        raise ValueError("points must be (N,3) with N>=2")
    deltas = np.diff(points, axis=0)
    # Forward difference at each sample; last copies previous.
    vel = np.vstack((deltas, deltas[-1:]))
    if window > 1:
        kernel = np.ones(int(window), dtype=float)
        kernel /= kernel.sum()
        pad = len(kernel) // 2
        vel_pad = np.pad(vel, ((pad, pad), (0, 0)), mode="edge")
        vel = np.stack(
            [np.convolve(vel_pad[:, axis], kernel, mode="valid") for axis in range(3)],
            axis=1,
        )
    norms = np.linalg.norm(vel, axis=1, keepdims=True)
    tangents = np.zeros_like(vel)
    moving = norms.ravel() > 1e-12
    tangents[moving] = vel[moving] / norms[moving]
    # Fill stationary samples with previous moving tangent (or +x).
    last = np.array([1.0, 0.0, 0.0])
    for index in range(len(tangents)):
        if moving[index]:
            last = tangents[index]
        else:
            tangents[index] = last
    return tangents


def _smooth_camera_angles(azim_deg, elev_deg, window=11):
    """Unwrap azimuth and lightly smooth (azim, elev) to cut GIF whip-pan."""
    azim = np.unwrap(np.deg2rad(np.asarray(azim_deg, dtype=float)))
    elev = np.asarray(elev_deg, dtype=float)
    window = max(3, int(window) | 1)
    kernel = np.ones(window, dtype=float) / float(window)
    pad = window // 2
    azim = np.convolve(np.pad(azim, (pad, pad), mode="edge"), kernel, mode="valid")
    elev = np.convolve(np.pad(elev, (pad, pad), mode="edge"), kernel, mode="valid")
    return np.rad2deg(azim), elev


def render_forward_view_gif(
        path, trajectory, reference, hard_bound,
        stage_mean_bound=None,
        dt=SAMPLE_TIME_SECONDS,
        max_frames=FORWARD_GIF_MAX_FRAMES,
        fps=FORWARD_GIF_FPS,
        title="copyBG SMPC — forward chase view"):
    """GIF from behind the vehicle looking along its forward/heading direction.

    Camera is a chase/FPV hybrid: eye sits slightly behind and above the
    current position, looking along the path tangent (drone nose / progress
    direction). Full reference ghost + flown trail + hard box stay in view.
    """
    from matplotlib.animation import FuncAnimation, PillowWriter

    trajectory = np.asarray(trajectory, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if trajectory.ndim != 2 or trajectory.shape[1] != 3:
        raise ValueError("trajectory must be (N,3)")
    completed = int(trajectory.shape[0])
    if completed < 2:
        raise ValueError("need at least 2 trajectory samples for a GIF")
    hard_bound = float(hard_bound)
    max_frames = max(2, int(max_frames))
    fps = max(1, int(fps))

    # Subsample evenly across the completed run.
    if completed <= max_frames:
        frame_indices = np.arange(completed, dtype=int)
    else:
        frame_indices = np.unique(
            np.linspace(0, completed - 1, max_frames, dtype=int)
        )
    tangents = _smooth_unit_tangents(trajectory, window=31)
    span = max(4.5, 1.35 * hard_bound)
    mean_box = None
    if stage_mean_bound is not None:
        mean_box = float(np.min(np.asarray(stage_mean_bound, dtype=float)[-1]))

    # Precompute chase angles on the subsampled timeline, then smooth so
    # 90° box corners do not whip the camera every few frames.
    raw_azim = []
    raw_elev = []
    for index in frame_indices:
        eye_dir = -tangents[int(index)]
        elev_i = float(np.degrees(np.arcsin(np.clip(eye_dir[2], -0.55, 0.55))))
        # Prefer a stable high-ish elevation (less pitch rock = less dizziness).
        elev_i = float(np.clip(18.0 + 0.35 * elev_i, 14.0, 32.0))
        azim_i = float(np.degrees(np.arctan2(eye_dir[1], eye_dir[0])))
        raw_azim.append(azim_i)
        raw_elev.append(elev_i)
    smooth_azim, smooth_elev = _smooth_camera_angles(raw_azim, raw_elev, window=11)

    figure = plt.figure(figsize=(7.2, 6.4), facecolor="#0b1020")
    axis = figure.add_subplot(111, projection="3d", facecolor="#0b1020")
    figure.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.92)

    ref_line, = axis.plot(
        reference[:, 0], reference[:, 1], reference[:, 2],
        color="#6b7280", linestyle="--", linewidth=0.9, alpha=0.55,
        label="reference",
    )
    trail_line, = axis.plot([], [], [], color="#2dd4bf", linewidth=2.0,
                            label="flown")
    future_line, = axis.plot([], [], [], color="#38bdf8", linewidth=1.2,
                             alpha=0.85, label="ahead")
    drone_dot = axis.scatter([], [], [], color="#f8fafc", s=48, depthshade=False)
    heading_quiver = None
    hud = axis.text2D(
        0.02, 0.96, "", transform=axis.transAxes, color="#e5e7eb",
        fontsize=9, family="monospace", va="top",
    )
    _draw_bound_box(axis, hard_bound, "#f87171", "-", "hard box")
    if mean_box is not None and mean_box > 0.0:
        _draw_bound_box(axis, mean_box, "#fbbf24", ":", "SMPC mean box")
    axis.set_xlabel("x", color="#cbd5e1")
    axis.set_ylabel("y", color="#cbd5e1")
    axis.set_zlabel("z", color="#cbd5e1")
    axis.tick_params(colors="#94a3b8", labelsize=7)
    axis.xaxis.pane.fill = False
    axis.yaxis.pane.fill = False
    axis.zaxis.pane.fill = False
    axis.xaxis.pane.set_edgecolor("#1f2937")
    axis.yaxis.pane.set_edgecolor("#1f2937")
    axis.zaxis.pane.set_edgecolor("#1f2937")
    try:
        axis.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    legend = axis.legend(
        loc="upper right", fontsize=7, framealpha=0.35,
        facecolor="#111827", edgecolor="#334155", labelcolor="#e5e7eb",
    )
    title_artist = figure.suptitle(title, color="#f8fafc", fontsize=12)

    def _set_limits(center):
        axis.set_xlim(center[0] - span, center[0] + span)
        axis.set_ylim(center[1] - span, center[1] + span)
        axis.set_zlim(center[2] - span, center[2] + span)

    def _update(frame_no):
        nonlocal heading_quiver
        index = int(frame_indices[frame_no])
        position = trajectory[index]
        tangent = tangents[index]
        elev = float(smooth_elev[frame_no])
        azim = float(smooth_azim[frame_no])
        look = position + 0.45 * span * tangent
        _set_limits(look)
        axis.view_init(elev=elev, azim=azim)

        trail_line.set_data(trajectory[: index + 1, 0], trajectory[: index + 1, 1])
        trail_line.set_3d_properties(trajectory[: index + 1, 2])
        ahead = trajectory[index: min(completed, index + max(40, completed // 30))]
        if ahead.shape[0] >= 2:
            future_line.set_data(ahead[:, 0], ahead[:, 1])
            future_line.set_3d_properties(ahead[:, 2])
        else:
            future_line.set_data([], [])
            future_line.set_3d_properties([])

        drone_dot._offsets3d = (
            np.array([position[0]]),
            np.array([position[1]]),
            np.array([position[2]]),
        )
        if heading_quiver is not None:
            heading_quiver.remove()
        arrow_len = 0.22 * span
        heading_quiver = axis.quiver(
            position[0], position[1], position[2],
            tangent[0] * arrow_len, tangent[1] * arrow_len, tangent[2] * arrow_len,
            color="#fde047", linewidth=1.6, arrow_length_ratio=0.25,
        )
        time_s = index * float(dt)
        err = position - reference[min(index, len(reference) - 1)]
        hud.set_text(
            f"t={time_s:6.1f}s  step={index:5d}/{completed - 1}\n"
            f"pos=({position[0]:+.2f},{position[1]:+.2f},{position[2]:+.2f})\n"
            f"err=({err[0]:+.2f},{err[1]:+.2f},{err[2]:+.2f})  "
            f"|e|={np.linalg.norm(err):.2f}"
        )
        return (trail_line, future_line, drone_dot, hud, title_artist, legend, ref_line)

    animation = FuncAnimation(
        figure, _update, frames=len(frame_indices), interval=1000.0 / fps,
        blit=False,
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    animation.save(path, writer=PillowWriter(fps=fps))
    plt.close(figure)
    return {
        "path": str(path),
        "frames": int(len(frame_indices)),
        "fps": int(fps),
        "completed_steps": completed,
        "duration_seconds": float(len(frame_indices)) / float(fps),
    }


def render_matlab_style_playback_gif(
        path, trajectory, reference, hard_bound,
        stage_mean_bound=None,
        dt=SAMPLE_TIME_SECONDS,
        max_frames=PLAYBACK_GIF_MAX_FRAMES,
        frame_ms=PLAYBACK_GIF_FRAME_MS,
        title="copyBG SMPC playback"):
    """MATLAB-style fixed-camera trajectory playback → GIF.

    Mimics the usual MATLAB pattern:
      plot3 / hold on / axis equal / view(az,el)  once,
      then loop: update line+marker data, drawnow, getframe,
      finally imwrite(..., 'gif', 'DelayTime', ...).

    Camera NEVER moves — only the path grows. No chase/FPV spin.
    Layout: large 3D + XY top + z(t) (paper-style multi-panel).
    """
    from PIL import Image

    trajectory = np.asarray(trajectory, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if trajectory.ndim != 2 or trajectory.shape[1] != 3:
        raise ValueError("trajectory must be (N,3)")
    completed = int(trajectory.shape[0])
    if completed < 2:
        raise ValueError("need at least 2 trajectory samples for a GIF")
    hard_bound = float(hard_bound)
    max_frames = max(2, int(max_frames))
    frame_ms = max(40, int(frame_ms))

    if completed <= max_frames:
        frame_indices = np.arange(completed, dtype=int)
    else:
        frame_indices = np.unique(
            np.linspace(0, completed - 1, max_frames, dtype=int)
        )

    mean_box = None
    if stage_mean_bound is not None:
        mean_box = float(np.min(np.asarray(stage_mean_bound, dtype=float)[-1]))

    # Fixed axis box from full data + hard bound (MATLAB axis equal style).
    all_pts = np.vstack((reference, trajectory))
    lo = np.minimum(all_pts.min(axis=0), -hard_bound) - 0.35
    hi = np.maximum(all_pts.max(axis=0), hard_bound) + 0.35
    # Symmetric cube-ish limits so 3D doesn't stretch.
    mid = 0.5 * (lo + hi)
    half = 0.5 * float(np.max(hi - lo))
    lim_lo = mid - half
    lim_hi = mid + half
    t_axis = np.arange(completed, dtype=float) * float(dt)
    err_norm = np.linalg.norm(trajectory - reference[:completed], axis=1)

    # Light MATLAB-ish figure (white, grid, fixed view).
    figure = plt.figure(figsize=(10.0, 6.2), facecolor="white", dpi=110)
    grid = figure.add_gridspec(2, 2, width_ratios=[1.55, 1.0],
                               height_ratios=[1.0, 1.0],
                               wspace=0.28, hspace=0.32,
                               left=0.06, right=0.98, top=0.90, bottom=0.08)
    ax3d = figure.add_subplot(grid[:, 0], projection="3d", facecolor="white")
    ax_xy = figure.add_subplot(grid[0, 1], facecolor="white")
    ax_zt = figure.add_subplot(grid[1, 1], facecolor="white")
    figure.suptitle(title, fontsize=12, fontweight="bold", color="#111827")

    # --- 3D panel (fixed camera) ---
    ax3d.plot(reference[:, 0], reference[:, 1], reference[:, 2],
              color="#9ca3af", linestyle="--", linewidth=1.1, label="reference")
    trail_3d, = ax3d.plot([], [], [], color="#2563eb", linewidth=2.0, label="SMPC")
    marker_3d = ax3d.plot([trajectory[0, 0]], [trajectory[0, 1]], [trajectory[0, 2]],
                          "o", color="#dc2626", markersize=7, label="current")[0]
    _draw_bound_box(ax3d, hard_bound, "#ef4444", "-", "hard")
    if mean_box is not None and mean_box > 0.0:
        _draw_bound_box(ax3d, mean_box, "#f59e0b", ":", "chance mean")
    ax3d.set_xlim(lim_lo[0], lim_hi[0])
    ax3d.set_ylim(lim_lo[1], lim_hi[1])
    ax3d.set_zlim(lim_lo[2], lim_hi[2])
    ax3d.set_xlabel("x")
    ax3d.set_ylabel("y")
    ax3d.set_zlabel("z")
    # Fixed MATLAB-like view — never changes across frames.
    ax3d.view_init(elev=22, azim=-52)
    try:
        ax3d.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    ax3d.legend(loc="upper left", fontsize=7, framealpha=0.9)
    ax3d.grid(True, alpha=0.35)

    # --- XY top view ---
    ax_xy.plot(reference[:, 0], reference[:, 1], "--", color="#9ca3af",
               linewidth=1.0, label="ref")
    trail_xy, = ax_xy.plot([], [], color="#2563eb", linewidth=1.8, label="SMPC")
    marker_xy = ax_xy.plot([trajectory[0, 0]], [trajectory[0, 1]], "o",
                           color="#dc2626", markersize=6)[0]
    hard_sq = plt.Rectangle((-hard_bound, -hard_bound), 2 * hard_bound,
                            2 * hard_bound, fill=False, edgecolor="#ef4444",
                            linewidth=1.0, label="hard")
    ax_xy.add_patch(hard_sq)
    if mean_box is not None and mean_box > 0.0:
        mean_sq = plt.Rectangle((-mean_box, -mean_box), 2 * mean_box,
                                2 * mean_box, fill=False, edgecolor="#f59e0b",
                                linestyle=":", linewidth=1.0)
        ax_xy.add_patch(mean_sq)
    ax_xy.set_xlim(lim_lo[0], lim_hi[0])
    ax_xy.set_ylim(lim_lo[1], lim_hi[1])
    ax_xy.set_aspect("equal", adjustable="box")
    ax_xy.set_xlabel("x")
    ax_xy.set_ylabel("y")
    ax_xy.set_title("XY top", fontsize=10)
    ax_xy.grid(True, alpha=0.4)
    ax_xy.legend(loc="upper right", fontsize=7)

    # --- z(t) + |e|(t) ---
    ax_zt.plot(t_axis, reference[:completed, 2], "--", color="#9ca3af",
               linewidth=1.0, label="z_ref")
    trail_z, = ax_zt.plot([], [], color="#2563eb", linewidth=1.6, label="z")
    trail_e, = ax_zt.plot([], [], color="#059669", linewidth=1.2, label="|e|")
    marker_t = ax_zt.axvline(0.0, color="#dc2626", linewidth=1.0, alpha=0.85)
    ax_zt.set_xlim(0.0, float(t_axis[-1]))
    zpad = 0.15 * max(1.0, float(np.ptp(reference[:completed, 2])))
    ax_zt.set_ylim(float(reference[:completed, 2].min()) - zpad,
                   float(max(reference[:completed, 2].max(),
                             err_norm.max())) + zpad)
    ax_zt.set_xlabel("t [s]")
    ax_zt.set_ylabel("z / |e|")
    ax_zt.set_title("altitude & tracking error", fontsize=10)
    ax_zt.grid(True, alpha=0.4)
    ax_zt.legend(loc="upper right", fontsize=7)

    hud = figure.text(
        0.06, 0.02, "", fontsize=9, family="monospace", color="#111827",
        ha="left", va="bottom",
    )

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frames = []

    for frame_no, index in enumerate(frame_indices):
        index = int(index)
        pos = trajectory[index]
        trail_3d.set_data(trajectory[: index + 1, 0], trajectory[: index + 1, 1])
        trail_3d.set_3d_properties(trajectory[: index + 1, 2])
        marker_3d.set_data([pos[0]], [pos[1]])
        marker_3d.set_3d_properties([pos[2]])

        trail_xy.set_data(trajectory[: index + 1, 0], trajectory[: index + 1, 1])
        marker_xy.set_data([pos[0]], [pos[1]])

        trail_z.set_data(t_axis[: index + 1], trajectory[: index + 1, 2])
        trail_e.set_data(t_axis[: index + 1], err_norm[: index + 1])
        marker_t.set_xdata([t_axis[index], t_axis[index]])

        e = pos - reference[min(index, len(reference) - 1)]
        hud.set_text(
            f"t={t_axis[index]:6.1f}s   step={index:5d}/{completed - 1}   "
            f"pos=({pos[0]:+.2f},{pos[1]:+.2f},{pos[2]:+.2f})   "
            f"|e|={np.linalg.norm(e):.3f}"
        )

        # Keep fixed view every frame (matplotlib sometimes resets).
        ax3d.view_init(elev=22, azim=-52)
        figure.canvas.draw()
        rgba = np.asarray(figure.canvas.buffer_rgba())
        frames.append(Image.fromarray(rgba[:, :, :3].copy()))

    plt.close(figure)

    if not frames:
        raise RuntimeError("no frames captured for playback GIF")
    # Explicit DelayTime like MATLAB imwrite(..., 'DelayTime', frame_ms/1000)
    frames[0].save(
        path,
        save_all=True,
        append_images=frames[1:],
        duration=frame_ms,
        loop=0,
        optimize=False,
        disposal=2,
    )
    play_s = float(len(frames) * frame_ms) / 1000.0
    return {
        "path": str(path),
        "style": "matlab_fixed_camera_playback",
        "frames": int(len(frames)),
        "frame_ms": int(frame_ms),
        "fps": float(1000.0 / frame_ms),
        "play_seconds": play_s,
        "completed_steps": completed,
        "view": {"elev": 22, "azim": -52},
    }


def _npz_payload(reference, smpc, deterministic_mpc, residual_innovation,
                 wind_velocity_mps, wind_latent_increment, covariance,
                 noise, initialization, hard_bound):
    return {
        "reference_standardized": np.asarray(reference, dtype=float),
        "smpc_trajectory_standardized": np.asarray(smpc["S"], dtype=float),
        "deterministic_mpc_trajectory_standardized": np.asarray(
            deterministic_mpc["S"], dtype=float,
        ),
        "smpc_control": np.asarray(smpc["U"], dtype=float),
        "deterministic_mpc_control": np.asarray(
            deterministic_mpc["U"], dtype=float,
        ),
        "innovation_latent": np.asarray(residual_innovation, dtype=float),
        "wind_velocity_mps": np.asarray(wind_velocity_mps, dtype=float),
        "wind_position_delta_xyz_m": (
            np.asarray(wind_velocity_mps, dtype=float) * SAMPLE_TIME_SECONDS
        ),
        "wind_latent_increment": np.asarray(wind_latent_increment, dtype=float),
        "Sigma_eps": np.asarray(noise["Sigma_eps"], dtype=float),
        "Sigma_obs_proxy": np.asarray(noise["Sigma_obs_proxy"], dtype=float),
        "Sigma_eps_aug": np.asarray(covariance["Sigma_eps_aug"], dtype=float),
        "Sigma_wind_aug": np.asarray(covariance["Sigma_wind_aug"], dtype=float),
        "Sigma_total_aug": np.asarray(covariance["Sigma_total_aug"], dtype=float),
        "G_w": np.asarray(covariance["G_w"], dtype=float),
        "hard_bound_standardized": np.asarray([hard_bound], dtype=float),
        "stage_tightening_standardized": np.asarray(
            smpc["stage_tightening"], dtype=float,
        ),
        "stage_mean_bound_standardized": np.asarray(
            smpc["stage_mean_bound"], dtype=float,
        ),
        "smpc_active_history": np.asarray(smpc["active_history"], dtype=np.uint8),
        "deterministic_mpc_active_history": np.asarray(
            deterministic_mpc["active_history"], dtype=np.uint8,
        ),
        "smpc_qp_constraint_residual_history": np.asarray(
            smpc["maximum_qp_constraint_residual_history"], dtype=float,
        ),
        "deterministic_mpc_qp_constraint_residual_history": np.asarray(
            deterministic_mpc["maximum_qp_constraint_residual_history"], dtype=float,
        ),
        "initial_state": np.asarray(initialization["state"], dtype=float),
        "initial_control": np.asarray(initialization["control"], dtype=float),
        "dt_used_seconds": np.asarray([SAMPLE_TIME_SECONDS], dtype=float),
        "shared_noise_exact": np.asarray([1], dtype=np.uint8),
    }


def run_experiment(steps=None, sigma_wind_mps=SIGMA_WIND_MPS,
                   hard_bound=None,
                   pressures=None, reference_mode=DEFAULT_REFERENCE_MODE,
                   rising_turns=None,
                   output_dir=None,
                   make_forward_gif=True):
    if steps is None:
        steps = (
            SMOOTH_SQUARE_DEFAULT_STEPS
            if reference_mode == "smooth-square-helix"
            else GROWING_SMOOTH_SQUARE_DEFAULT_STEPS
            if reference_mode == "growing-smooth-square-helix"
            else DEFAULT_STEPS
        )
    steps = int(steps)
    if not 40 <= steps <= MAX_STEPS:
        raise ValueError(f"steps must be between 40 and {MAX_STEPS}")
    if reference_mode not in REFERENCE_MODES:
        raise ValueError(
            f"reference_mode must be one of {REFERENCE_MODES}, got {reference_mode!r}"
        )
    if rising_turns is None:
        if reference_mode == "smooth-square-helix":
            rising_turns = SMOOTH_SQUARE_RISING_TURNS
        elif reference_mode == "growing-smooth-square-helix":
            rising_turns = GROWING_SMOOTH_SQUARE_RISING_TURNS
        else:
            rising_turns = CIRCLE_HELIX_RISING_TURNS
    rising_turns = float(rising_turns)
    if rising_turns <= 0.0:
        raise ValueError("rising_turns must be positive")

    model, E, source_hard_bound, scales, noise = (
        load_identification_and_noise_objects()
    )
    # Auto hard bound: large boundary-riding modes use the identification
    # source box (~3.91, copyBF recipe); small helix modes keep 3.4.
    if hard_bound is None:
        if reference_mode in (
            "boundary-tour",
            "smooth-square-helix",
            "growing-smooth-square-helix",
        ):
            hard_bound = float(source_hard_bound)
        else:
            hard_bound = float(HARD_BOUND_STANDARDIZED)
    else:
        hard_bound = float(hard_bound)

    if pressures is None:
        if reference_mode == "smooth-square-helix":
            pressures = SMOOTH_SQUARE_HELIX_PRESSURES
        elif reference_mode == "growing-smooth-square-helix":
            pressures = GROWING_SMOOTH_SQUARE_END
        elif reference_mode == "square-spiral":
            pressures = SQUARE_SPIRAL_PRESSURES
        elif reference_mode == "circle-helix":
            pressures = CIRCLE_HELIX_PRESSURES
        elif reference_mode == "boundary-tour":
            pressures = BOUNDARY_TOUR_TARGETS
        else:
            pressures = REFERENCE_PRESSURES_STANDARDIZED
    pressures = np.asarray(pressures, dtype=float).reshape(3)
    if np.any(pressures >= hard_bound):
        raise ValueError("reference pressures must remain inside the hard bound")

    dc_feasibility = None
    growth_meta = None
    if reference_mode == "smooth-square-helix":
        end_for_dc = pressures
        dc_feasibility = smooth_square_dc_feasibility(
            model, E, pressures=end_for_dc,
        )
        if not dc_feasibility["feasible"]:
            raise ValueError(
                f"{reference_mode} end pressures are DC-infeasible under the implemented ±6 inputs: "
                f"worst u_inf={dc_feasibility['worst_u_inf']:.3f} at task "
                f"{dc_feasibility['worst_task'].tolist()} "
                f"(cap={dc_feasibility['u_inf_cap']}). "
                "Lower --pressures or use growing-smooth-square-helix / "
                "boundary-tour for pressure geometry. This is an experiment "
                "calibration gate, not a CRTE/copyAU change."
            )
    elif reference_mode == "growing-smooth-square-helix":
        # Sustained paths must satisfy the implemented ±6 standardized-command cap.
        if GROWING_PATH_FAMILY == "z-face-large-arc-cruise":
            dc_feasibility = z_face_large_arc_cruise_dc_feasibility(
                model, E,
                xy_radius=float(pressures[0]),
                z_face=float(pressures[2]),
                u_inf_cap=GROWING_DC_U_INF_CAP,
            )
        else:
            dc_feasibility = large_arc_diamond_dc_feasibility(
                model, E, pressures=pressures,
                center_offset=GROWING_ARC_CENTER_OFFSET,
                u_inf_cap=GROWING_DC_U_INF_CAP,
            )

    if reference_mode == "smooth-square-helix":
        reference = smooth_square_helix_reference(
            steps, pressures, rising_turns=rising_turns,
        )
        reference_path_desc = (
            "smooth square helix: bottom dwell, rising n=4 superellipse "
            f"(rounded-box) turns={rising_turns} dwell={SMOOTH_SQUARE_DWELL_FRACTION} "
            f"at DC-feasible pressures {tuple(float(v) for v in pressures)} "
            f"(corner u_inf≈{dc_feasibility['worst_u_inf']:.2f}), "
            f"top dwell — multi-axis orbit ({steps} samples, arc-length paced)"
        )
    elif reference_mode == "growing-smooth-square-helix":
        if np.any(pressures >= hard_bound):
            raise ValueError(
                "growing peak pressures must remain inside the hard bound"
            )
        reference = growing_smooth_square_helix_reference(
            steps=steps,
            pressures_end=pressures,
            pressures_start=GROWING_SMOOTH_SQUARE_START,
            rising_turns=rising_turns,
        )
        phase_ids, phase_counts = growing_envelope_phase_index(steps)
        phases = growing_smooth_square_phase_masks(
            reference,
            r_dc=GROWING_DC_LABEL_R,
            hard_bound=hard_bound,
            phase_ids=phase_ids,
            peak_pressures=pressures,
        )
        phase_names = list(phases.get("phase_names", ()))
        growth_meta = {
            "structure": (
                "continuous horizontal large-arc cruise on upper/lower z hard faces "
                "(approach→top cruise→transfer→bottom cruise→return)"
            ),
            "path_family": GROWING_PATH_FAMILY,
            "arc_center_offset": dc_feasibility.get("center_offset"),
            "max_joint_xy_unit": float(dc_feasibility.get("max_joint_xy_unit", np.nan)),
            "pressures_start": np.asarray(GROWING_SMOOTH_SQUARE_START, dtype=float),
            "pressures_peak": np.asarray(pressures, dtype=float),
            "z_face_absolute": float(pressures[2]),
            "horizontal_large_arc_radius": float(pressures[0]),
            "phase_turns": list(float(v) for v in GROWING_Z_FACE_TURNS),
            "phase_fractions": list(GROWING_PHASE_FRACTIONS),
            "phase_counts": [int(c) for c in phase_counts],
            "phase_names": phase_names,
            "r_xy_max": phases["r_xy_max"],
            "r_box_max": phases["r_box_max"],
            "r_linf_max": phases["r_linf_max"],
            "fraction_dc_feasible_samples": float(np.mean(phases["dc_feasible_mask"])),
            "fraction_beyond_dc_samples": float(np.mean(phases["beyond_dc_mask"])),
            "fraction_top_cruise": float(phases.get("fraction_top_cruise", 0.0)),
            "fraction_bottom_cruise": float(phases.get("fraction_bottom_cruise", 0.0)),
            "near_hard_fraction": float(np.mean(phases.get("near_hard_mask", [False]))),
            "near_hard_band": float(GROWING_NEAR_HARD_BAND),
            "near_hard_target_fraction": float(GROWING_NEAR_HARD_TARGET_FRACTION),
            "peak_dc_feasible": bool(dc_feasibility["feasible"]),
            "peak_path_u_inf": float(dc_feasibility["worst_u_inf"]),
            "peak_dc_note": (
                "both moving z-face cruises fit the implemented standardized-command cap "
                f"(minimax u_inf≈{dc_feasibility['worst_u_inf']:.2f} <= {GROWING_DC_U_INF_CAP}); "
                "the horizontal path keeps flying while z remains near task hard"
                if dc_feasibility["feasible"]
                else (
                    f"z-face cruise exceeds the standardized-command cap "
                    f"(u_inf≈{dc_feasibility['worst_u_inf']:.2f}); "
                    "do not accept this reference"
                )
            ),
        }
        reference_path_desc = (
            "z-face large-arc near-hard cruise: "
            f"start {tuple(float(v) for v in GROWING_SMOOTH_SQUARE_START)} → "
            f"z=+{pressures[2]:g} moving cruise → z=-{pressures[2]:g} moving cruise "
            f"(family={GROWING_PATH_FAMILY}, xy radius={pressures[0]:g}, "
            f"turns={tuple(float(v) for v in GROWING_Z_FACE_TURNS)}, "
            f"minimax worst u_inf≈{dc_feasibility['worst_u_inf']:.2f}, "
            f"dc_feasible={dc_feasibility['feasible']}, "
            f"max_joint_xy_unit≈{dc_feasibility.get('max_joint_xy_unit', float('nan')):.3f}) "
            f"→ return; {steps} samples — majority reference time near task hard"
        )
    elif reference_mode == "square-spiral":
        reference = square_helix_reference(
            steps, pressures, rising_turns=rising_turns,
        )
        reference_path_desc = (
            "square helix: bottom dwell, rising L-infinity square turns, "
            "top dwell (optional sustained corner-pressure mode)"
        )
    elif reference_mode == "circle-helix":
        reference = circle_helix_reference(
            steps, pressures, rising_turns=rising_turns,
        )
        reference_path_desc = (
            "circle helix: bottom dwell, rising circular/elliptical turns "
            "(rx cos θ, ry sin θ), top dwell — cylindrical tracking path"
        )
    elif reference_mode == "boundary-tour":
        reference = continuous_boundary_tour(steps, pressures)
        reference_path_desc = (
            "copyBF-style continuous boundary tour: C1 x→y→z→x axial legs "
            f"at targets {tuple(float(v) for v in pressures)}"
        )
    else:
        reference = continuous_near_boundary_reference(steps, pressures)
        reference_path_desc = "center -> +x -> +y -> +z -> center, C1 smooth"

    residual_innovation = innovation_sequence(model, steps, INNOVATION_SEED)
    wind_velocity_mps = physical_wind_velocity(
        steps, sigma_wind_mps, WIND_SEED,
    )
    task_map = E.T @ model["P"][:, :model["R"].shape[1]]
    wind_latent_increment, G_w = wind_velocity_to_latent_increment(
        wind_velocity_mps, task_map, scales["y_scale"][:3],
    )
    covariance = build_process_covariances(
        model, E, scales["y_scale"][:3], sigma_wind_mps,
    )
    if not np.allclose(G_w, covariance["G_w"], atol=1e-14):
        raise RuntimeError("Wind realization and covariance use different G_w maps")

    # Feasible start: large boundary anchors may have no input-bounded
    # equilibrium. Match copyBF — park just inside the terminal mean bound.
    initial_task = np.asarray(reference[0], dtype=float).reshape(3)
    try:
        initialization = initialization_diagnostics(model, E, initial_task)
    except ValueError:
        probe = CachedBoundaryController(
            model, E, hard_bound, covariance["Sigma_eps_aug"],
            horizon=HORIZON_STEPS, alpha_joint=ALPHA_JOINT,
            control_interval=CONTROL_INTERVAL_STEPS,
        )
        terminal_mean = np.asarray(probe.stage_mean_bound[-1], dtype=float)
        initial_task = np.array([
            float(terminal_mean[0] - 0.05),
            0.0,
            0.0,
        ])
        if initial_task[0] <= 0.0:
            raise RuntimeError(
                "Terminal mean bound too tight for a feasible boundary start"
            )
        initialization = initialization_diagnostics(model, E, initial_task)

    smpc, deterministic_mpc = run_controller_pair(
        model, E, hard_bound, reference, residual_innovation,
        wind_latent_increment, covariance["Sigma_eps_aug"],
        initial_task=initial_task,
    )
    if smpc["fallback_count"] or deterministic_mpc["fallback_count"]:
        raise RuntimeError("Fallback is forbidden in copyBG")
    if smpc["maximum_qp_constraint_residual"] > 1e-6:
        raise RuntimeError("SMPC QP constraint residual exceeds 1e-6")
    if deterministic_mpc["maximum_qp_constraint_residual"] > 1e-6:
        raise RuntimeError("Deterministic MPC QP constraint residual exceeds 1e-6")

    if growth_meta is not None:
        growth_meta["reference_near_hard_occupancy"] = near_hard_occupancy(
            reference, hard_bound,
        )
        growth_meta["smpc_near_hard_occupancy"] = near_hard_occupancy(
            smpc["S"], hard_bound,
        )
        growth_meta["deterministic_mpc_near_hard_occupancy"] = near_hard_occupancy(
            deterministic_mpc["S"], hard_bound,
        )
        growth_meta["smpc_near_hard_target_met"] = bool(
            growth_meta["smpc_near_hard_occupancy"]["fraction"] >= 0.50
            and growth_meta["smpc_near_hard_occupancy"]["violation_steps"] == 0
        )

    results_dir = Path(output_dir) if output_dir else HERE / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    stem = "copyBG_pelican_physical_wind_smpc"
    paths = {
        "json": results_dir / f"{stem}.json",
        "npz": results_dir / f"{stem}.npz",
        "mat": results_dir / f"{stem}.mat",
        "smpc_png": results_dir / f"{stem}_smpc.png",
        "deterministic_mpc_png": results_dir / f"{stem}_deterministic_mpc.png",
    }
    plot_controller_diagnostics(
        paths["smpc_png"], "SMPC", smpc, reference, hard_bound,
        smpc["stage_mean_bound"], sigma_wind_mps,
    )
    plot_controller_diagnostics(
        paths["deterministic_mpc_png"], "deterministic MPC",
        deterministic_mpc, reference, hard_bound,
        smpc["stage_mean_bound"], sigma_wind_mps,
    )
    paths["noise_compare_png"] = (
        results_dir / f"{stem}_noise_actual_vs_estimated.png"
    )
    noise_compare_meta = plot_noise_actual_vs_estimated(
        paths["noise_compare_png"],
        residual_innovation=residual_innovation,
        wind_velocity_mps=wind_velocity_mps,
        wind_latent_increment=wind_latent_increment,
        Sigma_eps=noise["Sigma_eps"],
        Sigma_wind_aug=covariance["Sigma_wind_aug"],
        Sigma_total_aug=covariance["Sigma_total_aug"],
        sigma_wind_mps=sigma_wind_mps,
        dt_seconds=SAMPLE_TIME_SECONDS,
    )
    gif_meta = None
    if make_forward_gif and int(smpc["completed_steps"]) >= 200:
        # Default product GIF: MATLAB-style fixed-camera playback (not chase).
        paths["smpc_playback_gif"] = results_dir / f"{stem}_smpc_playback.gif"
        gif_meta = render_matlab_style_playback_gif(
            paths["smpc_playback_gif"],
            trajectory=smpc["S"],
            reference=reference,
            hard_bound=hard_bound,
            stage_mean_bound=smpc["stage_mean_bound"],
            dt=SAMPLE_TIME_SECONDS,
            title=(
                f"copyBG SMPC playback | {reference_mode} "
                f"turns={rising_turns:g} wind={sigma_wind_mps:.3f}"
            ),
        )
        # Also keep a copy under the older forward filename so verify/docs
        # that expect a gif in results/ still find one.
        paths["smpc_forward_gif"] = results_dir / f"{stem}_smpc_forward.gif"
        try:
            import shutil
            shutil.copy2(paths["smpc_playback_gif"], paths["smpc_forward_gif"])
        except Exception:
            pass
    np.savez_compressed(
        paths["npz"], **_npz_payload(
            reference, smpc, deterministic_mpc, residual_innovation,
            wind_velocity_mps, wind_latent_increment, covariance,
            noise, initialization, hard_bound,
        )
    )
    mat_payload = build_mat_payload(
        reference, smpc, deterministic_mpc, residual_innovation,
        wind_velocity_mps, wind_latent_increment, covariance, noise,
        scales, hard_bound, smpc["stage_tightening"], initialization,
    )
    scipy.io.savemat(paths["mat"], mat_payload, do_compression=True)

    summary = {
        "experiment": "copyBG_pelican_physical_wind_smpc",
        "validation_type": "Pelican-data-identified model-in-the-loop",
        "is_real_flight_validation": False,
        "requested_steps": steps,
        "sample_rate_hz": 1.0 / SAMPLE_TIME_SECONDS,
        "dt_used_seconds": SAMPLE_TIME_SECONDS,
        "control_interval_steps": CONTROL_INTERVAL_STEPS,
        "horizon_steps": HORIZON_STEPS,
        "horizon_seconds": (
            HORIZON_STEPS * CONTROL_INTERVAL_STEPS * SAMPLE_TIME_SECONDS
        ),
        "Q_weight": Q_WEIGHT,
        "Ru": RU,
        "alpha_joint": ALPHA_JOINT,
        "input_bound_normalized": [-6.0, 6.0],
        "hard_bound_standardized": hard_bound,
        "source_model_hard_bound_standardized": source_hard_bound,
        "reference_mode": reference_mode,
        "reference_pressures_standardized": pressures,
        "rising_turns": rising_turns,
        "smooth_square_dwell_fraction": (
            SMOOTH_SQUARE_DWELL_FRACTION
            if reference_mode == "smooth-square-helix"
            else GROWING_SMOOTH_SQUARE_DWELL_FRACTION
            if reference_mode == "growing-smooth-square-helix"
            else None
        ),
        "reference_path": reference_path_desc,
        "growing_smooth_square": (
            None if growth_meta is None else _json_ready(growth_meta)
        ),
        "smooth_square_dc_feasibility": (
            None
            if dc_feasibility is None
            else _json_ready({
                "pressures": dc_feasibility["pressures"],
                "family": dc_feasibility.get("family", "superellipse"),
                "corner_power": dc_feasibility.get("corner_power"),
                "center_offset": dc_feasibility.get("center_offset"),
                "max_joint_xy_unit": dc_feasibility.get("max_joint_xy_unit"),
                "u_inf_cap": dc_feasibility["u_inf_cap"],
                "worst_u_inf": dc_feasibility["worst_u_inf"],
                "worst_task": dc_feasibility["worst_task"],
                "worst_u": dc_feasibility["worst_u"],
                "feasible": dc_feasibility["feasible"],
                "note": dc_feasibility["note"],
            })
        ),
        "forward_gif": gif_meta,
        "playback_gif": gif_meta,
        "sigma_wind_mps": sigma_wind_mps,
        "position_unit": "m",
        "wind_seed": WIND_SEED,
        "innovation_seed": INNOVATION_SEED,
        "training_segments": noise["training_segments"],
        "training_samples": noise["training_samples"],
        "transition_count": noise["transition_count"],
        "uses_true_Sigma_n": noise["uses_true_Sigma_n"],
        "obs_tracked_leakage_norm": noise["obs_tracked_leakage_norm"],
        "Sigma_eps_trace": float(np.trace(noise["Sigma_eps"])),
        "Sigma_obs_proxy_trace": float(np.trace(noise["Sigma_obs_proxy"])),
        "Sigma_wind_aug_trace": float(np.trace(covariance["Sigma_wind_aug"])),
        "Sigma_total_aug_trace": float(np.trace(covariance["Sigma_total_aug"])),
        "noise_actual_vs_estimated": _json_ready(noise_compare_meta),
        "shared_noise_exact": 1,
        "disturbance_sha256": smpc["disturbance_sha256"],
        "initialization": _json_ready({
            "dynamic_residual": initialization["dynamic_residual"],
            "task_residual": initialization["task_residual"],
            "maximum_absolute_input": initialization["maximum_absolute_input"],
            "control": initialization["control"],
        }),
        "smpc": _controller_summary(smpc),
        "deterministic_mpc": _controller_summary(deterministic_mpc),
        "smpc_vs_mpc": {
            "active_qp_steps_smpc": int(smpc["active_qp_steps"]),
            "active_qp_steps_mpc": int(deterministic_mpc["active_qp_steps"]),
            "positive_dual_qp_steps_smpc": int(smpc["positive_dual_qp_steps"]),
            "positive_dual_qp_steps_mpc": int(
                deterministic_mpc["positive_dual_qp_steps"]
            ),
            "hard_violation_steps_smpc": int(smpc["violation_steps"]),
            "hard_violation_steps_mpc": int(
                deterministic_mpc["violation_steps"]
            ),
            "minimum_hard_margin_smpc": float(smpc["minimum_hard_margin"]),
            "minimum_hard_margin_mpc": float(
                deterministic_mpc["minimum_hard_margin"]
            ),
            "rmse_smpc": np.asarray(smpc["rmse"], dtype=float),
            "rmse_mpc": np.asarray(deterministic_mpc["rmse"], dtype=float),
            "chance_expected_active": bool(
                reference_mode in (
                    "boundary-tour",
                    "tour",
                    "growing-smooth-square-helix",
                )
            ),
            "note": (
                "Chance constraints bind when the planned trajectory presses "
                "stage_mean_bound. Default v4.0 mode flies horizontal large arcs "
                "on z=±3.57 for a majority of the run; its steady input demand "
                "is checked against the inherited normalized ±6 experiment cap, "
                "not a claimed Pelican hardware rating. "
                "Pure boundary-tour remains available for single-axis legs."
            ),
        },
        "interpretation": (
            "The deterministic MPC trace ends at its first non-optimal QP; "
            "no fallback or post-failure trajectory is synthesized. "
            "Default narrative is approach→upper z-face cruise→transfer→lower "
            "z-face cruise→return. Near-hard occupancy is counted only while "
            "inside the task hard box. SMPC value is read from active_qp / "
            "duals / hard-margin vs deterministic MPC."
        ),
        "artifacts": {key: path.name for key, path in paths.items()},
    }
    paths["json"].write_text(
        json.dumps(_json_ready(summary), indent=2, allow_nan=False),
        encoding="ascii",
    )
    return summary


def _main():
    parser = argparse.ArgumentParser(
        description="Run copyBG Pelican physical-wind SMPC/MPC comparison.",
    )
    parser.add_argument(
        "--steps", type=int, default=None,
        help=(
            f"closed-loop samples @ {1.0 / SAMPLE_TIME_SECONDS:.0f} Hz; "
            f"default {GROWING_SMOOTH_SQUARE_DEFAULT_STEPS} for "
            f"growing-smooth-square-helix (default mode), "
            f"{SMOOTH_SQUARE_DEFAULT_STEPS} for smooth-square-helix, "
            f"{DEFAULT_STEPS} for boundary-tour (max {MAX_STEPS}). "
            "Longer = slower path on the same geometry."
        ),
    )
    parser.add_argument("--sigma-wind-mps", type=float, default=SIGMA_WIND_MPS)
    parser.add_argument(
        "--hard-bound", type=float, default=None,
        help=(
            "hard box in standardized task coords; default auto = source "
            "model bound for growing/smooth-square/boundary-tour, else "
            f"{HARD_BOUND_STANDARDIZED}"
        ),
    )
    parser.add_argument(
        "--reference", choices=REFERENCE_MODES,
        default=DEFAULT_REFERENCE_MODE,
        help=(
            "reference path: growing-smooth-square-helix (default: "
            "moving large-arc cruises on the upper/lower z hard faces), "
            "boundary-tour (single-axis SMPC legs), "
            "smooth-square-helix (fixed multi-axis DC box), circle-helix, "
            "square-spiral, or center-return tour"
        ),
    )
    parser.add_argument(
        "--rising-turns", type=float, default=None,
        help=(
            "helix revolutions per grow/shrink leg (growing) or rising "
            f"segment; default {GROWING_SMOOTH_SQUARE_RISING_TURNS} (growing) / "
            f"{SMOOTH_SQUARE_RISING_TURNS} (smooth-square) / "
            f"{CIRCLE_HELIX_RISING_TURNS} (other helices)"
        ),
    )
    parser.add_argument(
        "--pressures", type=float, nargs=3, default=None,
        metavar=("PX", "PY", "PZ"),
        help="reference pressures; default per mode",
    )
    parser.add_argument("--output-dir", type=Path, default=HERE / "results")
    arguments = parser.parse_args()
    summary = run_experiment(
        steps=arguments.steps,
        sigma_wind_mps=arguments.sigma_wind_mps,
        hard_bound=arguments.hard_bound,
        pressures=arguments.pressures,
        reference_mode=arguments.reference,
        rising_turns=arguments.rising_turns,
        output_dir=arguments.output_dir,
    )
    print(json.dumps({
        "json": summary["artifacts"]["json"],
        "smpc_completed_steps": summary["smpc"]["completed_steps"],
        "smpc_active_qp_steps": summary["smpc"]["active_qp_steps"],
        "smpc_hard_violation_rate": summary["smpc"]["hard_violation_rate"],
        "deterministic_mpc_completed_steps": summary[
            "deterministic_mpc"
        ]["completed_steps"],
        "deterministic_mpc_qp_failure_step": summary[
            "deterministic_mpc"
        ]["qp_failure_step"],
    }, indent=2))


if __name__ == "__main__":
    _main()
