"""Continuous 120 s spiral approaching each axis-specific tightened bound.

The external wind channel is a per-step i.i.d. Gaussian process noise
``w_k ~ N(0, sigma_w^2 I_ell)`` injected in latent coordinates, in the same
spirit as copyT process-noise injection.

Two independent random streams feed the closed loop:

* ``residual_noise ~ N(0, Sigma_eps[:ell,:ell])`` -- the latent innovation
  estimated from the Pelican dataset. This drives the controller's only
  accuracy-statement about Sigma_eps; it is not a measurement-noise
  description and is not claimed to equal any specific scalar sigma_e.
  SMPC and deterministic MPC both see the exact same realization.
* ``process_wind`` -- the per-step i.i.d. Gaussian disturbance injected in
  latent coordinates. Sampled with a deliberately independent seed so its
  base RNG sequence is not the residual stream's.

Both controllers receive the same complete ``residual_noise`` and
``process_wind`` arrays; only the controller differs.
"""

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io


HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "run_long_3d_boundary_disturbance.py"
_spec = importlib.util.spec_from_file_location("copybf_base", BASE_PATH)
BF = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(BF)


# External wind channel: independent Gaussian process noise in latent
# coordinates. Driven by a deliberately distinct RNG seed so it does not
# share the residual stream's underlying sequence.
SIGMA_W = 0.045
WIND_SEED = 7007

# Identified residual innovation from the Pelican dataset. These two seeds
# are intentionally different -- they control two statistically independent
# RNG streams that share realization across SMPC vs deterministic MPC but
# not across (residual_noise, process_wind) themselves.
RESIDUAL_SEED = 7


def tightened_boundary_spiral(tightened_bound, steps=12000,
                              radial_margin=0.08, vertical_margin=0.08,
                              rx=None, ry=None, z_lo=None, z_hi=None):
    """Bottom dwell, rising elliptical helix, and top dwell without resets.

    Defaults follow the terminal tightened box. Explicit rx/ry/z_lo/z_hi let
    the caller account for model reachability (input saturation limits the
    xy-plane ellipse while z can approach its tightened bound).
    """
    tightened_bound = np.asarray(tightened_bound, dtype=float)
    steps = int(steps)
    bottom_steps = steps // 5
    top_steps = steps // 5
    rise_steps = steps - bottom_steps - top_steps
    if rx is None:
        rx = tightened_bound[0] - float(radial_margin)
    if ry is None:
        ry = tightened_bound[1] - float(radial_margin)
    if z_lo is None:
        z_lo = -(tightened_bound[2] - float(vertical_margin))
    if z_hi is None:
        z_hi = tightened_bound[2] - float(vertical_margin)
    z_lo = float(z_lo)
    z_hi = float(z_hi)

    # One bottom turn, three rising turns, one top turn. Endpoint handling keeps
    # position continuous at both joins and includes all axis extrema.
    theta_bottom = np.linspace(0.0, 2.0 * np.pi, bottom_steps, endpoint=False)
    theta_rise = np.linspace(2.0 * np.pi, 8.0 * np.pi, rise_steps, endpoint=False)
    theta_top = np.linspace(8.0 * np.pi, 10.0 * np.pi, top_steps, endpoint=True)
    theta = np.concatenate((theta_bottom, theta_rise, theta_top))
    z = np.concatenate((
        np.full(bottom_steps, z_lo),
        np.linspace(z_lo, z_hi, rise_steps, endpoint=False),
        np.full(top_steps, z_hi),
    ))
    return np.column_stack((rx * np.cos(theta), ry * np.sin(theta), z))


def process_wind_noise(steps, ell, sigma_w=SIGMA_W, seed=WIND_SEED):
    """Per-step independent Gaussian process wind in latent coordinates.

    ``w_k ~ N(0, sigma_w^2 I_ell)`` for each k, sampled with the fixed seed
    so the closed-loop response is reproducible. This mirrors the copyT
    process-noise channel; it is a separate RNG stream from the residual
    innovation because ``WIND_SEED != RESIDUAL_SEED``.
    """
    steps = int(steps)
    ell = int(ell)
    rng = np.random.default_rng(int(seed))
    return sigma_w * rng.standard_normal(size=(steps, ell))


def residual_innovation_sequence(model, ell, steps, seed=RESIDUAL_SEED):
    """Per-step latent innovation N(0, Sigma_eps[:ell,:ell]) from Pelican data.

    The covariance is the lifted innovation covariance estimated from the
    Pelican dataset. This stream feeds the closed loop independently of the
    external wind channel and does NOT represent measurement noise.
    """
    steps = int(steps)
    ell = int(ell)
    rng = np.random.default_rng(int(seed))
    return rng.multivariate_normal(
        np.zeros(ell), model["Sigma_eps"][:ell, :ell], size=steps
    )


def _core_metric_bundle(result, controller_name):
    """Return the flat scalar metrics required in MAT/NPZ exports."""
    return {
        f"{controller_name}_violation_steps": int(result["violation_steps"]),
        f"{controller_name}_minimum_hard_margin": float(
            result["minimum_hard_margin"]
        ),
        f"{controller_name}_fallback_count": int(result["fallback_count"]),
        f"{controller_name}_input_saturation_steps": int(
            result["input_saturation_steps"]
        ),
        f"{controller_name}_active_qp_steps": int(result["active_qp_steps"]),
        f"{controller_name}_rmse_x": float(result["rmse"][0]),
        f"{controller_name}_rmse_y": float(result["rmse"][1]),
        f"{controller_name}_rmse_z": float(result["rmse"][2]),
    }


def build_mat_payload(reference, smpc, mpc, noise, wind, hard_bound,
                      tightened_bound, sigma_w, seed, wind_seed):
    """Compose the MATLAB export payload.

    Pure helper: no filesystem side effects, no simulation. MAT arrays are
    cast to float/int shapes that ``scipy.io.savemat`` accepts. The contract
    includes both controllers' core metrics, the disturbance realizations,
    and provenance (seed + wind_seed).
    """
    return {
        "reference": np.asarray(reference, dtype=float),
        "smpc_trajectory": np.asarray(smpc["S"], dtype=float),
        "smpc_control": np.asarray(smpc["U"], dtype=float),
        "deterministic_mpc_trajectory": np.asarray(mpc["S"], dtype=float),
        "deterministic_mpc_control": np.asarray(mpc["U"], dtype=float),
        "residual_noise": np.asarray(noise, dtype=float),
        "process_wind": np.asarray(wind, dtype=float),
        "hard_bound": np.asarray([float(hard_bound)], dtype=float),
        "tightened_bound": np.asarray(tightened_bound, dtype=float),
        "sigma_w": np.asarray([float(sigma_w)], dtype=float),
        "seed": np.asarray([int(seed)], dtype=float),
        "wind_seed": np.asarray([int(wind_seed)], dtype=float),
        **_core_metric_bundle(smpc, "smpc"),
        **_core_metric_bundle(mpc, "deterministic_mpc"),
    }


def save_mat_results(path, payload):
    """Persist ``payload`` to ``path`` as a compressed .mat file."""
    scipy.io.savemat(str(path), payload, do_compression=True)


def _draw_box(ax, limits, color, style, label):
    limits = np.asarray(limits, dtype=float)
    corners = np.array([
        [sx * limits[0], sy * limits[1], sz * limits[2]]
        for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)
    ])
    first = True
    for i, a in enumerate(corners):
        for j, b in enumerate(corners):
            if j > i and np.sum(a != b) == 1:
                ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]],
                        color=color, ls=style, lw=0.8, alpha=0.58,
                        label=label if first else None)
                first = False


def compute_axis_limits(reference, smpc, mpc, hard_bound, tightened_bound,
                        padding=0.05):
    """Return per-axis ``(lo, hi)`` limits containing every finite trajectory
    point of ``reference``/``smpc``/``mpc`` with ``padding`` extra room on each
    side (relative to the per-axis full span).

    The hard box (``+-hard_bound``) and the tightened box (axis-asymmetric)
    are still drawn by the caller; this helper only prevents axis cropping
    when the controllers blow past the hard bound by orders of magnitude
    (the worst observed peaks in this scenario exceed 700 in standardized
    coordinates, while ``hard_bound ~ 3.91``).

    The returned limits are symmetric around zero per axis (``[-L, L]``) so
    that the box wireframes drawn on top of the trajectories remain visually
    centered. The minimum limit is at least the per-axis tightened bound
    value so that even on unrealistically well-behaved data the box stays
    visible.

    Pure function: no matplotlib, no side effects, no RNG.
    """
    arrays = [
        np.asarray(reference, dtype=float),
        np.asarray(smpc["S"], dtype=float),
        np.asarray(mpc["S"], dtype=float),
    ]
    finite_mask = np.ones(arrays[0].shape[0], dtype=bool)
    for arr in arrays:
        finite_mask &= np.all(np.isfinite(arr), axis=1)
    if not finite_mask.any():
        finite_mask = np.ones(arrays[0].shape[0], dtype=bool)

    finite = [arr[finite_mask] for arr in arrays]
    abs_max = max(
        float(np.abs(arr).max()) for arr in finite if arr.size > 0
    )
    box_min = max(
        float(np.abs(np.asarray(hard_bound)).ravel()[0]),
        float(np.max(np.abs(tightened_bound))),
        1.0,
    )
    lower = max(abs_max * (1.0 + float(padding)), box_min)
    # Per-axis symmetric limits around zero: row index = axis,
    # columns = [lo, hi]. Returned as a ``(3, 2)`` array so the caller
    # can index ``limits[axis, 0]`` / ``limits[axis, 1]`` for ax.set.
    return np.array(
        [[-lower, lower]] * 3, dtype=float,
    )


def _violation_metrics(trajectory, hard_bound):
    """Return ``(violation_steps, minimum_hard_margin)`` from a closed-loop
    trajectory. Mirrors ``BF.simulate`` exactly so the figure caption and
    the on-disk JSON stay aligned without re-running the controller."""
    task = np.asarray(trajectory, dtype=float)
    hard = float(np.asarray(hard_bound).ravel()[0])
    violations = int(np.sum(np.any(np.abs(task) > hard, axis=1)))
    margin = float(np.min(hard - np.abs(task)))
    return violations, margin


def plot_spiral(reference, smpc, mpc, hard_bound, tightened_bound, path,
                sigma_w=SIGMA_W, seed=WIND_SEED, total_steps=None):
    fig = plt.figure(figsize=(10.5, 8.6))
    ax = fig.add_subplot(111, projection="3d")
    limits = compute_axis_limits(
        reference, smpc, mpc, hard_bound, tightened_bound,
    )
    _draw_box(ax, np.full(3, hard_bound), "#1a1a1a", "-.", "hard bound")
    _draw_box(ax, tightened_bound, "#7b2cbf", "--",
              "terminal tightened mean bound")
    ax.plot(reference[:, 0], reference[:, 1], reference[:, 2],
            color="#1a1a1a", ls="--", lw=1.35, label="spiral reference")
    ax.plot(mpc["S"][:, 0], mpc["S"][:, 1], mpc["S"][:, 2],
            color="#c44536", ls=":", lw=1.75, label="deterministic MPC")
    ax.plot(smpc["S"][:, 0], smpc["S"][:, 1], smpc["S"][:, 2],
            color="#087e8b", lw=1.15, label="SMPC")
    ax.scatter(*smpc["S"][0], color="#087e8b", marker="o", s=32,
               label="start at bottom")
    ax.scatter(*smpc["S"][-1], color="#087e8b", marker="x", s=42,
               label="end at top")
    ax.set(xlim=(limits[0, 0], limits[0, 1]),
           ylim=(limits[1, 0], limits[1, 1]),
           zlim=(limits[2, 0], limits[2, 1]),
           xlabel="standardized x",
           ylabel="standardized y", zlabel="standardized z")
    # Box-aspect 1:1:1 keeps the per-axis ranges honest: with these limits
    # the visible cube is ~[1400]^3, so the +-3.9 hard box looks like a
    # tiny dot at the center. That is the visual cue the user wants.
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=23, azim=-58)
    smpc_steps, smpc_margin = _violation_metrics(smpc["S"], hard_bound)
    mpc_steps, mpc_margin = _violation_metrics(mpc["S"], hard_bound)
    if total_steps is None:
        total_steps = int(np.asarray(smpc["S"]).shape[0])
    ax.set_title(
        f"copyBF: 120 s bottom-to-top tightened-bound spiral with Gaussian "
        f"process wind (sigma_w={sigma_w:.4f}, seed={seed})"
    )
    fig.text(
        0.5, 0.005,
        f"SMPC: {smpc_steps}/{total_steps} steps violate hard bound, "
        f"min hard margin {smpc_margin:+.2f}  |  "
        f"MPC: {mpc_steps}/{total_steps} steps violate hard bound, "
        f"min hard margin {mpc_margin:+.2f}",
        ha="center", va="bottom", fontsize=8, color="#333333",
    )
    ax.legend(loc="upper left", fontsize=8)
    fig.subplots_adjust(bottom=0.13)
    fig.savefig(path, dpi=190)
    plt.close(fig)


def redraw_png_from_mat(mat_path, png_path=None, sigma_w=SIGMA_W,
                        wind_seed=WIND_SEED):
    """Reload the closed-loop result from ``mat_path`` and rewrite the PNG.

    No simulation. Loads ``reference``, ``smpc_trajectory``,
    ``deterministic_mpc_trajectory``, ``hard_bound``, and ``tightened_bound``
    from the MAT artifact, then calls :func:`plot_spiral` with the original
    trajectory arrays. Used to refresh the figure after the axis-limits bug
    fix without re-running the 12000-step closed loop.
    """
    mat_path = Path(mat_path)
    if png_path is None:
        png_path = mat_path.with_suffix(".png")
    payload = scipy.io.loadmat(str(mat_path))
    smpc_arr = np.asarray(payload["smpc_trajectory"], dtype=float)
    mpc_arr = np.asarray(payload["deterministic_mpc_trajectory"], dtype=float)
    reference = np.asarray(payload["reference"], dtype=float)
    hard_bound = float(np.asarray(payload["hard_bound"]).ravel()[0])
    tightened_bound = np.asarray(
        payload["tightened_bound"], dtype=float
    ).ravel()
    total_steps = int(smpc_arr.shape[0])
    # Mirror the metric dict the simulator emits so plot_spiral's caption
    # reads the same numbers as JSON/MAT.
    smpc = {"S": smpc_arr, "U": np.asarray(
        payload["smpc_control"], dtype=float
    )}
    mpc = {"S": mpc_arr, "U": np.asarray(
        payload["deterministic_mpc_control"], dtype=float
    )}
    plot_spiral(reference, smpc, mpc, hard_bound, tightened_bound, png_path,
                sigma_w=sigma_w, seed=wind_seed, total_steps=total_steps)
    return png_path


def run(seed=RESIDUAL_SEED, steps=12000, sigma_w=SIGMA_W, wind_seed=WIND_SEED):
    model, E, hard_bound, _ = BF.BE.load_model()
    controller = BF.InstrumentedZAwareSMPC(model, E, task_bound=hard_bound)
    tightening = BF.terminal_tightening(controller, E)
    tightened = hard_bound - tightening
    reference = tightened_boundary_spiral(
        tightened, steps, rx=2.7, ry=2.8, z_lo=-3.75, z_hi=3.75,
    )
    ell = model["A"].shape[0] // 2

    # External wind channel: per-step independent Gaussian process noise in
    # latent coordinates. Drawn from its own RNG seed (WIND_SEED != seed) so
    # this stream shares no underlying base sequence with the residual noise.
    wind = process_wind_noise(steps, ell, sigma_w=sigma_w, seed=wind_seed)

    # Identified residual innovation N(0, Sigma_eps) from the Pelican
    # dataset. Treated strictly as the latent innovation estimated from
    # data; it is not the SMPC measurement noise and is independent of
    # the wind channel by construction (separate RNG seed).
    noise = residual_innovation_sequence(model, ell, steps, seed=seed)

    # Both controllers see the exact same realizations below.
    smpc = BF.simulate(model, E, hard_bound, reference, noise, wind, 0.05,
                       period_steps=1200, initial_task=reference[0])
    mpc = BF.simulate(model, E, hard_bound, reference, noise, wind, None,
                      period_steps=1200, initial_task=reference[0])

    result = {
        "scope": (
            "continuous 120 s frozen-model-in-the-loop tightened-bound "
            "spiral with Gaussian process wind; not real flight"
        ),
        "seed": seed,
        "wind_seed": int(wind_seed),
        "sigma_w": float(sigma_w),
        "wind_channel": (
            "per-step independent Gaussian process noise w_k ~ "
            f"N(0, {sigma_w}^2 I_ell) injected in latent coordinates"
        ),
        "residual_channel": (
            "latent innovation N(0, Sigma_eps[:ell,:ell]) estimated from "
            "the Pelican dataset; not measurement noise; SMPC and MPC see "
            "the same realization"
        ),
        "steps": steps, "duration_seconds": steps / 100.0,
        "state_resets": 0, "turns": 5,
        "hard_bound": hard_bound,
        "terminal_tightened_bound": tightened.tolist(),
        "reference_min": reference.min(axis=0).tolist(),
        "reference_max": reference.max(axis=0).tolist(),
        "tightened_bound_at_terminal": tightened.tolist(),
        "smpc": BF._serializable_metrics(smpc),
        "deterministic_mpc": BF._serializable_metrics(mpc),
    }
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    (out / "copyBF_continuous_3d_tightened_spiral.json").write_text(
        json.dumps(result, indent=2), encoding="ascii"
    )
    np.savez(
        out / "copyBF_continuous_3d_tightened_spiral.npz",
        reference=reference,
        residual_noise=noise,
        process_wind=wind,
        smpc_trajectory=smpc["S"],
        smpc_control=smpc["U"],
        deterministic_mpc_trajectory=mpc["S"],
        deterministic_mpc_control=mpc["U"],
    )

    mat_payload = build_mat_payload(
        reference, smpc, mpc, noise, wind, hard_bound, tightened,
        sigma_w, seed, wind_seed,
    )
    save_mat_results(
        out / "copyBF_continuous_3d_tightened_spiral.mat", mat_payload,
    )

    plot_spiral(reference, smpc, mpc, hard_bound, tightened,
                out / "copyBF_continuous_3d_tightened_spiral.png",
                sigma_w=sigma_w, seed=wind_seed, total_steps=steps)
    print(json.dumps({
        "smpc_active": smpc["active_qp_steps"],
        "smpc_fallback": smpc["fallback_count"],
        "mpc_fallback": mpc["fallback_count"],
        "smpc_violations": smpc["violation_steps"],
        "mpc_violations": mpc["violation_steps"],
        "smpc_input_saturation_steps": smpc["input_saturation_steps"],
        "mpc_input_saturation_steps": mpc["input_saturation_steps"],
        "sigma_w": sigma_w,
        "residual_seed": seed,
        "wind_seed": wind_seed,
    }))
    return result


def _main():
    parser = argparse.ArgumentParser(
        description=(
            "Run the copyBF tightened-bound spiral closed-loop, or reload "
            "an existing MAT artifact and redraw the 3D PNG without "
            "re-running the 12000-step simulation."
        ),
    )
    parser.add_argument(
        "--redraw-png", action="store_true",
        help="Load the existing MAT and rewrite the PNG only.",
    )
    parser.add_argument(
        "--mat-path", type=Path, default=None,
        help="Override the MAT artifact path used with --redraw-png.",
    )
    parser.add_argument(
        "--png-output", type=Path, default=None,
        help="Override the PNG output path used with --redraw-png.",
    )
    args = parser.parse_args()
    if args.redraw_png:
        mat_path = (
            args.mat_path
            if args.mat_path is not None
            else HERE / "results" / "copyBF_continuous_3d_tightened_spiral.mat"
        )
        png_path = redraw_png_from_mat(
            mat_path, png_path=args.png_output,
            sigma_w=SIGMA_W, wind_seed=WIND_SEED,
        )
        print(f"redrew PNG: {png_path}")
        return png_path
    return run()


if __name__ == "__main__":
    _main()
