"""Continuous 120 s 3D boundary-tour experiment built on copyBF."""

import importlib.util
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
BF_PATH = HERE / "run_long_3d_boundary_disturbance.py"
_spec = importlib.util.spec_from_file_location("copybf_base", BF_PATH)
BF = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(BF)


def continuous_boundary_tour(steps=12000, targets=(3.60, 3.55, 3.75)):
    """Closed x->y->z->x tour with C1-smooth anchor transitions."""
    anchors = np.array([
        [targets[0], 0.0, 0.0],
        [0.0, targets[1], 0.0],
        [0.0, 0.0, targets[2]],
        [targets[0], 0.0, 0.0],
    ])
    t = np.linspace(0.0, 3.0, int(steps), endpoint=True)
    segment = np.minimum(np.floor(t).astype(int), 2)
    local = t - segment
    blend = local * local * (3.0 - 2.0 * local)
    reference = (
        (1.0 - blend[:, None]) * anchors[segment]
        + blend[:, None] * anchors[segment + 1]
    )
    return reference


def directional_zero_impulse_gusts(task_map, reference, period_steps=1200,
                                    width=20, amplitude=0.0002):
    task_gusts = np.zeros_like(reference)
    for cycle_start in range(0, len(reference), int(period_steps)):
        start = cycle_start + int(0.5 * period_steps)
        stop = min(start + 2 * width, len(reference), cycle_start + period_steps)
        if stop - start < 4:
            continue
        half = (stop - start) // 2
        window = np.sin(np.linspace(0.0, np.pi, half, endpoint=True))
        window /= np.max(window)
        pulse = np.concatenate((window, -window))
        direction = reference[start] / max(np.linalg.norm(reference[start]), 1e-12)
        task_gusts[start:start + len(pulse)] += amplitude * pulse[:, None] * direction
    latent = (np.linalg.pinv(task_map) @ task_gusts.T).T
    return latent, task_gusts


def plot_trajectory(reference, smpc, mpc, hard_bound, tightened, path):
    fig = plt.figure(figsize=(10.5, 8.4))
    ax = fig.add_subplot(111, projection="3d")
    # Draw the positive-octant edges relevant to this boundary tour.
    for limits, color, style, label in [
        (np.full(3, hard_bound), "#1a1a1a", "-.", "hard bound"),
        (np.asarray(tightened), "#7b2cbf", "--", "terminal tightened mean bound"),
    ]:
        corners = np.array([[sx * limits[0], sy * limits[1], sz * limits[2]]
                            for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)])
        first = True
        for i, a in enumerate(corners):
            for j, b in enumerate(corners):
                if j > i and np.sum(a != b) == 1:
                    ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]],
                            color=color, ls=style, lw=0.8, alpha=0.55,
                            label=label if first else None)
                    first = False
    ax.plot(reference[:, 0], reference[:, 1], reference[:, 2],
            color="#1a1a1a", ls="--", lw=1.5, label="reference")
    ax.plot(mpc["S"][:, 0], mpc["S"][:, 1], mpc["S"][:, 2],
            color="#c44536", ls=":", lw=1.8, label="deterministic MPC")
    ax.plot(smpc["S"][:, 0], smpc["S"][:, 1], smpc["S"][:, 2],
            color="#087e8b", lw=1.2, label="SMPC")
    ax.scatter(*smpc["S"][0], color="#087e8b", marker="o", s=35, label="start/end")
    ax.set(xlim=(-hard_bound, hard_bound), ylim=(-hard_bound, hard_bound),
           zlim=(-hard_bound, hard_bound), xlabel="standardized x",
           ylabel="standardized y", zlabel="standardized z")
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=25, azim=-58)
    ax.set_title("copyBF: continuous 120 s 3D boundary tour with external gusts")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=190)
    plt.close(fig)


def run(seed=7, steps=12000):
    model, E, hard_bound, _ = BF.BE.load_model()
    controller = BF.InstrumentedZAwareSMPC(model, E, task_bound=hard_bound)
    tightening = BF.terminal_tightening(controller, E)
    tightened = hard_bound - tightening
    reference = continuous_boundary_tour(steps)
    ell = model["A"].shape[0] // 2
    task_map = E.T @ model["P"][:, :ell]
    gusts, task_gusts = directional_zero_impulse_gusts(task_map, reference)
    noise = np.random.default_rng(seed).multivariate_normal(
        np.zeros(ell), model["Sigma_eps"][:ell, :ell], size=steps
    )
    initial_task = np.array([tightened[0] - 0.05, 0.0, 0.0])
    smpc = BF.simulate(model, E, hard_bound, reference, noise, gusts, 0.05,
                       period_steps=1200, initial_task=initial_task)
    mpc = BF.simulate(model, E, hard_bound, reference, noise, gusts, None,
                      period_steps=1200, initial_task=initial_task)
    result = {
        "scope": "continuous 120 s frozen-model-in-the-loop 3D tour; not real flight",
        "seed": seed, "steps": steps, "duration_seconds": steps / 100.0,
        "state_resets": 0, "hard_bound": hard_bound,
        "terminal_tightened_bound": tightened.tolist(),
        "smpc": BF._serializable_metrics(smpc),
        "deterministic_mpc": BF._serializable_metrics(mpc),
    }
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    (out / "copyBF_continuous_3d_boundary_tour.json").write_text(
        json.dumps(result, indent=2), encoding="ascii"
    )
    np.savez(out / "copyBF_continuous_3d_boundary_tour.npz",
             reference=reference, smpc=smpc["S"], deterministic_mpc=mpc["S"],
             external_task_gusts=task_gusts)
    plot_trajectory(reference, smpc, mpc, hard_bound, tightened,
                    out / "copyBF_continuous_3d_boundary_tour.png")
    print(json.dumps({
        "smpc_active": smpc["active_qp_steps"],
        "smpc_fallback": smpc["fallback_count"],
        "mpc_fallback": mpc["fallback_count"],
        "smpc_violations": smpc["violation_steps"],
        "mpc_violations": mpc["violation_steps"],
    }))
    return result


if __name__ == "__main__":
    run()
