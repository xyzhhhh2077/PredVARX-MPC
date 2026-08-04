"""Fixed-anchor CRTE identification model in a noisy closed-loop simulation.

The model is identified only from Waterloo Pelican flights 1-36. Held-out
flights provide initial conditions and reference trajectories; their measured
next states are never injected into the simulated plant.
"""

import argparse
import json
from pathlib import Path

import cvxpy as cp
import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
DATA_MAT = HERE.parent / "copyAY_pelican_soft_preference" / "data" / "copyAY_pelican_dataset.mat"
MODEL_MAT = HERE / "results" / "copyBB_pelican_fixed_anchor_data.mat"


def fixed_position_anchor(p=10):
    """Prior controlled axes: standardized Pelican x and y position."""
    if p < 2:
        raise ValueError("At least two outputs are required")
    E = np.zeros((p, 2))
    E[0, 0] = 1.0
    E[1, 1] = 1.0
    return E


def task_cost_matrix(E, weight=80.0):
    """Original copyAU task-only state penalty."""
    return float(weight) * E @ E.T


def mstep_lift(A, B, d):
    Ad = np.linalg.matrix_power(A, d)
    Bd = np.zeros_like(B)
    Ai = np.eye(A.shape[0])
    for _ in range(d):
        Bd += Ai @ B
        Ai = Ai @ A
    return Ad, Bd


def lift_noise_covariance(A, Sigma_eps, d):
    """Covariance accumulated over d physical plant steps."""
    Sigma_d = np.zeros_like(Sigma_eps)
    Ai = np.eye(A.shape[0])
    for _ in range(d):
        Sigma_d += Ai @ Sigma_eps @ Ai.T
        Ai = Ai @ A
    return (Sigma_d + Sigma_d.T) / 2.0


def plant_step_with_noise(model, z, u, eps):
    return (
        model["A"] @ z
        + model["B"] @ (np.asarray(u).reshape(-1, 1) - model["u_mean"])
        + eps
    )


def load_model_and_data():
    with h5py.File(MODEL_MAT, "r") as f:
        m = f["model"]
        model = {
            key: np.asarray(m[key]).T
            for key in ("A", "B", "P", "R", "Sigma_eps", "y_mean", "u_mean")
        }
        E = np.asarray(f["Eanchor"]).T
        y_offset = np.asarray(f["y_offset"]).ravel()
        y_scale = np.asarray(f["y_scale"]).ravel()
        u_offset = np.asarray(f["u_offset"]).ravel()
        u_scale = np.asarray(f["u_scale"]).ravel()
        h = float(np.asarray(f["task_bound"]).ravel()[0])
    data = sio.loadmat(DATA_MAT)
    scales = {
        "y_offset": y_offset,
        "y_scale": y_scale,
        "u_offset": u_offset,
        "u_scale": u_scale,
    }
    return model, E, scales, h, data


class FixedAnchorSMPC:
    def __init__(self, model, E, d=10, horizon=18, q_weight=80.0, ru=0.18,
                 task_bound=3.5, alpha_joint=0.05):
        self.model = model
        self.E = E
        self.d = d
        self.N = horizon
        self.nu = model["B"].shape[1]
        self.ny = model["P"].shape[0]
        self.Ad, self.Bd = mstep_lift(model["A"], model["B"], d)
        self.Sigma_d = lift_noise_covariance(model["A"], model["Sigma_eps"], d)
        self.Q = task_cost_matrix(E, q_weight)
        self.Ru = float(ru) * np.eye(self.nu)
        self.u_min = np.full(self.nu, -6.0)
        self.u_max = np.full(self.nu, 6.0)
        self.task_bound = float(task_bound)
        self.alpha_joint = float(alpha_joint)
        self.y_mean = model["y_mean"].reshape(-1, 1)
        self.u_mean = model["u_mean"].reshape(-1, 1)

        self.M = []
        self.G = []
        for j in range(1, self.N + 1):
            self.M.append(model["P"] @ np.linalg.matrix_power(self.Ad, j))
            Gj = np.zeros((self.ny, self.N * self.nu))
            for i in range(j):
                block = model["P"] @ np.linalg.matrix_power(self.Ad, j - 1 - i) @ self.Bd
                Gj[:, i * self.nu:(i + 1) * self.nu] = block
            self.G.append(Gj)

        self.Sigma_z = []
        Sz = np.zeros_like(self.Ad)
        for _ in range(self.N):
            Sz = self.Ad @ Sz @ self.Ad.T + self.Sigma_d
            self.Sigma_z.append((Sz + Sz.T) / 2.0)

    def step(self, y_std, s_reference):
        z = self.model["R"].T @ (np.asarray(y_std).reshape(-1, 1) - self.y_mean)
        U0 = np.tile(self.u_mean, (self.N, 1))
        U = cp.Variable((self.N * self.nu, 1))
        r_full = self.E @ np.asarray(s_reference).reshape(-1, 1)
        cost = cp.quad_form(U - U0, np.kron(np.eye(self.N), self.Ru))
        constraints = [np.tile(self.u_min, self.N).reshape(-1, 1) <= U,
                       U <= np.tile(self.u_max, self.N).reshape(-1, 1)]
        risk_each = self.alpha_joint / (2.0 * self.E.shape[1] * self.N)
        z_quantile = norm.ppf(1.0 - risk_each)
        for j in range(self.N):
            mu = self.M[j] @ z + self.y_mean + self.G[j] @ (U - U0)
            cost += cp.quad_form(mu - r_full, self.Q)
            Sigma_y = self.model["P"] @ self.Sigma_z[j] @ self.model["P"].T
            for axis in range(self.E.shape[1]):
                hq = self.E[:, axis].reshape(1, -1)
                variance = max(float((hq @ Sigma_y @ hq.T)[0, 0]), 1e-12)
                tight = z_quantile * np.sqrt(variance)
                constraints += [hq @ mu <= self.task_bound - tight,
                                hq @ mu >= -self.task_bound + tight]
        problem = cp.Problem(cp.Minimize(cost), constraints)
        problem.solve(solver=cp.OSQP, verbose=False, warm_start=True)
        if problem.status not in ("optimal", "optimal_inaccurate"):
            return None, problem.status
        return np.clip(U.value[:self.nu].ravel(), self.u_min, self.u_max), problem.status


def standardized_segment(data, scales, segid, steps):
    seg = data["segment_id"].ravel()
    idx = np.flatnonzero(seg == segid)[:steps]
    if len(idx) < steps:
        raise ValueError(f"segment {segid} has only {len(idx)} samples")
    y_std = (data["y"][:, idx] - scales["y_offset"][:, None]) / scales["y_scale"][:, None]
    u_std = (data["u"][:, idx] - scales["u_offset"][:, None]) / scales["u_scale"][:, None]
    return y_std, u_std


def simulate_segment(model, E, scales, data, task_bound, segid, steps, seed,
                     d, horizon, q_weight, ru):
    y_real, u_expert = standardized_segment(data, scales, segid, steps)
    S_ref = (E.T @ y_real).T
    z0 = model["R"].T @ (y_real[:, [0]] - model["y_mean"])
    z_ctrl = z0.copy()
    z_open = z0.copy()
    controller = FixedAnchorSMPC(model, E, d=d, horizon=horizon,
                                 q_weight=q_weight, ru=ru,
                                 task_bound=task_bound)
    rng = np.random.default_rng(seed)
    S_ctrl = np.zeros_like(S_ref)
    S_open = np.zeros_like(S_ref)
    U_ctrl = np.zeros((steps, model["B"].shape[1]))
    u_ctrl = model["u_mean"].ravel()
    u_open = model["u_mean"].ravel()
    fallback = 0
    for k in range(steps):
        y_ctrl = model["y_mean"] + model["P"] @ z_ctrl
        y_open = model["y_mean"] + model["P"] @ z_open
        S_ctrl[k] = (E.T @ y_ctrl).ravel()
        S_open[k] = (E.T @ y_open).ravel()
        if k % d == 0:
            u_new, status = controller.step(y_ctrl.ravel(), S_ref[k])
            if u_new is None:
                fallback += 1
                u_ctrl = model["u_mean"].ravel()
            else:
                u_ctrl = u_new
        U_ctrl[k] = u_ctrl
        eps = rng.multivariate_normal(
            np.zeros(model["A"].shape[0]), model["Sigma_eps"]
        ).reshape(-1, 1)
        z_ctrl = plant_step_with_noise(model, z_ctrl, u_ctrl, eps)
        z_open = plant_step_with_noise(model, z_open, u_open, eps)

    rmse = np.sqrt(np.mean((S_ctrl - S_ref) ** 2, axis=0))
    rmse_open = np.sqrt(np.mean((S_open - S_ref) ** 2, axis=0))
    metrics = {
        "rmse_smpc": rmse.tolist(),
        "rmse_openloop": rmse_open.tolist(),
        "improvement_pct": (100.0 * (rmse_open - rmse) / rmse_open).tolist(),
        "reference_std": S_ref.std(axis=0).tolist(),
        "nrmse_smpc": (rmse / np.maximum(S_ref.std(axis=0), 1e-12)).tolist(),
        "input_std_smpc": U_ctrl.std(axis=0).tolist(),
        "input_std_expert": u_expert.T.std(axis=0).tolist(),
        "saturation_pct": float(100.0 * np.mean(np.max(np.abs(U_ctrl), axis=1) >= 5.999)),
        "violation_steps": int(np.sum(np.max(np.abs(S_ctrl), axis=1) > task_bound)),
        "fallback_count": int(fallback),
    }
    traces = {"S_ref": S_ref, "S_smpc": S_ctrl, "S_open": S_open, "U_smpc": U_ctrl}
    return metrics, traces


def run_experiment(segs, steps=1200, seed=7, d=10, horizon=18,
                   q_weight=80.0, ru=0.18, prefix="copyBB_fixed_anchor_closed_loop"):
    model, E, scales, task_bound, data = load_model_and_data()
    per_seg = {}
    all_traces = {}
    for segid in segs:
        metrics, traces = simulate_segment(
            model, E, scales, data, task_bound, segid, steps, seed,
            d, horizon, q_weight, ru,
        )
        per_seg[str(segid)] = metrics
        all_traces[str(segid)] = traces
        print(f"seg{segid}: RMSE={np.round(metrics['rmse_smpc'], 3)} "
              f"open={np.round(metrics['rmse_openloop'], 3)} "
              f"improve={np.round(metrics['improvement_pct'], 1)}% "
              f"sat={metrics['saturation_pct']:.1f}% fallback={metrics['fallback_count']}")

    out = {
        "scope": "real-data-trained frozen model; noisy model-in-the-loop control simulation",
        "training_flights": list(range(1, 37)),
        "held_out_reference_flights": [int(x) for x in segs],
        "fixed_anchor": "standardized x/y position channels",
        "d": d, "horizon": horizon, "q_weight": q_weight, "ru": ru,
        "task_bound": task_bound, "steps": steps, "seed": seed,
        "per_seg": per_seg,
    }
    results = HERE / "results"
    results.mkdir(exist_ok=True)
    (results / f"{prefix}.json").write_text(json.dumps(out, indent=2), encoding="ascii")
    np.savez(results / f"{prefix}.npz", traces=all_traces, metrics=out)
    controller = FixedAnchorSMPC(
        model, E, d=d, horizon=horizon, q_weight=q_weight, ru=ru,
        task_bound=task_bound,
    )
    terminal_tightening = np.empty(E.shape[1])
    risk_each = controller.alpha_joint / (2.0 * E.shape[1] * horizon)
    z_quantile = norm.ppf(1.0 - risk_each)
    Sigma_y_terminal = model["P"] @ controller.Sigma_z[-1] @ model["P"].T
    for axis in range(E.shape[1]):
        hq = E[:, axis].reshape(1, -1)
        variance = max(float((hq @ Sigma_y_terminal @ hq.T)[0, 0]), 1e-12)
        terminal_tightening[axis] = z_quantile * np.sqrt(variance)
    plot_results(
        all_traces, segs, results / f"{prefix}.png", task_bound,
        task_bound - terminal_tightening,
    )
    return out


def plot_results(traces, segs, path, task_bound, terminal_mean_bound):
    fig, axes = plt.subplots(len(segs), 2, figsize=(13, 2.5 * len(segs)), sharex=True)
    axes = np.atleast_2d(axes)
    for row, segid in enumerate(segs):
        t = np.arange(len(traces[str(segid)]["S_ref"])) / 100.0
        for axis in range(2):
            ax = axes[row, axis]
            ax.plot(t, traces[str(segid)]["S_ref"][:, axis], "k", lw=1.2, label="held-out reference")
            ax.plot(t, traces[str(segid)]["S_smpc"][:, axis], color="#087e8b", lw=1.0, label="fixed-anchor SMPC")
            ax.plot(t, traces[str(segid)]["S_open"][:, axis], color="#c44536", lw=0.8, alpha=0.75, label="open loop")
            hard_label = "hard bound" if row == 0 and axis == 0 else None
            mean_label = "1.8 s tightened mean bound" if row == 0 and axis == 0 else None
            ax.axhline(task_bound, color="#333333", ls="-.", lw=0.9, label=hard_label)
            ax.axhline(-task_bound, color="#333333", ls="-.", lw=0.9)
            ax.axhline(terminal_mean_bound[axis], color="#7b2cbf", ls="--", lw=0.9, label=mean_label)
            ax.axhline(-terminal_mean_bound[axis], color="#7b2cbf", ls="--", lw=0.9)
            ax.grid(True, alpha=0.25)
            ax.set_ylabel(f"flight {segid}\nstd. {'x' if axis == 0 else 'y'}")
            if row == 0:
                ax.set_title(f"Fixed {'x' if axis == 0 else 'y'} position task")
            if row == len(segs) - 1:
                ax.set_xlabel("simulation time (s)")
    axes[0, 0].legend(loc="best", ncol=2, fontsize=7)
    fig.suptitle("copyBB: real-data-trained model, noisy model-in-the-loop SMPC", fontsize=13)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--segs", nargs="+", type=int, default=[37, 38, 39, 40, 41])
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--d", type=int, default=10)
    parser.add_argument("--horizon", type=int, default=18)
    parser.add_argument("--q-weight", type=float, default=80.0)
    parser.add_argument("--ru", type=float, default=0.18)
    parser.add_argument("--prefix", default="copyBB_fixed_anchor_closed_loop")
    args = parser.parse_args()
    run_experiment(args.segs, args.steps, args.seed, args.d, args.horizon,
                   args.q_weight, args.ru, args.prefix)


if __name__ == "__main__":
    main()
