"""Three-dimensional fixed-anchor Pelican SMPC trajectory experiment."""

import argparse
import json
from pathlib import Path

import cvxpy as cp
import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
MODEL_MAT = HERE / "results" / "copyBD_pelican_3d_fixed_anchor_data.mat"


def fixed_position_anchor_3d(p=10):
    if p < 3:
        raise ValueError("At least three outputs are required")
    E = np.zeros((p, 3))
    E[:3, :3] = np.eye(3)
    return E


def task_cost_matrix(E, weight=80.0):
    return float(weight) * E @ E.T


def per_face_risk(alpha_joint, task_axes, horizon):
    return float(alpha_joint) / (2.0 * int(task_axes) * int(horizon))


def moving_reference_3d(steps, hard_bound, x_amplitude=0.25,
                        lateral_amplitude=0.55, vertical_amplitude=0.32,
                        cycles=1.0):
    """Trackable closed 3D curve with nonzero motion in every task axis."""
    del hard_bound
    phase = np.linspace(0.0, 2.0 * np.pi * cycles, int(steps), endpoint=True)
    x = float(x_amplitude) * np.cos(phase)
    y = float(lateral_amplitude) * np.sin(phase)
    z = float(vertical_amplitude) * np.sin(2.0 * phase)
    return np.column_stack([x, y, z])


def future_reference_horizon(reference, k, d, horizon):
    """Reference samples aligned with the N lifted predicted states."""
    reference = np.asarray(reference)
    indices = np.minimum(
        int(k) + int(d) * np.arange(1, int(horizon) + 1),
        len(reference) - 1,
    )
    return reference[indices]


def mstep_lift(A, B, d):
    Ad = np.linalg.matrix_power(A, d)
    Bd = np.zeros_like(B)
    Ai = np.eye(A.shape[0])
    for _ in range(d):
        Bd += Ai @ B
        Ai = Ai @ A
    return Ad, Bd


def lift_noise_covariance(A, Sigma_eps, d):
    Sigma_d = np.zeros_like(Sigma_eps)
    Ai = np.eye(A.shape[0])
    for _ in range(d):
        Sigma_d += Ai @ Sigma_eps @ Ai.T
        Ai = Ai @ A
    return (Sigma_d + Sigma_d.T) / 2.0


def load_model():
    with h5py.File(MODEL_MAT, "r") as f:
        m = f["model"]
        model = {
            key: np.asarray(m[key]).T
            for key in ("A", "B", "P", "R", "Sigma_eps", "y_mean", "u_mean")
        }
        E = np.asarray(f["Eanchor"]).T
        hard_bound = float(np.asarray(f["task_bound"]).ravel()[0])
        scales = {
            "y_offset": np.asarray(f["y_offset"]).ravel(),
            "y_scale": np.asarray(f["y_scale"]).ravel(),
        }
    return model, E, hard_bound, scales


class FixedAnchorSMPC3D:
    def __init__(self, model, E, d=10, horizon=18, q_weight=80.0, ru=0.18,
                 task_bound=3.5, alpha_joint=0.05):
        self.model = model
        self.E = E
        self.d = int(d)
        self.N = int(horizon)
        self.nu = model["B"].shape[1]
        self.ny = model["P"].shape[0]
        self.Ad, self.Bd = mstep_lift(model["A"], model["B"], self.d)
        self.Sigma_d = lift_noise_covariance(model["A"], model["Sigma_eps"], self.d)
        self.Q = task_cost_matrix(E, q_weight)
        self.Ru = float(ru) * np.eye(self.nu)
        self.u_min = np.full(self.nu, -6.0)
        self.u_max = np.full(self.nu, 6.0)
        self.task_bound = float(task_bound)
        self.alpha_joint = alpha_joint
        self.y_mean = model["y_mean"].reshape(-1, 1)
        self.u_mean = model["u_mean"].reshape(-1, 1)

        self.M, self.G = [], []
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

    def step(self, y_std, s_reference_horizon):
        z = self.model["R"].T @ (np.asarray(y_std).reshape(-1, 1) - self.y_mean)
        U0 = np.tile(self.u_mean, (self.N, 1))
        U = cp.Variable((self.N * self.nu, 1))
        references = np.asarray(s_reference_horizon)
        if references.shape != (self.N, self.E.shape[1]):
            raise ValueError("Reference horizon must have shape (N, task_axes)")
        cost = cp.quad_form(U - U0, np.kron(np.eye(self.N), self.Ru))
        constraints = [np.tile(self.u_min, self.N).reshape(-1, 1) <= U,
                       U <= np.tile(self.u_max, self.N).reshape(-1, 1)]
        z_quantile = None
        if self.alpha_joint is not None:
            risk_each = per_face_risk(self.alpha_joint, self.E.shape[1], self.N)
            z_quantile = norm.ppf(1.0 - risk_each)
        for j in range(self.N):
            mu = self.M[j] @ z + self.y_mean + self.G[j] @ (U - U0)
            r_full = self.E @ references[j].reshape(-1, 1)
            cost += cp.quad_form(mu - r_full, self.Q)
            Sigma_y = self.model["P"] @ self.Sigma_z[j] @ self.model["P"].T
            for axis in range(self.E.shape[1]):
                hq = self.E[:, axis].reshape(1, -1)
                tight = 0.0
                if z_quantile is not None:
                    variance = max(float((hq @ Sigma_y @ hq.T)[0, 0]), 1e-12)
                    tight = z_quantile * np.sqrt(variance)
                constraints += [hq @ mu <= self.task_bound - tight,
                                hq @ mu >= -self.task_bound + tight]
        problem = cp.Problem(cp.Minimize(cost), constraints)
        problem.solve(solver=cp.OSQP, verbose=False, warm_start=True)
        if problem.status not in ("optimal", "optimal_inaccurate"):
            return None, problem.status
        return np.clip(U.value[:self.nu].ravel(), self.u_min, self.u_max), problem.status


def simulate(model, E, hard_bound, reference, disturbances, alpha_joint,
             d=10, horizon=18, q_weight=80.0, ru=0.18):
    controller = FixedAnchorSMPC3D(
        model, E, d=d, horizon=horizon, q_weight=q_weight, ru=ru,
        task_bound=hard_bound, alpha_joint=alpha_joint,
    )
    steps = len(reference)
    y0 = model["y_mean"] + E @ reference[0].reshape(-1, 1)
    z = model["R"].T @ (y0 - model["y_mean"])
    S = np.zeros((steps, 3))
    U = np.zeros((steps, model["B"].shape[1]))
    u = model["u_mean"].ravel()
    fallback = 0
    qp_count = 0
    for k in range(steps):
        y = model["y_mean"] + model["P"] @ z
        S[k] = (E.T @ y).ravel()
        if k % d == 0:
            reference_horizon = future_reference_horizon(reference, k, d, horizon)
            u_new, _ = controller.step(y.ravel(), reference_horizon)
            qp_count += 1
            if u_new is None:
                fallback += 1
                u = model["u_mean"].ravel()
            else:
                u = u_new
        U[k] = u
        z = model["A"] @ z + model["B"] @ (u.reshape(-1, 1) - model["u_mean"])
        z += disturbances[k].reshape(-1, 1)
    violations = np.any(np.abs(S) > hard_bound, axis=1)
    return {
        "S": S,
        "U": U,
        "violation_steps": int(violations.sum()),
        "fallback_count": int(fallback),
        "qp_count": int(qp_count),
        "rmse": np.sqrt(np.mean((S - reference) ** 2, axis=0)),
    }


def plot_results(reference, smpc, mpc, hard_bound, path):
    t = np.arange(len(reference)) / 100.0
    labels = ("x", "y", "z")
    fig = plt.figure(figsize=(15, 9))
    grid = fig.add_gridspec(3, 2, width_ratios=[1.45, 1.0])
    for axis, label in enumerate(labels):
        ax = fig.add_subplot(grid[axis, 0])
        ax.plot(t, reference[:, axis], "k--", lw=1.2, label="3D reference")
        ax.plot(t, mpc["S"][:, axis], color="#c44536", lw=1.8, ls=":",
                label="deterministic MPC", zorder=2)
        ax.plot(t, smpc["S"][:, axis], color="#087e8b", lw=1.0,
                label="SMPC", zorder=3)
        ax.axhline(hard_bound, color="#333333", ls="-.", lw=0.8)
        ax.axhline(-hard_bound, color="#333333", ls="-.", lw=0.8)
        ax.set_ylabel(f"standardized {label}")
        ax.grid(True, alpha=0.25)
        if axis == 0:
            ax.legend(loc="best", ncol=3, fontsize=8)
        if axis == 2:
            ax.set_xlabel("simulation time (s)")
    ax3 = fig.add_subplot(grid[:, 1], projection="3d")
    ax3.plot(reference[:, 0], reference[:, 1], reference[:, 2], "k--", lw=1.5,
             label="reference")
    ax3.plot(mpc["S"][:, 0], mpc["S"][:, 1], mpc["S"][:, 2],
             color="#c44536", lw=2.0, ls=":", label="deterministic MPC")
    ax3.plot(smpc["S"][:, 0], smpc["S"][:, 1], smpc["S"][:, 2],
             color="#087e8b", lw=1.1, label="SMPC")
    ax3.set_xlabel("standardized x")
    ax3.set_ylabel("standardized y")
    ax3.set_zlabel("standardized z")
    ax3.set_title("Three-dimensional position trajectory")
    ax3.legend(loc="best", fontsize=8)
    ax3.view_init(elev=24, azim=-62)
    fig.suptitle("copyBD: real-data-trained 3D fixed-anchor model-in-the-loop control")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run_experiment(seeds=(7, 19, 31), steps=1200):
    model, E, hard_bound, scales = load_model()
    reference = moving_reference_3d(steps, hard_bound)
    results = {}
    representative = None
    for seed in seeds:
        rng = np.random.default_rng(seed)
        disturbances = rng.multivariate_normal(
            np.zeros(model["A"].shape[0]), model["Sigma_eps"], size=steps
        )
        smpc = simulate(model, E, hard_bound, reference, disturbances, 0.05)
        mpc = simulate(model, E, hard_bound, reference, disturbances, None)
        results[str(seed)] = {
            "smpc_violation_steps": smpc["violation_steps"],
            "mpc_violation_steps": mpc["violation_steps"],
            "smpc_fallback_count": smpc["fallback_count"],
            "mpc_fallback_count": mpc["fallback_count"],
            "smpc_rmse": smpc["rmse"].tolist(),
            "mpc_rmse": mpc["rmse"].tolist(),
        }
        if representative is None:
            representative = (smpc, mpc)
        print(f"seed {seed}: SMPC violations={smpc['violation_steps']} "
              f"MPC={mpc['violation_steps']} fallback={smpc['fallback_count']}/{mpc['fallback_count']}")
    out = {
        "scope": "3D fixed-anchor noisy model-in-the-loop simulation; not real flight",
        "fixed_anchor": "standardized x/y/z position channels",
        "training_flights": list(range(1, 37)),
        "seeds": list(seeds),
        "steps": int(steps),
        "d": 10,
        "horizon": 18,
        "q_weight": 80.0,
        "ru": 0.18,
        "alpha_joint": 0.05,
        "hard_bound": hard_bound,
        "physical_reference_ranges": {
            axis: [float(np.min(reference[:, i] * scales["y_scale"][i] + scales["y_offset"][i])),
                   float(np.max(reference[:, i] * scales["y_scale"][i] + scales["y_offset"][i]))]
            for i, axis in enumerate(("x_mm", "y_mm", "z_mm"))
        },
        "per_seed": results,
    }
    results_dir = HERE / "results"
    results_dir.mkdir(exist_ok=True)
    (results_dir / "copyBD_3d_boundary_trajectory.json").write_text(
        json.dumps(out, indent=2), encoding="ascii"
    )
    np.savez(results_dir / "copyBD_3d_boundary_trajectory.npz",
             reference=reference, smpc=representative[0]["S"],
             mpc=representative[1]["S"], smpc_u=representative[0]["U"],
             mpc_u=representative[1]["U"])
    plot_results(reference, representative[0], representative[1], hard_bound,
                 results_dir / "copyBD_3d_boundary_trajectory.png")
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", nargs="+", type=int, default=[7, 19, 31])
    parser.add_argument("--steps", type=int, default=1200)
    args = parser.parse_args()
    run_experiment(args.seeds, args.steps)


if __name__ == "__main__":
    main()
