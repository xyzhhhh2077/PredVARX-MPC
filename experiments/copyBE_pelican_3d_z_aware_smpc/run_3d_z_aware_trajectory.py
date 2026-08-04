"""copyBE second-order z-aware fixed-xyz Pelican SMPC experiment."""

import argparse
import importlib.util
import json
from pathlib import Path

import cvxpy as cp
import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
BASE_PATH = HERE.parent / "copyBD_pelican_3d_fixed_anchor_smpc" / "run_3d_boundary_trajectory.py"
_spec = importlib.util.spec_from_file_location("copybd_base", BASE_PATH)
base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(base)
MODEL_MAT = HERE / "results" / "copyBE_pelican_3d_z_aware_data.mat"


def collective_motor_direction(m=4):
    return np.ones((int(m), 1)) / np.sqrt(float(m))


def second_order_companion(A1, A2, B):
    n = A1.shape[0]
    A_aug = np.block([[A1, A2], [np.eye(n), np.zeros((n, n))]])
    B_aug = np.vstack([B, np.zeros((n, B.shape[1]))])
    return A_aug, B_aug


def augment_output_loading(P):
    return np.hstack([P, np.zeros_like(P)])


def load_model():
    with h5py.File(MODEL_MAT, "r") as f:
        m = f["model"]
        model = {key: np.asarray(m[key]).T for key in (
            "A", "B", "P", "R", "Sigma_eps_base", "y_mean", "u_mean"
        )}
        E = np.asarray(f["Eanchor"]).T
        hard_bound = float(np.asarray(f["task_bound"]).ravel()[0])
        scales = {
            "y_offset": np.asarray(f["y_offset"]).ravel(),
            "y_scale": np.asarray(f["y_scale"]).ravel(),
        }
    n = model["A"].shape[0]
    ell = n // 2
    Sigma = np.zeros((n, n))
    Sigma[:ell, :ell] = model.pop("Sigma_eps_base")
    model["Sigma_eps"] = Sigma
    return model, E, hard_bound, scales


class ZAwareSMPC:
    def __init__(self, model, E, d=10, horizon=18, q_weight=80.0,
                 ru=0.18, task_bound=3.0, alpha_joint=0.05):
        self.model = model
        self.E = E
        self.d = int(d)
        self.N = int(horizon)
        self.nu = model["B"].shape[1]
        self.Ad, self.Bd = base.mstep_lift(model["A"], model["B"], self.d)
        self.Sigma_d = base.lift_noise_covariance(
            model["A"], model["Sigma_eps"], self.d
        )
        self.Q = base.task_cost_matrix(E, q_weight)
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
            Gj = np.zeros((model["P"].shape[0], self.N * self.nu))
            for i in range(j):
                Gj[:, i * self.nu:(i + 1) * self.nu] = (
                    model["P"] @ np.linalg.matrix_power(self.Ad, j - 1 - i) @ self.Bd
                )
            self.G.append(Gj)
        self.Sigma_z = []
        Sz = np.zeros_like(self.Ad)
        for _ in range(self.N):
            Sz = self.Ad @ Sz @ self.Ad.T + self.Sigma_d
            self.Sigma_z.append((Sz + Sz.T) / 2.0)

    def step(self, z_aug, reference_horizon):
        z_aug = np.asarray(z_aug).reshape(-1, 1)
        references = np.asarray(reference_horizon)
        if references.shape != (self.N, self.E.shape[1]):
            raise ValueError("Reference horizon must have shape (N, task_axes)")
        U0 = np.tile(self.u_mean, (self.N, 1))
        U = cp.Variable((self.N * self.nu, 1))
        cost = cp.quad_form(U - U0, np.kron(np.eye(self.N), self.Ru))
        constraints = [np.tile(self.u_min, self.N).reshape(-1, 1) <= U,
                       U <= np.tile(self.u_max, self.N).reshape(-1, 1)]
        z_quantile = None
        if self.alpha_joint is not None:
            risk = base.per_face_risk(self.alpha_joint, self.E.shape[1], self.N)
            z_quantile = norm.ppf(1.0 - risk)
        for j in range(self.N):
            mu = self.M[j] @ z_aug + self.y_mean + self.G[j] @ (U - U0)
            target = self.E @ references[j].reshape(-1, 1)
            cost += cp.quad_form(mu - target, self.Q)
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
    controller = ZAwareSMPC(model, E, d, horizon, q_weight, ru,
                            hard_bound, alpha_joint)
    ell = model["A"].shape[0] // 2
    y0 = model["y_mean"] + E @ reference[0].reshape(-1, 1)
    z0 = model["R"].T @ (y0 - model["y_mean"])
    z_aug = np.vstack([z0, z0])
    steps = len(reference)
    S = np.zeros((steps, 3))
    U = np.zeros((steps, model["B"].shape[1]))
    u = model["u_mean"].ravel()
    fallback = 0
    for k in range(steps):
        y = model["y_mean"] + model["P"] @ z_aug
        S[k] = (E.T @ y).ravel()
        if k % d == 0:
            refs = base.future_reference_horizon(reference, k, d, horizon)
            u_new, _ = controller.step(z_aug, refs)
            if u_new is None:
                fallback += 1
                u = model["u_mean"].ravel()
            else:
                u = u_new
        U[k] = u
        z_aug = model["A"] @ z_aug + model["B"] @ (
            u.reshape(-1, 1) - model["u_mean"]
        )
        z_aug[:ell] += disturbances[k].reshape(-1, 1)
    violations = np.any(np.abs(S) > hard_bound, axis=1)
    return {"S": S, "U": U, "violation_steps": int(violations.sum()),
            "fallback_count": fallback,
            "rmse": np.sqrt(np.mean((S - reference) ** 2, axis=0))}


def plot_results(reference, smpc, mpc, hard_bound, path):
    base.plot_results(reference, smpc, mpc, hard_bound, path)


def run_experiment(seeds=(7, 19, 31), steps=1200):
    model, E, hard_bound, scales = load_model()
    reference = base.moving_reference_3d(steps, hard_bound)
    results = {}
    representative = None
    ell = model["A"].shape[0] // 2
    for seed in seeds:
        rng = np.random.default_rng(seed)
        disturbances = rng.multivariate_normal(
            np.zeros(ell), model["Sigma_eps"][:ell, :ell], size=steps
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
        "scope": "second-order z-aware 3D frozen-model-in-the-loop simulation",
        "model_order": 2, "seeds": list(seeds), "steps": int(steps),
        "d": 10, "horizon": 18, "q_weight": 80.0, "ru": 0.18,
        "alpha_joint": 0.05, "hard_bound": hard_bound, "per_seed": results,
        "physical_reference_ranges": {
            axis: [float(np.min(reference[:, i] * scales["y_scale"][i] + scales["y_offset"][i])),
                   float(np.max(reference[:, i] * scales["y_scale"][i] + scales["y_offset"][i]))]
            for i, axis in enumerate(("x_mm", "y_mm", "z_mm"))
        },
    }
    r = HERE / "results"
    r.mkdir(exist_ok=True)
    (r / "copyBE_3d_z_aware_trajectory.json").write_text(
        json.dumps(out, indent=2), encoding="ascii"
    )
    np.savez(r / "copyBE_3d_z_aware_trajectory.npz", reference=reference,
             smpc=representative[0]["S"], mpc=representative[1]["S"],
             smpc_u=representative[0]["U"], mpc_u=representative[1]["U"])
    plot_results(reference, representative[0], representative[1], hard_bound,
                 r / "copyBE_3d_z_aware_trajectory.png")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, nargs="+", default=[7, 19, 31])
    parser.add_argument("--steps", type=int, default=1200)
    args = parser.parse_args()
    run_experiment(tuple(args.seeds), args.steps)
