"""Closed-loop SMPC on gym-pybullet-drones CF2X using the identified copyAZ model.

Loads results/copyAZ_pybullet_model.mat (identified by copyAZ_pybullet_identify.m
through the copyAY pipeline), then runs a centered SMPC QP (Python port of
centered_smpc_step.m) against the live CtrlAviary simulator.

Reference: 3-segment position setpoints (x/y moves ~0.5 m, z hold at takeoff
height), expressed in the task space s = E' y_std.

Baseline: same reference with constant hover RPM (open loop). The CF2X has no
onboard attitude stabilizer, so the baseline drifts -- the SMPC must use the
identified model to hold position, which is the phase-2 method validation.

Usage: python run_pybullet_smpc_closed_loop.py [--steps 1200] [--seed 7]
"""
import argparse
import json
import os
import time

import cvxpy as cp
import numpy as np
import scipy.io as sio

from gym_pybullet_drones.envs.CtrlAviary import CtrlAviary
from gym_pybullet_drones.utils.enums import DroneModel, Physics

HOVER_RPM = 22000.0
HERE = os.path.dirname(os.path.abspath(__file__))


class SMPC:
    """Centered linear-Gaussian SMPC QP, port of centered_smpc_step.m."""

    def __init__(self, model, opt):
        self.model = model
        self.opt = opt
        self.ny = model["P"].shape[0]
        self.nu = model["B"].shape[1]
        self.N = opt["N"]
        nq = opt["H"].shape[0]
        self.nq = nq
        A, B, P, Sigma_eps = model["A"], model["B"], model["P"], model["Sigma_eps"]
        self.M = [P @ np.linalg.matrix_power(A, j) for j in range(1, self.N + 1)]
        self.G = []
        for j in range(1, self.N + 1):
            Gj = np.zeros((self.ny, self.N * self.nu))
            for i in range(j):
                Gj[:, i * self.nu:(i + 1) * self.nu] = \
                    P @ np.linalg.matrix_power(A, j - 1 - i) @ B
            self.G.append(Gj)
        # covariance propagation
        self.Sigma_z_list = []
        Sz = np.zeros_like(A)
        for _ in range(self.N):
            Sz = A @ Sz @ A.T + Sigma_eps
            self.Sigma_z_list.append(Sz)
        self.y_mean = model["y_mean"].reshape(-1, 1)
        self.u_mean = model["u_mean"].reshape(-1, 1)

    def step(self, y_std, r_full):
        """Solve one QP. Returns (u_first_std, info)."""
        o = self.opt
        N, nu = self.N, self.nu
        y_mean, u_mean = self.y_mean, self.u_mean
        z = self.model["R"].T @ (y_std.reshape(-1, 1) - y_mean)
        U0 = np.tile(u_mean, (N, 1))
        U = cp.Variable((N * nu, 1))
        # J = sum_j ||M_j z + y_mean - G_j U - r||_Q^2 + ||U - U0||_Ru^2
        cost = 0
        constraints = []
        from scipy.stats import norm
        for j in range(self.N):
            ej0 = self.M[j] @ z + y_mean - self.G[j] @ U0 - r_full.reshape(-1, 1)
            cost += cp.quad_form(self.G[j] @ U - self.G[j] @ U0 + ej0, o["Q"])
            # chance constraints on tracked axes only (AL geometry contract)
            Sigma_y = self.model["P"] @ self.Sigma_z_list[j] @ self.model["P"].T
            for q in range(self.nq):
                hq = o["H"][q, :].reshape(1, -1)
                mu_j = self.M[j] @ z + y_mean - self.G[j] @ (U - U0)
                hmu = hq @ mu_j
                hSh = hq @ Sigma_y @ hq.T
                risk = o["alpha_joint"] / (2 * self.nq * self.N)
                kappa = norm.ppf(1 - risk) * np.sqrt(hSh)
                constraints += [hmu <= o["h"][q] - kappa,
                                hmu >= -o["h"][q] + kappa]
        cost += cp.quad_form(U - U0, o["RuBar"])
        prob = cp.Problem(cp.Minimize(cost), constraints)
        prob.solve(solver=cp.OSQP, verbose=False)
        if prob.status not in ("optimal", "optimal_inaccurate"):
            return None, {"status": prob.status, "cost": np.nan}
        Uv = U.value.reshape(-1)
        # input bounds
        u_std = np.clip(Uv[:nu], o["u_min"], o["u_max"])
        return u_std, {"status": prob.status, "cost": prob.value}


def build_reference(n_steps, q, seed=7):
    """3-segment reference on the position-dominated task axis s1 only.

    s2 is the motor-speed combination axis (see E): no reference is defined
    there (soft-preference, zero), matching the copyAR tracked-axis design.
    """
    rng = np.random.default_rng(seed)
    ref = np.zeros((q, n_steps))
    seg = n_steps // 3
    levels = np.array([0.35, -0.3, 0.15])
    for j in range(3):
        ref[0, j * seg:(j + 1) * seg] = levels[j]
    return ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=1200)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--ctrl_freq", type=int, default=100)
    args = ap.parse_args()

    # ---- load identified model (MATLAB -v7.3 -> HDF5) ----
    import h5py
    with h5py.File(os.path.join(HERE, "results", "copyAZ_pybullet_model.mat"),
                   "r") as f:
        m = f["model"]
        model = {k: np.asarray(m[k]).T for k in
                 ("A", "B", "P", "R", "Sigma_eps")}
        model["y_mean"] = np.asarray(m["y_mean"]).T
        model["u_mean"] = np.asarray(m["u_mean"]).T
        E = np.asarray(f["E"]).T
        scales_raw = {k: np.asarray(f["scales"][k]).ravel() for k in
                      ("y_offset", "y_scale", "u_offset", "u_scale")}
    scales = {"y_offset": scales_raw["y_offset"],
              "y_scale": scales_raw["y_scale"],
              "u_offset": scales_raw["u_offset"],
              "u_scale": scales_raw["u_scale"]}
    q = E.shape[1]

    # ---- controller options ----
    Qw = np.diag([1, 1, 1, 1, 1, 1, 0.05, 0.05, 0.05, 0.05])
    Q = Qw / 6.0
    Ru = 1e-4 * np.eye(4)
    u_min = (0.85 * HOVER_RPM - scales["u_offset"]) / scales["u_scale"]
    u_max = (1.15 * HOVER_RPM - scales["u_offset"]) / scales["u_scale"]
    task_limit = 0.8
    opt = {"N": 18, "Q": Q, "Ru": Ru, "RuBar": np.kron(np.eye(18), Ru),
           "u_min": u_min, "u_max": u_max,
           "H": E.T, "h": task_limit * np.ones(q), "alpha_joint": 0.05}
    smpc = SMPC(model, opt)

    T = args.steps
    reference = build_reference(T, q, args.seed)

    env = CtrlAviary(drone_model=DroneModel.CF2X, num_drones=1,
                     physics=Physics.PYB, pyb_freq=4 * args.ctrl_freq,
                     ctrl_freq=args.ctrl_freq, gui=False, user_debug_gui=False)
    env.reset()
    # settle at hover
    for _ in range(2 * args.ctrl_freq):
        env.step(np.array([[HOVER_RPM] * 4]))

    y_rec = np.zeros((10, T)); u_rec = np.zeros((4, T)); s_rec = np.zeros((q, T))
    exit_ok = np.zeros(T, dtype=bool)
    t0 = time.time()
    for k in range(T):
        obs, *_ = env.step(np.array([[HOVER_RPM] * 4]))
        o = obs.reshape(-1)
        yk = np.hstack([o[0:3], o[7:10], o[16:20]]).reshape(-1, 1)
        y_std = (yk - scales["y_offset"].reshape(-1, 1)) / \
            scales["y_scale"].reshape(-1, 1)
        y_rec[:, k] = yk.ravel()
        s_rec[:, k] = (E.T @ y_std).ravel()
        desired = reference[:, min(k + 1, T - 1)]
        Gram = E.T @ E
        r_full = E @ np.linalg.solve(Gram, desired)
        u_std, info = smpc.step(y_std.ravel(), r_full.ravel())
        if u_std is None:
            u_std = np.tile(model["u_mean"].ravel() / 1.0, 1)
            u_std = np.clip(u_std, u_min, u_max)
        else:
            exit_ok[k] = True
        u_cmd = scales["u_offset"] + scales["u_scale"] * u_std
        u_rec[:, k] = u_cmd
        env.step(u_cmd.reshape(1, 4))
    env.close()

    # ---- baseline: open-loop hover on the same reference task ----
    env2 = CtrlAviary(drone_model=DroneModel.CF2X, num_drones=1,
                      physics=Physics.PYB, pyb_freq=4 * args.ctrl_freq,
                      ctrl_freq=args.ctrl_freq, gui=False, user_debug_gui=False)
    env2.reset()
    for _ in range(2 * args.ctrl_freq):
        env2.step(np.array([[HOVER_RPM] * 4]))
    y_base = np.zeros((10, T)); s_base = np.zeros((q, T))
    for k in range(T):
        obs, *_ = env2.step(np.array([[HOVER_RPM] * 4]))
        o = obs.reshape(-1)
        yk = np.hstack([o[0:3], o[7:10], o[16:20]]).reshape(-1, 1)
        y_std = (yk - scales["y_offset"].reshape(-1, 1)) / \
            scales["y_scale"].reshape(-1, 1)
        y_base[:, k] = yk.ravel()
        s_base[:, k] = (E.T @ y_std).ravel()
    env2.close()

    # ---- metrics (warmup 200 steps) ----
    warm = slice(200, T)
    err = s_rec[:, warm] - reference[:, warm]
    errb = s_base[:, warm] - reference[:, warm]
    rmse = np.sqrt(np.mean(err ** 2, axis=1))
    rmse_b = np.sqrt(np.mean(errb ** 2, axis=1))
    improvement = 1 - rmse / rmse_b
    qp_success = exit_ok[200:].mean()
    sat = np.mean(np.abs(u_rec) >= (0.15 * HOVER_RPM), axis=1)

    print("--- closed-loop result ---")
    print("SMPC RMSE:", np.round(rmse, 4), " baseline:", np.round(rmse_b, 4))
    print("improvement:", np.round(100 * improvement, 1), "%")
    print("QP success:", qp_success, " input saturation:", np.round(sat, 3))

    out = os.path.join(HERE, "results")
    os.makedirs(out, exist_ok=True)
    sio.savemat(os.path.join(out, "copyAZ_pybullet_closed_loop_data.mat"),
                {"y": y_rec, "u": u_rec, "s": s_rec, "reference": reference,
                 "y_base": y_base, "s_base": s_base, "rmse": rmse,
                 "rmse_baseline": rmse_b, "improvement": improvement,
                 "qp_success": qp_success, "input_saturation": sat})
    meta = {"steps": T, "ctrl_freq": args.ctrl_freq, "seed": args.seed,
            "N": opt["N"], "task_limit": task_limit,
            "note": "phase-2 method validation: CF2X != Pelican plant"}
    with open(os.path.join(out, "copyAZ_pybullet_closed_loop.json"), "w") as fp:
        json.dump(meta, fp, indent=2)
    print("saved", os.path.join(out, "copyAZ_pybullet_closed_loop_data.mat"))


if __name__ == "__main__":
    main()
