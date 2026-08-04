"""B: closed-loop SMPC on the REAL identified Pelican model (model in-the-loop).

Plant = copyBA model identified from real Waterloo Pelican flight data
(results/copyBA_pelican_position_task_data.mat). The plant is simulated as

    z_{k+1} = A z_k + B (u_k - u_mean) + eps_k,   eps ~ N(0, Sigma_eps)
    y_std   = y_mean + P z_k

i.e. the real-data model itself is the controlled object (with its estimated
process noise). The SMPC QP is the Python port of centered_smpc_step.m
(already verified 1-step against MATLAB).

Multirate control: the QP is solved every d steps (d=10 -> 10 Hz decision
rate on the 100 Hz plant) so the N=18 prediction horizon covers ~1.8 s,
matching the Pelican position time scale (4 integrator modes in A).

Reference: piecewise task-axis setpoints s_ref in task space s = E' y_std.
The reference segments are chosen inside the training task-axis range.

Initial state: real validation segment 37 start point (s0 = (0.01, -1.22)),
a real flight state, not the model equilibrium.

Baseline: open-loop hold u = u_mean (no feedback), same plant, same initial
state, showing the plant drifts without control.

Usage: python run_real_model_closed_loop.py [--steps 1200] [--seed 7] [--d 10]
"""

import argparse
import json
import os

import cvxpy as cp
import numpy as np
import scipy.io as sio
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_MAT = os.path.join(HERE, "..", "copyAY_pelican_soft_preference", "data", "copyAY_pelican_dataset.mat")


def load_copyba_model():
    """Load the copyBA real-Pelican model (h5py transposes MATLAB matrices)."""
    with h5py.File(os.path.join(HERE, "results", "copyBA_pelican_position_task_data.mat"), "r") as f:
        m = f["model"]
        model = {
            "A": np.asarray(m["A"]).T,
            "B": np.asarray(m["B"]).T,
            "P": np.asarray(m["P"]).T,
            "R": np.asarray(m["R"]).T,
            "Sigma_eps": np.asarray(m["Sigma_eps"]).T,
            "y_mean": np.asarray(m["y_mean"]).T,
            "u_mean": np.asarray(m["u_mean"]).T,
        }
        Etask = np.asarray(f["Etask"]).T
    # standardization statistics = training segments 1:36 (same as MATLAB)
    D = sio.loadmat(DATA_MAT)
    y, u, seg = D["y"], D["u"], D["segment_id"].ravel()
    tr = np.isin(seg, np.arange(1, 37))
    scales = {
        "y_offset": y[:, tr].mean(1),
        "y_scale": y[:, tr].std(1),
        "u_offset": u[:, tr].mean(1),
        "u_scale": u[:, tr].std(1),
    }
    return model, Etask, scales


class SMPC:
    """Centered linear-Gaussian SMPC QP, port of centered_smpc_step.m.

    Works on the multirate plant model (A_d, B_d = d-step lifted model); the
    horizon j then counts decision steps, each of physical length d.
    The lifted noise covariance Sigma_d accumulates the d physical steps:
    Sigma_d = sum_{i=0}^{d-1} A^i Sigma_eps A'^i, so a decision-step
    prediction carries the full d-step disturbance variance.
    """

    def __init__(self, model, opt, d=1):
        self.model = model
        self.opt = opt
        self.ny = model["P"].shape[0]
        self.nu = model["B"].shape[1]
        self.N = opt["N"]
        nq = opt["H"].shape[0]
        self.nq = nq
        A, B, P, Sigma_eps = model["A"], model["B"], model["P"], model["Sigma_eps"]
        # lifted (d-step) disturbance covariance for the decision step
        Sigma_d = np.zeros_like(Sigma_eps)
        Ai = np.eye(A.shape[0])
        for _ in range(d):
            Sigma_d += Ai @ Sigma_eps @ Ai.T
            Ai = Ai @ A
        self.M = [P @ np.linalg.matrix_power(A, j) for j in range(1, self.N + 1)]
        self.G = []
        for j in range(1, self.N + 1):
            Gj = np.zeros((self.ny, self.N * self.nu))
            for i in range(j):
                Gj[:, i * self.nu:(i + 1) * self.nu] = \
                    P @ np.linalg.matrix_power(A, j - 1 - i) @ B
            self.G.append(Gj)
        self.Sigma_z_list = []
        Sz = np.zeros_like(A)
        for _ in range(self.N):
            Sz = A @ Sz @ A.T + Sigma_d
            self.Sigma_z_list.append(Sz)
        self.y_mean = model["y_mean"].reshape(-1, 1)
        self.u_mean = model["u_mean"].reshape(-1, 1)

    def step(self, y_std, r_full):
        o = self.opt
        N, nu = self.N, self.nu
        y_mean, u_mean = self.y_mean, self.u_mean
        z = self.model["R"].T @ (y_std.reshape(-1, 1) - y_mean)
        U0 = np.tile(u_mean, (N, 1))
        U = cp.Variable((N * nu, 1))
        cost = 0
        constraints = []
        from scipy.stats import norm
        for j in range(self.N):
            ej0 = self.M[j] @ z + y_mean - self.G[j] @ U0 - r_full.reshape(-1, 1)
            cost += cp.quad_form(self.G[j] @ U + ej0, o["Q"])
            Sigma_y = self.model["P"] @ self.Sigma_z_list[j] @ self.model["P"].T
            for q in range(self.nq):
                hq = o["H"][q, :].reshape(1, -1)
                mu_j = self.M[j] @ z + y_mean + self.G[j] @ (U - U0)
                kappa = norm.ppf(1.0 - o["alpha_joint"] / (2.0 * self.nq * N)) * \
                    float(np.sqrt((hq @ Sigma_y @ hq.T)[0, 0]))
                constraints.append(hq @ mu_j <= o["h"][q] - kappa)
                constraints.append(hq @ mu_j >= -o["h"][q] + kappa)
        constraints += [o["u_min"] <= U, U <= o["u_max"]]
        cost += cp.quad_form(U - U0, o["RuBar"])
        prob = cp.Problem(cp.Minimize(cost), constraints)
        prob.solve(solver=cp.OSQP, verbose=False)
        if prob.status not in ("optimal", "optimal_inaccurate"):
            return None, prob.status
        Uv = U.value
        u_first = np.clip(Uv[:nu].ravel(), o["u_min"], o["u_max"])
        return u_first, prob.status


def r_full_of(Etask, s_ref):
    """Map task-axis reference to full y_std reference: r = E (E'E)^{-1} s."""
    return Etask @ np.linalg.solve(Etask.T @ Etask, np.asarray(s_ref, dtype=float))


def mstep_lift(A, B, d):
    """d-step lifted model: z_{k+d} = A_d z_k + B_d (u - u_mean), u held."""
    Ad = np.linalg.matrix_power(A, d)
    Bd = np.zeros_like(B)
    Ai = np.eye(A.shape[0])
    for _ in range(d):
        Bd += Ai @ B
        Ai = Ai @ A
    return Ad, Bd


def plant_step(model, z, u, rng):
    """One physical-step plant evolution with process noise."""
    eps = rng.multivariate_normal(np.zeros(model["A"].shape[0]), model["Sigma_eps"]).reshape(-1, 1)
    z = model["A"] @ z + model["B"] @ (u.reshape(-1, 1) - model["u_mean"]) + eps
    return z


def run_closed_loop(steps, seed, d, h, qmul, init_real, ref_levels, out_prefix):
    model, Etask, scales = load_copyba_model()
    p, m, q = 10, 4, 2
    ell = model["A"].shape[0]

    # multirate lifted model for the controller
    Ad, Bd = mstep_lift(model["A"], model["B"], d)
    ctrl_model = dict(model)
    ctrl_model["A"] = Ad
    ctrl_model["B"] = Bd

    Qw = np.diag([1.0] * 6 + [0.05] * 4)
    Q = Qw / 6.0 * qmul
    N = 18
    opt = {
        "N": N, "Q": Q,
        "Ru": 1e-4 * np.eye(m),
        "RuBar": np.kron(np.eye(N), 1e-4 * np.eye(m)),
        "u_min": np.full(m, -6.0), "u_max": np.full(m, 6.0),
        "H": Etask.T, "h": np.full(q, h), "alpha_joint": 0.05,
    }
    smpc = SMPC(ctrl_model, opt, d=d)

    rng = np.random.default_rng(seed)
    z = np.zeros((ell, 1))
    if init_real:
        D = sio.loadmat(DATA_MAT)
        y, seg = D["y"], D["segment_id"].ravel()
        idx = np.where(seg == 37)[0]
        y0 = y[:, idx[0]]
        y0_std = (y0 - scales["y_offset"]) / scales["y_scale"]
        z = model["R"].T @ (y0_std.reshape(-1, 1) - model["y_mean"])
    y_std = model["y_mean"] + model["P"] @ z

    # reference timeline (physical steps)
    S_ref = np.zeros((steps, q))
    for k in range(steps):
        seg_i = min(k // (steps // len(ref_levels)), len(ref_levels) - 1)
        S_ref[k] = ref_levels[seg_i]

    S_smpc = np.zeros((steps, q))
    U_smpc = np.zeros((steps, m))
    S_ol = np.zeros((steps, q))
    U_ol = np.zeros((steps, m))
    u_ol = model["u_mean"].ravel()  # open-loop hold baseline

    z_ol = z.copy()
    y_ol = y_std.copy()

    for k in range(steps):
        s_now = np.asarray(Etask.T @ y_std).ravel()
        S_smpc[k] = s_now
        if k % d == 0:
            r_full = r_full_of(Etask, S_ref[k])
            u, status = smpc.step(y_std.ravel(), r_full)
            if u is None:
                print(f"k={k}: QP {status}")
                u = np.zeros(m)
        U_smpc[k] = u
        z = plant_step(model, z, u, rng)
        y_std = model["y_mean"] + model["P"] @ z

        # open-loop baseline on its own plant (same noise stream drawn separately)
        U_ol[k] = u_ol
        z_ol = plant_step(model, z_ol, u_ol, rng)
        y_ol = model["y_mean"] + model["P"] @ z_ol
        S_ol[k] = np.asarray(Etask.T @ y_ol).ravel()

    # metrics
    rmse_s = np.sqrt(np.mean((S_smpc - S_ref) ** 2, axis=0))
    rmse_ol = np.sqrt(np.mean((S_ol - S_ref) ** 2, axis=0))
    sat = int(np.mean(np.max(np.abs(U_smpc) > 5.999, axis=1)) * 100)
    viol = 0
    for k in range(steps):
        if np.max(np.abs(S_smpc[k])) > h:
            viol += 1
    metrics = {
        "plant": "copyBA real-Pelican identified model (in-the-loop)",
        "task_axis": "position-dominant (pos 1.0 / euler 0.3 / motors 0.05, pref 0.99)",
        "d": d, "N": N, "h": h, "qmul": qmul,
        "init": "real seg37 start" if init_real else "model equilibrium",
        "steps": steps, "seed": seed,
        "rmse_smpc": [float(x) for x in rmse_s],
        "rmse_openloop": [float(x) for x in rmse_ol],
        "violation_steps": viol,
        "input_saturation_pct": sat,
    }
    with open(os.path.join(HERE, "results", out_prefix + ".json"), "w") as f:
        json.dump(metrics, f, indent=1)
    np.savez(os.path.join(HERE, "results", out_prefix + ".npz"),
             S_smpc=S_smpc, S_ref=S_ref, U_smpc=U_smpc, S_ol=S_ol, U_ol=U_ol)
    return metrics


def run_multi_segment(segs, steps, seed, d, h, qmul, Ru, out_prefix):
    """Closed loop per real flight segment: reference = the segment's own
    real task-axis trajectory, initial state = real segment start.
    Also reports SMPC vs expert input magnitudes."""
    model, Etask, scales = load_copyba_model()
    p, m, q = 10, 4, 2
    ell = model["A"].shape[0]
    Ad, Bd = mstep_lift(model["A"], model["B"], d)
    ctrl_model = dict(model)
    ctrl_model["A"] = Ad
    ctrl_model["B"] = Bd
    Qw = np.diag([1.0] * 6 + [0.05] * 4)
    Q = Qw / 6.0 * qmul
    N = 18
    opt = {
        "N": N, "Q": Q,
        "Ru": Ru * np.eye(m),
        "RuBar": np.kron(np.eye(N), Ru * np.eye(m)),
        "u_min": np.full(m, -6.0), "u_max": np.full(m, 6.0),
        "H": Etask.T, "h": np.full(q, h), "alpha_joint": 0.05,
    }
    smpc = SMPC(ctrl_model, opt, d=d)
    D = sio.loadmat(DATA_MAT)
    y, u_all, seg_all = D["y"], D["u"], D["segment_id"].ravel()
    per_seg = {}
    S_smpc_all, S_real_all, U_smpc_all = [], [], []
    for segid in segs:
        idx = np.where(seg_all == segid)[0][:steps]
        if len(idx) < 100:
            print(f"seg {segid}: only {len(idx)} samples, skip")
            continue
        S_real = np.asarray(Etask.T @ ((y[:, idx] - scales["y_offset"].reshape(-1, 1))
                                       / scales["y_scale"].reshape(-1, 1))).T  # (steps, q)
        U_expert = ((u_all[:, idx] - scales["u_offset"].reshape(-1, 1))
                    / scales["u_scale"].reshape(-1, 1))  # (m, steps)
        y0_std = (y[:, idx[0]] - scales["y_offset"]) / scales["y_scale"]
        z = model["R"].T @ (y0_std.reshape(-1, 1) - model["y_mean"])
        y_std = model["y_mean"] + model["P"] @ z
        rng = np.random.default_rng(seed)
        S_s, U_s = [], []
        for k in range(steps):
            if k % d == 0:
                r_full = r_full_of(Etask, S_real[k])
                u, status = smpc.step(y_std.ravel(), r_full)
                if u is None:
                    print(f"seg{segid} k={k}: QP {status}")
                    u = np.zeros(m)
            z = plant_step(model, z, u, rng)
            y_std = model["y_mean"] + model["P"] @ z
            S_s.append(np.asarray(Etask.T @ y_std).ravel())
            U_s.append(u)
        S_s = np.array(S_s)
        U_s = np.array(U_s)
        rmse = np.sqrt(np.mean((S_s - S_real) ** 2, axis=0))
        sat = int(np.mean(np.max(np.abs(U_s) > 5.999, axis=1)) * 100)
        viol = int(np.sum(np.max(np.abs(S_s), axis=1) > h))
        per_seg[int(segid)] = {
            "rmse_smpc": [float(x) for x in rmse],
            "sat_pct": sat, "viol": viol,
            "u_std_smpc": [float(x) for x in U_s.std(0)],
            "u_std_expert": [float(x) for x in U_expert.std(1)],
            "s_std": [float(x) for x in S_real.std(0)],
        }
        S_smpc_all.append(S_s)
        S_real_all.append(S_real)
        U_smpc_all.append(U_s)
        print(f"seg{segid}: rmse=({rmse[0]:.3f},{rmse[1]:.3f}) sat={sat}% viol={viol} "
              f"u_std smpc={np.round(U_s.std(0),2)} expert={np.round(U_expert.std(1),2)}",
              flush=True)
    out = {"segs": [int(s) for s in segs], "steps": steps, "seed": seed, "d": d,
           "h": h, "qmul": qmul, "Ru": Ru, "per_seg": per_seg}
    with open(os.path.join(HERE, "results", out_prefix + ".json"), "w") as f:
        json.dump(out, f, indent=1)
    np.savez(os.path.join(HERE, "results", out_prefix + ".npz"),
             per_seg=per_seg, segs=np.asarray([int(s) for s in segs]),
             S_smpc_all=S_smpc_all, S_real_all=S_real_all, U_smpc_all=U_smpc_all)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", type=int, default=1200)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--d", type=int, default=10, help="decision interval (physical steps per QP call)")
    ap.add_argument("--h", type=float, default=2.0, help="task-axis bound")
    ap.add_argument("--qmul", type=float, default=10.0, help="Q multiplier")
    ap.add_argument("--init-equilibrium", action="store_true", help="start at model equilibrium instead of real seg37")
    ap.add_argument("--prefix", default="copyBA_real_model_closed_loop")
    args = ap.parse_args()

    ref_levels = [[0.8, 0.0], [-0.8, 0.0], [0.4, 0.0]]
    m = run_closed_loop(args.steps, args.seed, args.d, args.h, args.qmul,
                        not args.init_equilibrium, ref_levels, args.prefix)
    print(json.dumps(m, indent=1))


if __name__ == "__main__":
    main()
