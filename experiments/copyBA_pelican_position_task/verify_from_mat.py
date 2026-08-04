"""Independent verification of copyBA closed-loop results from .mat sources.

Recomputes from scratch (no trust in the run script's in-memory objects):
1. model params from results/copyBA_pelican_position_task_data.mat
2. identity R'P = I, E physical weights (position-dominant?)
3. real task-axis trajectories s = E' y_std for seg37-41 from the dataset
4. closed-loop RMSE / violations / saturation / u-std vs expert, recomputed
   from the saved npz trajectories and the .mat dataset
"""
import json
import os
import sys

import h5py
import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
DATA_MAT = os.path.join(HERE, "..", "copyAY_pelican_soft_preference", "data", "copyAY_pelican_dataset.mat")


def load_model_mat(path):
    """h5py reader for matlab v7.3 model file (mirrors load_copyba_model)."""
    with h5py.File(path, "r") as f:
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
        Etask = np.asarray(f["Etask"]).T  # (q, ny)
    return model, Etask


def main():
    print("=" * 70)
    print("1) MODEL .mat (h5py, transposed like load_copyba_model)")
    model, E = load_model_mat(os.path.join(RES, "copyBA_pelican_position_task_data.mat"))
    A = model["A"]
    B = model["B"]
    P = model["P"]
    R = model["R"]
    Se = model["Sigma_eps"]
    y_mean = np.asarray(model["y_mean"]).ravel()
    u_mean = np.asarray(model["u_mean"]).ravel()
    ell, ny = A.shape[0], P.shape[1]
    m_in = B.shape[1]
    print(f"shapes: A{A.shape} B{B.shape} P{P.shape} R{R.shape} E{E.shape} "
          f"Sigma_eps{Se.shape}  ell={ell} m={m_in} ny={ny}")
    print(f"R'P - I max abs: {np.abs(R.T @ P - np.eye(ell)).max():.3e}")

    # E physical weights: project E columns onto raw physical channels.
    print("-" * 70)
    print("2) TASK AXIS PHYSICAL WEIGHTS (position dominant?)")
    Dm = sio.loadmat(DATA_MAT)
    y_all = np.asarray(Dm["y"])
    seg_all = np.asarray(Dm["segment_id"]).ravel()
    u_all = np.asarray(Dm["u"])
    tr = np.isin(seg_all, np.arange(1, 37))
    y_off = y_all[:, tr].mean(1)
    y_scl = y_all[:, tr].std(1)
    u_off = u_all[:, tr].mean(1)
    u_scl = u_all[:, tr].std(1)
    names = ["px", "py", "pz", "e1", "e2", "e3", "m1", "m2", "m3", "m4"]
    # E is (ny, q) after load_copyba_model transpose; each column = task axis
    for j in range(E.shape[1]):
        w = np.abs(E[:, j]) * y_scl  # physical units (m, rad, rpm)
        print(f"axis s{j+1}: " + ", ".join(f"{n}={w[i]:.2f}" for i, n in enumerate(names)))

    print("-" * 70)
    print("3) INTEGRATOR MODES of A")
    ev = np.linalg.eigvals(A)
    n_int = int(np.sum(np.abs(ev - 1.0) < 1e-6))
    print(f"eigenvalues: {np.round(ev, 4)}  |z-1|<1e-6 count: {n_int}")
    print(f"A spectral radius: {np.max(np.abs(ev)):.4f}")

    print("-" * 70)
    print("4) REAL TASK TRAJECTORIES FROM DATASET (independent recompute)")
    steps = 1200
    for s in [37, 38, 39, 40, 41]:
        idx = np.where(seg_all == s)[0]
        print(f"seg {s}: total samples {len(idx)}")
        idx = idx[:steps]
        y_std = (y_all[:, idx] - y_off[:, None]) / y_scl[:, None]
        s_real = (E.T @ y_std).T  # (steps, q)
        print(f"  s range: s1 [{s_real[:, 0].min():+.2f},{s_real[:, 0].max():+.2f}] "
              f"s2 [{s_real[:, 1].min():+.2f},{s_real[:, 1].max():+.2f}] "
              f"std ({s_real[:, 0].std():.3f},{s_real[:, 1].std():.3f})")

    print("-" * 70)
    print("5) CLOSED-LOOP RESULTS RECOMPUTED (npz trajectories vs .mat reference)")
    npz = np.load(os.path.join(RES, "copyBA_multi_segment_Ru006.npz"), allow_pickle=True)
    per_seg = npz["per_seg"].item()
    S_s = npz["S_smpc_all"]
    S_r = npz["S_real_all"]
    U_s = npz["U_smpc_all"]
    h = 2.0
    for i, s in enumerate([37, 38, 39, 40, 41]):
        idx = np.where(seg_all == s)[0][:steps]
        y_std = (y_all[:, idx] - y_off[:, None]) / y_scl[:, None]
        s_real_ck = (E.T @ y_std).T
        U_expert = (u_all[:, idx] - u_off[:, None]) / u_scl[:, None]
        ref_max = np.abs(S_r[i] - s_real_ck).max()
        rmse = np.sqrt(np.mean((S_s[i] - S_r[i]) ** 2, axis=0))
        sat = int(np.mean(np.max(np.abs(U_s[i]) > 5.999, axis=1)) * 100)
        viol = int(np.sum(np.max(np.abs(S_s[i]), axis=1) > h))
        print(f"seg {s}: ref-recompute max diff {ref_max:.2e} | "
              f"RMSE recomputed ({rmse[0]:.3f},{rmse[1]:.3f}) "
              f"| json says ({per_seg[s]['rmse_smpc'][0]:.3f},{per_seg[s]['rmse_smpc'][1]:.3f})")
        print(f"  u std recomputed SMPC {np.round(U_s[i].std(0), 3)} expert {np.round(U_expert.std(1), 3)}")
        print(f"  sat {sat}% viol {viol}")

    print("-" * 70)
    print("6) validation block in .mat (offline task R2)")
    with h5py.File(os.path.join(RES, "copyBA_pelican_position_task_data.mat"), "r") as f:
        v = f["validation"]
        print("  task_rmse:", np.asarray(v["task_rmse"]).ravel())
        print("  task_r2:", np.asarray(v["task_r2"]).ravel())
        print("  task_persistence_rmse:", np.asarray(v["task_persistence_rmse"]).ravel())


if __name__ == "__main__":
    main()
