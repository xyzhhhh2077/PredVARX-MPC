"""Plot copyBA multi-segment in-the-loop results: track real flights 37-41.

Two panels per segment: task-axis tracking (SMPC vs real expert flight) and
input magnitude (SMPC vs expert) with per-segment RMSE / u-std annotation.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    with open(os.path.join(HERE, "results", "copyBA_multi_segment_Ru006.json")) as f:
        out = json.load(f)
    segs = out["segs"]
    npz = np.load(os.path.join(HERE, "results", "copyBA_multi_segment_Ru006.npz"),
                  allow_pickle=True)
    per_seg = npz["per_seg"].item()
    S_smpc_all = npz["S_smpc_all"]
    S_real_all = npz["S_real_all"]
    U_smpc_all = npz["U_smpc_all"]

    # expert inputs per segment for comparison (standardized)
    import scipy.io as sio
    import sys
    sys.path.insert(0, HERE)
    from run_real_model_closed_loop import load_copyba_model, DATA_MAT
    _, _, scales = load_copyba_model()
    Dm = sio.loadmat(DATA_MAT)
    u_all, seg_all = Dm["u"], Dm["segment_id"].ravel()
    U_expert_all = []
    for s in segs:
        idx = np.where(seg_all == s)[0][:1200]
        U_expert_all.append(
            (u_all[:, idx] - scales["u_offset"].reshape(-1, 1))
            / scales["u_scale"].reshape(-1, 1))  # (m, steps)

    fig, axes = plt.subplots(len(segs), 2, figsize=(14, 3.4 * len(segs)))
    t = np.arange(1200) / 100.0
    for i, s in enumerate(segs):
        d = per_seg[s]
        ax = axes[i, 0]
        ax.plot(t, S_real_all[i][:, 0], "k--", lw=1.4, label="real flight s1")
        ax.plot(t, S_smpc_all[i][:, 0], "-", lw=1.1, color="tab:blue", label="SMPC s1")
        ax.plot(t, S_real_all[i][:, 1], "k:", lw=1.4, alpha=0.6, label="real flight s2")
        ax.plot(t, S_smpc_all[i][:, 1], "-", lw=1.1, color="tab:red", label="SMPC s2")
        ax.set_ylabel("task axis s")
        ax.set_title(f"seg {s}: tracking (RMSE "
                     f"({d['rmse_smpc'][0]:.3f},{d['rmse_smpc'][1]:.3f}), "
                     f"viol {d['viol']})")
        ax.legend(loc="best", fontsize=8)
        ax.grid(alpha=0.3)
        if i == len(segs) - 1:
            ax.set_xlabel("time (s)")

        ax2 = axes[i, 1]
        ax2.plot(t, U_smpc_all[i][:, 0], lw=1.1, color="tab:blue", label="SMPC u1")
        ax2.plot(t, U_expert_all[i][0], lw=0.8, ls=":", color="tab:blue", alpha=0.55)
        ax2.plot(t, U_smpc_all[i][:, 1], lw=1.1, color="tab:orange", label="SMPC u2")
        ax2.plot(t, U_expert_all[i][1], lw=0.8, ls=":", color="tab:orange", alpha=0.55)
        ax2.plot(t, U_smpc_all[i][:, 2], lw=1.1, color="tab:green", label="SMPC u3")
        ax2.plot(t, U_expert_all[i][2], lw=0.8, ls=":", color="tab:green", alpha=0.55)
        ax2.plot(t, U_smpc_all[i][:, 3], lw=1.1, color="tab:purple", label="SMPC u4")
        ax2.plot(t, U_expert_all[i][3], lw=0.8, ls=":", color="tab:purple", alpha=0.55)
        ax2.plot([], [], "k-", lw=1.1, label="SMPC")
        ax2.plot([], [], "k:", lw=0.8, alpha=0.55, label="expert")
        ax2.set_ylabel("u (std)")
        ax2.set_title(f"seg {s}: inputs (u std SMPC {np.round(d['u_std_smpc'], 2)}, "
                      f"expert {np.round(d['u_std_expert'], 2)})")
        ax2.grid(alpha=0.3)
        if i == len(segs) - 1:
            ax2.set_xlabel("time (s)")
        if i == 0:
            ax2.legend(loc="best", fontsize=7, ncol=2)

    fig.suptitle("copyBA: SMPC on real identified Pelican model tracks real flights 37-41\n"
                 "d=10 (10 Hz), N=18, Ru=0.06 (input matched to expert), h=2.0, "
                 "Sigma_d = d-step lifted noise; dashed = real expert flight", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_png = os.path.join(HERE, "results", "copyBA_multi_segment_Ru006_fig.png")
    fig.savefig(out_png, dpi=140)
    print("saved", out_png)


if __name__ == "__main__":
    main()
