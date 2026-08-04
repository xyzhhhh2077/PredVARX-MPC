"""Plot copyBA in-the-loop closed-loop: track real seg37 flight vs expert."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_MAT = os.path.join(HERE, "..", "copyAY_pelican_soft_preference", "data", "copyAY_pelican_dataset.mat")


def main():
    D = np.load(os.path.join(HERE, "results", "copyBA_track_real_seg37.npz"))
    S_smpc, S_ref, U_smpc = D["S_smpc"], D["S_ref"], D["U_smpc"]
    met = json.load(open(os.path.join(HERE, "results", "copyBA_track_real_seg37_metrics.json")))

    # real expert inputs (standardized)
    import sys
    sys.path.insert(0, HERE)
    from run_real_model_closed_loop import load_copyba_model
    _, _, scales = load_copyba_model()
    Dm = sio.loadmat(DATA_MAT)
    u_real = Dm["u"]
    seg = Dm["segment_id"].ravel()
    idx = np.where(seg == 37)[0][:1200]
    U_real = (u_real[:, idx] - scales["u_offset"].reshape(-1, 1)) / scales["u_scale"].reshape(-1, 1)

    fig, axes = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True)
    t = np.arange(S_smpc.shape[0]) / 100.0

    for i, nm in enumerate(["s1 (position-dominant)", "s2"]):
        ax = axes[i]
        ax.plot(t, S_ref[:, i], "k--", lw=1.6, label="real flight seg37 (expert, reference)")
        ax.plot(t, S_smpc[:, i], "-", lw=1.2, color="tab:blue", label="SMPC on real identified model")
        ax.set_ylabel(nm)
        ax.legend(loc="best", fontsize=9)
        ax.grid(alpha=0.3)
        rmse = met["rmse_smpc"][i]
        ax.set_title(f"RMSE {rmse:.3f}", loc="right", fontsize=9)

    ax = axes[2]
    for j in range(4):
        ax.plot(t, U_smpc[:, j], lw=1.1, color=f"C{j}", label=f"SMPC u{j+1}")
        ax.plot(t, U_real[j], lw=0.7, ls=":", color=f"C{j}", alpha=0.5)
    ax.plot([], [], "k-", lw=1.1, label="SMPC")
    ax.plot([], [], "k:", lw=0.7, alpha=0.5, label="expert (real)")
    ax.set_ylabel("u (std units)")
    ax.set_xlabel("time (s)")
    ax.legend(loc="best", fontsize=8, ncol=2)
    ax.grid(alpha=0.3)

    fig.suptitle("copyBA: SMPC closed loop on REAL identified Pelican model, tracking real flight seg37\n"
                 f"task axis = position-dominant E; d={met['d']} (10 Hz decisions), Ru={met['ru']}, "
                 f"u saturation {met['sat_pct']:.0f}%", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(HERE, "results", "copyBA_track_real_seg37_fig.png")
    fig.savefig(out, dpi=140)
    print("saved", out)


if __name__ == "__main__":
    main()
