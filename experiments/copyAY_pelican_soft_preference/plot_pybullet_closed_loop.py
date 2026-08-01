"""Plot copyAZ pybullet closed-loop results: SMPC vs open-loop baseline."""
import os
import scipy.io as sio
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
d = sio.loadmat(os.path.join(HERE, "results", "copyAZ_pybullet_closed_loop_data.mat"))
s = d["s"]; ref = d["reference"]; s_b = d["s_base"]
y = d["y"]; y_b = d["y_base"]; u = d["u"]
T = s.shape[1]; q = s.shape[0]
warm = slice(200, T)

fig, axes = plt.subplots(4, 1, figsize=(11, 13), sharex=True)
tt = np.arange(T) / 100.0

for i in range(q):
    ax = axes[i]
    ax.plot(tt, ref[i], "k--", lw=1.2, label="reference (task)")
    ax.plot(tt, s[i], "tab:blue", lw=1.0, label="SMPC closed-loop")
    ax.plot(tt, s_b[i], "tab:red", lw=0.9, alpha=0.8, label="open-loop hover")
    ax.set_ylabel(f"task axis s{i+1} (std)")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)
    ax.axvspan(tt[200], tt[-1], color="gray", alpha=0.08)
ax = axes[2]
for i, lab in enumerate(["x", "y", "z"]):
    ax.plot(tt, y[i], lw=1.0, label=f"SMPC {lab}")
    ax.plot(tt, y_b[i], lw=0.8, ls="--", alpha=0.7, label=f"open {lab}")
ax.set_ylabel("position (m)")
ax.legend(loc="upper right", fontsize=8, ncol=2)
ax.grid(alpha=0.3)
ax.axvspan(tt[200], tt[-1], color="gray", alpha=0.08)
ax = axes[3]
for i in range(4):
    ax.plot(tt, u[i], lw=0.9, label=f"RPM {i+1}")
ax.set_ylabel("motor RPM cmd")
ax.legend(loc="upper right", fontsize=8, ncol=4)
ax.grid(alpha=0.3)
ax.set_xlabel("time (s)")
ax.axvspan(tt[200], tt[-1], color="gray", alpha=0.08)
axes[0].set_title(
    f"copyAZ closed loop on CF2X (gym-pybullet-drones): SMPC vs open-loop hover\n"
    f"task-axis RMSE SMPC={np.round(np.sqrt(np.mean((s[:,warm]-ref[:,warm])**2,1)),3)}, "
    f"baseline={np.round(np.sqrt(np.mean((s_b[:,warm]-ref[:,warm])**2,1)),3)}; "
    f"QP success={d['qp_success'][0,0]:.3f}", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "results", "copyAZ_pybullet_closed_loop_fig.png"), dpi=130)
print("saved copyAZ_pybullet_closed_loop_fig.png")
