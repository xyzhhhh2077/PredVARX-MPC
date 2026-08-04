"""Build the face-petal-stress comparison figure + REPORT.md."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent


def main():
    j = json.loads((OUT / "copyBH_pelican_u_cap_ablation.json").read_text(encoding="utf-8"))
    z = np.load(OUT / "copyBH_pelican_u_cap_ablation.npz")
    ref = z["reference_standardized"]
    smpc = z["smpc_trajectory_standardized"]
    mpc = z["deterministic_mpc_trajectory_standardized"]
    act = np.asarray(z["smpc_active_history"]).astype(bool).ravel()
    hard = float(np.asarray(z["hard_bound_standardized"], dtype=float).reshape(-1)[0])
    mean_bound = np.asarray(z["stage_mean_bound_standardized"], float)
    term = mean_bound[-1]
    t = np.arange(len(ref)) * float(np.asarray(z["dt_used_seconds"], dtype=float).reshape(-1)[0])

    fig = plt.figure(figsize=(12.5, 8.2))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.35, 1.0, 0.7], hspace=0.35, wspace=0.25)

    ax3 = fig.add_subplot(gs[0, 0], projection="3d")
    ax3.plot(ref[:, 0], ref[:, 1], ref[:, 2], "k--", lw=1.0, alpha=0.75, label="reference")
    ax3.plot(smpc[:, 0], smpc[:, 1], smpc[:, 2], color="#1f77b4", lw=1.1, label="SMPC")
    ax3.plot(mpc[:, 0], mpc[:, 1], mpc[:, 2], color="#d62728", lw=0.9, alpha=0.85, label="det. MPC")
    ax3.scatter(*ref[0], c="g", s=30, depthshade=False)
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    ax3.set_zlabel("z")
    ax3.set_title("3D face-petal stress")
    ax3.legend(loc="upper left", fontsize=8)

    axxy = fig.add_subplot(gs[0, 1])
    axxy.plot(ref[:, 0], ref[:, 1], "k--", lw=1.0, alpha=0.7, label="ref")
    axxy.plot(smpc[:, 0], smpc[:, 1], color="#1f77b4", lw=1.0, label="SMPC")
    axxy.plot(mpc[:, 0], mpc[:, 1], color="#d62728", lw=0.9, alpha=0.85, label="MPC")
    hs = hard
    axxy.plot(
        [-hs, hs, hs, -hs, -hs],
        [-hs, -hs, hs, hs, -hs],
        color="0.3", ls=":", lw=1.0, label="hard box",
    )
    axxy.set_aspect("equal", adjustable="box")
    axxy.set_xlabel("x")
    axxy.set_ylabel("y")
    axxy.set_title("top view (x-y)")
    axxy.legend(fontsize=8, loc="upper right")
    axxy.grid(True, alpha=0.25)

    ax = fig.add_subplot(gs[1, :])
    cols = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    for i, name in enumerate("xyz"):
        ax.plot(t, ref[:, i], color=cols[i], ls="--", lw=0.9, alpha=0.55)
        ax.plot(t, smpc[:, i], color=cols[i], lw=1.15, label=f"{name} SMPC")
        ax.plot(t, mpc[:, i], color=cols[i], ls=":", lw=1.0, alpha=0.9)
    ax.axhline(hard, color="k", ls="--", lw=1.0, label="hard")
    ax.axhline(-hard, color="k", ls="--", lw=1.0)
    for i, tb in enumerate(term):
        ax.axhline(tb, color=cols[i], ls=":", lw=0.8, alpha=0.5)
        ax.axhline(-tb, color=cols[i], ls=":", lw=0.8, alpha=0.5)
    ax.set_xlim(t[0], t[-1])
    ax.set_ylabel("task (std)")
    ax.set_title(
        "task vs time (solid=SMPC, dotted=MPC, dashed=ref; pale dotted=terminal mean bound)"
    )
    ax.legend(ncol=5, fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.25)

    axb = fig.add_subplot(gs[2, :])
    if act.size == len(ref):
        act_s = act
        tt = t
    elif act.size and len(ref) % act.size == 0:
        rep = len(ref) // act.size
        act_s = np.repeat(act, rep)[: len(ref)]
        tt = t
    else:
        act_s = act.astype(float)
        tt = np.linspace(t[0], t[-1], act.size, endpoint=True)
    axb.fill_between(
        tt, 0, act_s.astype(float), step="pre", color="#1f77b4", alpha=0.35,
        label="SMPC chance active",
    )
    err_s = np.linalg.norm(smpc - ref, axis=1)
    err_m = np.linalg.norm(mpc - ref, axis=1)
    axb2 = axb.twinx()
    axb2.plot(t, err_s, color="#1f77b4", lw=1.0, label="||e|| SMPC")
    axb2.plot(t, err_m, color="#d62728", lw=0.9, alpha=0.85, label="||e|| MPC")
    axb.set_xlabel("t [s]")
    axb.set_ylabel("chance active")
    axb2.set_ylabel("||task-ref||")
    axb.set_xlim(t[0], t[-1])
    axb.set_ylim(-0.05, 1.15)
    axb.set_title("chance-QP activation (SMPC) and tracking error")
    h1, l1 = axb.get_legend_handles_labels()
    h2, l2 = axb2.get_legend_handles_labels()
    axb.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
    axb.grid(True, alpha=0.25)

    fig.suptitle(
        (
            "copyBH face-petal-stress | cap=+/-7 | hard={:.3f} | "
            "mean_bound~[{:.2f},{:.2f},{:.2f}] | "
            "SMPC activeQP={}/{}, min hard margin {:.3f} vs MPC {:.3f}"
        ).format(
            hard, term[0], term[1], term[2],
            j["smpc"]["active_qp_steps"], j["smpc"]["qp_count"],
            j["smpc"]["minimum_hard_margin_standardized"],
            j["deterministic_mpc"]["minimum_hard_margin_standardized"],
        ),
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig_path = OUT / "copyBH_face_petal_stress_comparison.png"
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)
    print("wrote", fig_path)

    lines = []
    lines.append("# copyBH face-petal-stress REPORT")
    lines.append("")
    lines.append("Reference-only multi-face stress test on frozen copyBG controller.")
    lines.append("")
    lines.append("## Setup")
    lines.append("- reference: face-petal-stress (center -> face -> on-face slide -> center x6)")
    lines.append(f"- pressures: {j.get('reference_pressures_standardized')}")
    lines.append(f"- hard: {hard:.6f}")
    lines.append(f"- terminal chance mean bound: {term.tolist()}")
    lines.append(f"- input cap: {j.get('input_bound_normalized')}")
    lines.append(
        f"- steps: {j.get('requested_steps')} @ 100 Hz, wind sigma={j.get('sigma_wind_mps')} m/s"
    )
    lines.append(
        f"- seeds: wind={j.get('wind_seed')}, innovation={j.get('innovation_seed')}"
    )
    lines.append("")
    lines.append("## Results (SMPC vs deterministic MPC, shared noise/wind)")
    lines.append("")
    lines.append("| metric | SMPC | det. MPC |")
    lines.append("|---|---:|---:|")
    for label, key in [
        ("completed_steps", "completed_steps"),
        ("hard_violation_steps", "hard_violation_steps"),
        ("hard_violation_rate", "hard_violation_rate"),
        ("active_qp_steps", "active_qp_steps"),
        ("qp_failure_count", "qp_failure_count"),
        ("fallback_count", "fallback_count"),
        ("input_saturation_steps", "input_saturation_steps"),
        ("min_hard_margin", "minimum_hard_margin_standardized"),
    ]:
        lines.append(
            f"| {label} | {j['smpc'][key]} | {j['deterministic_mpc'][key]} |"
        )
    lines.append(
        "| peak |x|,|y|,|z| | "
        f"{np.round(j['smpc']['peak_absolute_position_standardized'], 3).tolist()} | "
        f"{np.round(j['deterministic_mpc']['peak_absolute_position_standardized'], 3).tolist()} |"
    )
    lines.append(
        f"| RMSE xyz | {np.round(j['smpc']['rmse_standardized'], 4).tolist()} | "
        f"{np.round(j['deterministic_mpc']['rmse_standardized'], 4).tolist()} |"
    )
    lines.append("")
    lines.append("## Readout")
    lines.append("- Both controllers finish 18000 steps with **0 hard violations**.")
    lines.append(
        f"- SMPC chance QP active on **{j['smpc']['active_qp_steps']}** decision steps; "
        "MPC active=0 (no chance layer)."
    )
    lines.append(
        f"- SMPC keeps a larger minimum hard margin "
        f"(**{j['smpc']['minimum_hard_margin_standardized']:.4f}** vs "
        f"**{j['deterministic_mpc']['minimum_hard_margin_standardized']:.4f}**), mainly on x."
    )
    lines.append(
        f"- Peak |x|: SMPC **{j['smpc']['peak_absolute_position_standardized'][0]:.3f}** < "
        f"MPC **{j['deterministic_mpc']['peak_absolute_position_standardized'][0]:.3f}** "
        "(SMPC pulls back under risk constraints)."
    )
    lines.append(
        "- This is a **directional stress demo**, not a sustained near-hard cruise "
        "(that remains v4.0 z-face)."
    )
    lines.append("")
    lines.append("## Artifacts")
    lines.append("- `copyBH_face_petal_stress_comparison.png`")
    lines.append("- `copyBH_pelican_u_cap_ablation_smpc.png`")
    lines.append("- `copyBH_pelican_u_cap_ablation_deterministic_mpc.png`")
    lines.append("- `copyBH_pelican_u_cap_ablation.{json,npz,mat}`")
    lines.append("- `preview_reference.png`")
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("wrote REPORT.md")
    print(
        "active qp fraction of decisions",
        j["smpc"]["active_qp_steps"] / j["smpc"]["qp_count"],
    )


if __name__ == "__main__":
    main()
