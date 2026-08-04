"""One-shot jitter diagnosis for copyBG results (experiment layer only)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat

BASE = Path(__file__).resolve().parent
DIRS = {
    "growing": BASE / "results_growing",
    "boundary": BASE / "results",
    "smooth_dc": BASE / "results_smooth_dc_feasible",
}


def _arr(m, *names):
    keys = [k for k in m.keys() if not k.startswith("_")]
    for name in names:
        if name in m:
            a = np.asarray(m[name], dtype=float)
            if a.ndim == 2 and a.shape[0] < a.shape[1] and a.shape[0] <= 8:
                a = a.T
            return a
    for name in names:
        for k in keys:
            if name.lower() in k.lower():
                a = np.asarray(m[k], dtype=float)
                if a.ndim == 2 and a.shape[0] < a.shape[1] and a.shape[0] <= 8:
                    a = a.T
                return a
    return None


def analyze(name: str, d: Path) -> None:
    jp = next(d.glob("*.json"), None)
    mp = next(d.glob("*.mat"), None)
    print(f"\n=== {name} ({d.name}) ===")
    if jp is None or mp is None:
        print("missing artifacts")
        return
    s = json.loads(jp.read_text(encoding="utf-8"))
    m = loadmat(mp, squeeze_me=True, struct_as_record=False)
    keys = [k for k in m.keys() if not k.startswith("_")]
    print("mat keys:", keys)
    print(
        "mode", s.get("reference_mode"),
        "steps", s.get("requested_steps"),
        "sigma_wind", s.get("sigma_wind_mps"),
    )
    sm = s.get("smpc", {})
    print(
        "RMSE", sm.get("rmse"),
        "active_qp", sm.get("active_qp_steps"),
        "sat", sm.get("input_saturation_rate"),
        "min_hard_margin", sm.get("minimum_hard_margin"),
    )
    print("smpc_vs_mpc", s.get("smpc_vs_mpc"))
    gm = s.get("growing_smooth_square")
    if gm:
        print("growth", gm)

    ref = _arr(m, "reference", "ref_task", "Y_ref", "reference_task")
    y = _arr(m, "y_task_smpc", "trajectory_smpc", "task_smpc", "Y_smpc")
    u = _arr(m, "u_smpc", "control_smpc", "U_smpc", "inputs_smpc")
    if ref is None:
        print("no reference in mat")
        return
    speed = np.linalg.norm(np.diff(ref, axis=0), axis=1)
    d2 = np.linalg.norm(np.diff(ref, n=2, axis=0), axis=1)
    print(
        f"ref {ref.shape} speed mean/p95/max "
        f"{speed.mean():.5f}/{np.percentile(speed, 95):.5f}/{speed.max():.5f} "
        f"std={speed.std():.5f}"
    )
    print(
        f"ref |Δ²| mean/p95/max "
        f"{d2.mean():.6f}/{np.percentile(d2, 95):.6f}/{d2.max():.6f}"
    )
    # corner vs side: for rounded square, |x|~|y| near corners
    rbox = np.maximum(np.abs(ref[:, 0]), np.abs(ref[:, 1]))
    thr = 0.12 * np.maximum(rbox, 1e-6)
    near_c = (np.abs(np.abs(ref[:, 0]) - np.abs(ref[:, 1])) < thr) & (
        rbox > 0.55 * np.max(rbox)
    )
    side = ~near_c & (rbox > 0.4 * np.max(rbox))
    if near_c[:-2].any() and side[:-2].any():
        print(
            "ref |Δ²| near-corner vs side:",
            f"{d2[near_c[:-2]].mean():.6f} vs {d2[side[:-2]].mean():.6f}",
        )

    if y is not None:
        n = min(len(y), len(ref))
        e = y[:n] - ref[:n]
        de = np.linalg.norm(np.diff(e, axis=0), axis=1)
        dy = np.linalg.norm(np.diff(y[:n], axis=0), axis=1)
        print(
            "track RMS", np.sqrt((e ** 2).mean(0)),
            f"|Δe| mean/p95 {de.mean():.5f}/{np.percentile(de, 95):.5f}",
            f"|Δy| p95 {np.percentile(dy, 95):.5f}",
        )
        # residual after removing slow trend (high-pass via diff)
        print(
            "err high-freq proxy: std(Δe_axis)=",
            np.std(np.diff(e, axis=0), axis=0),
        )
        if near_c[: n - 1].any():
            print(
                "|Δe| corner vs other:",
                f"{de[near_c[: n - 1]].mean():.5f} vs "
                f"{de[~near_c[: n - 1]].mean():.5f}",
            )

    if u is not None:
        du = np.max(np.abs(np.diff(u, axis=0)), axis=1)
        hold = float((du < 1e-12).mean())
        print(
            f"u {u.shape} hold_frac={hold:.3f} "
            f"|Δu|_∞ mean/p95 {du.mean():.4f}/{np.percentile(du, 95):.4f} "
            f"sat={(np.max(np.abs(u), axis=1) >= 5.999).mean():.3f} "
            f"u∞ p95={np.percentile(np.max(np.abs(u), axis=1), 95):.3f}"
        )
        # control updates only every d=10 → nonzero du every 10 steps
        nz = np.where(du > 1e-12)[0]
        if len(nz) > 5:
            gaps = np.diff(nz)
            print(
                f"control update gap median/mean {np.median(gaps):.1f}/"
                f"{gaps.mean():.2f} (expect ~10)"
            )


if __name__ == "__main__":
    for name, d in DIRS.items():
        analyze(name, d)
