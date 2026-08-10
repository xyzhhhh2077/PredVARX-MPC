#!/usr/bin/env python3
"""
Generate a scientific, publication-ready comparison figure for the
20-seed fair benchmark (BJ vs BL vs Frozen).

Design philosophy
-----------------
* OLD figure (aggregate_comparison.png): 3 subplots side-by-side mixing
  unequal sample sizes (Frozen n=10, BJ/BL n=16), no direction annotation,
  no paired structure.  Reviewers / supervisors could not read it.

* NEW figure: built around the BJ-BL strictly-paired difference (n=15),
  with per-seed paired dot-difference plots, violin+strip density overlays,
  paired t-test p-values and Cohen's d, explicit direction (↓/↑ better),
  and stratified sample-size labels.

Reads (read-only):
  results/per_seed.json
  results/aggregate.json
Writes:
  results/scientific_comparison.png
"""

import json
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy import stats

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"
PER_SEED = RESULTS_DIR / "per_seed.json"
AGGREGATE = RESULTS_DIR / "aggregate.json"
OUT_PNG = RESULTS_DIR / "scientific_comparison.png"

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
with open(PER_SEED) as f:
    per_seed = json.load(f)
with open(AGGREGATE) as f:
    aggregate = json.load(f)

seed_records = per_seed["seeds"]

# ---------------------------------------------------------------------------
# Extract per-seed metric values; only complete runs contribute
# ---------------------------------------------------------------------------
CONTROLLERS = ["Frozen", "BJ", "BL"]


def get_metric(ctrl_block, key):
    """Return the metric value or None if missing/incomplete."""
    if not ctrl_block.get("complete", False):
        return None
    if key == "rmse_xyz":
        rmse = ctrl_block.get("rmse_xyz")
        if rmse is None:
            return None
        # Euclidean RMSE across (x, y, z)
        return float(np.sqrt(np.sum(np.array(rmse) ** 2)))
    if key == "rmse_x":
        rmse = ctrl_block.get("rmse_xyz")
        return float(rmse[0]) if rmse else None
    if key == "rmse_y":
        rmse = ctrl_block.get("rmse_xyz")
        return float(rmse[1]) if rmse else None
    if key == "rmse_z":
        rmse = ctrl_block.get("rmse_xyz")
        return float(rmse[2]) if rmse else None
    val = ctrl_block.get(key)
    return float(val) if val is not None else None


# Collect per-seed values into dicts keyed by seed index
def collect(key):
    out = {}
    for idx, rec in enumerate(seed_records):
        ctrls = rec["controllers"]
        row = {}
        for c in CONTROLLERS:
            v = get_metric(ctrls[c], key)
            if v is not None:
                row[c] = v
        if row:
            out[idx] = row
    return out

metrics_raw = {
    "rmse_xyz": collect("rmse_xyz"),
    "cum_cost": collect("cumulative_realized_stage_cost"),
    "min_margin": collect("minimum_hard_margin"),
}

# Sample sizes (complete runs per controller per metric)
def n_complete(key):
    return {c: sum(1 for s in metrics_raw[key].values() if c in s) for c in CONTROLLERS}

# ---------------------------------------------------------------------------
# Paired BJ-BL selection (both complete, same seed)
# ---------------------------------------------------------------------------
def paired_indices(key):
    """Return list of (seed_idx, bj_val, bl_val) where BOTH BJ and BL are complete."""
    pairs = []
    for idx, row in metrics_raw[key].items():
        if "BJ" in row and "BL" in row:
            pairs.append((idx, row["BJ"], row["BL"]))
    return pairs

paired_rmse = paired_indices("rmse_xyz")
paired_cost = paired_indices("cum_cost")
paired_margin = paired_indices("min_margin")

# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------
def paired_stats(bj_vals, bl_vals):
    diff = np.array(bl_vals) - np.array(bj_vals)   # BL - BJ  (positive => BL worse for cost/rmse)
    mean_diff = np.mean(diff)
    std_diff = np.std(diff, ddof=1)
    n = len(diff)
    se = std_diff / np.sqrt(n)
    ci95 = 1.96 * se
    t_stat, p_val = stats.ttest_rel(np.array(bl_vals), np.array(bj_vals))
    # Cohen's d for paired: mean_diff / std_diff
    cohen_d = mean_diff / std_diff if std_diff > 0 else 0.0
    return {
        "n": n,
        "mean_diff": mean_diff,
        "std_diff": std_diff,
        "ci95": ci95,
        "t": t_stat,
        "p": p_val,
        "d": cohen_d,
    }

def fmt_p(p):
    if p < 0.001:
        return "p<0.001"
    elif p < 0.01:
        return f"p={p:.3f}"
    elif p < 0.05:
        return f"p={p:.3f}"
    else:
        return f"p={p:.2f}"

# ---------------------------------------------------------------------------
# Figure layout: 2 rows x 3 columns
#   Row 1: paired dot-difference (left=BJ, right=BL, connecting lines) per metric
#   Row 2: violin + strip overlay of paired difference (BL-BJ) per metric
# ---------------------------------------------------------------------------
sns.set_theme(style="whitegrid", palette="colorblind", font_scale=1.0)
plt.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 200,
    "font.family": "DejaVu Sans",
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.titleweight": "bold",
})

palettes = sns.color_palette("colorblind")
C_BJ = palettes[0]      # blue
C_BL = palettes[1]      # orange
C_FROZEN = palettes[7]  # gray
C_DIFF = palettes[2]     # green

fig = plt.figure(figsize=(16, 9))
gs = gridspec.GridSpec(
    2, 3, figure=fig,
    hspace=0.42, wspace=0.32,
    left=0.07, right=0.97, top=0.92, bottom=0.10,
)

# Metric display config
# direction: "lower_better" means ↓ better; "higher_better" means ↑ better
metric_config = [
    {
        "key": "rmse_xyz",
        "paired": paired_rmse,
        "title": "Tracking RMSE (‖x,y,z‖)",
        "ylabel": "RMSE [m]",
        "direction": "↓ better",
        "direction_xy": (-0.15, 0.98),
    },
    {
        "key": "cum_cost",
        "paired": paired_cost,
        "title": "Cumulative Realized Stage Cost",
        "ylabel": "Cost",
        "direction": "↓ better",
        "direction_xy": (-0.15, 0.98),
    },
    {
        "key": "min_margin",
        "paired": paired_margin,
        "title": "Minimum Hard Safety Margin",
        "ylabel": "Margin [m]",
        "direction": "↑ better",
        "direction_xy": (1.15, 0.98),
    },
]

# Sample-size annotation string
n_rmse = n_complete("rmse_xyz")
n_cost = n_complete("cum_cost")
n_margin = n_complete("min_margin")
n_paired = {m["key"]: len(m["paired"]) for m in metric_config}

# ---- Row 1: Paired dot-difference plot ----
for col, mc in enumerate(metric_config):
    ax = fig.add_subplot(gs[0, col])
    pairs = mc["paired"]
    bj_vals = [p[1] for p in pairs]
    bl_vals = [p[2] for p in pairs]
    n = len(pairs)

    # Sort by BJ value for visual clarity
    order = np.argsort(bj_vals)
    bj_sorted = np.array(bj_vals)[order]
    bl_sorted = np.array(bl_vals)[order]

    x_left = 0.0
    x_right = 1.0

    # Connect lines
    for i in range(n):
        color = C_BL if bl_sorted[i] < bj_sorted[i] else C_BJ
        alpha = 0.45 if bl_sorted[i] < bj_sorted[i] else 0.25
        lw = 0.9 if abs(bj_sorted[i]-bl_sorted[i]) > 1e-6 else 0.4
        ax.plot([x_left, x_right], [bj_sorted[i], bl_sorted[i]],
                color=color, alpha=alpha, linewidth=lw, zorder=1)

    # BJ dots
    jit_l = np.random.RandomState(42).uniform(-0.03, 0.03, n)
    ax.scatter(np.full(n, x_left) + jit_l, bj_sorted,
               s=55, color=C_BJ, edgecolor="white", linewidth=0.6, zorder=3,
               label=f"BJ (n={n_complete(mc['key'])['BJ']})")
    # BL dots
    jit_r = np.random.RandomState(43).uniform(-0.03, 0.03, n)
    ax.scatter(np.full(n, x_right) + jit_r, bl_sorted,
               s=55, color=C_BL, edgecolor="white", linewidth=0.6, zorder=3,
               label=f"BL (n={n_complete(mc['key'])['BL']})")

    # Mean markers
    ax.scatter([x_left], [np.mean(bj_sorted)], marker="D", s=100, color=C_BJ,
               edgecolor="black", linewidth=1.2, zorder=4, label="BJ mean")
    ax.scatter([x_right], [np.mean(bl_sorted)], marker="D", s=100, color=C_BL,
               edgecolor="black", linewidth=1.2, zorder=4, label="BL mean")

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["BJ", "BL"])
    ax.set_ylabel(mc["ylabel"])
    ax.set_title(mc["title"])
    # Direction annotation
    ax.annotate(mc["direction"], xy=mc["direction_xy"], xycoords="axes fraction",
                fontsize=10, fontweight="bold", color="0.3",
                ha="left" if mc["direction_xy"][0] < 0 else "right",
                va="top",
                bbox=dict(boxstyle="round,pad=0.2", fc="lightyellow", ec="0.7", lw=0.5))
    ax.legend(fontsize=8, loc="best", framealpha=0.8)

    # Paired n annotation
    ax.text(0.5, -0.14, f"paired n = {n}", transform=ax.transAxes,
            ha="center", fontsize=9, style="italic", color="0.4")

# ---- Row 2: Paired difference (BL − BJ) violin + strip ----
for col, mc in enumerate(metric_config):
    ax = fig.add_subplot(gs[1, col])
    pairs = mc["paired"]
    bj_vals = np.array([p[1] for p in pairs])
    bl_vals = np.array([p[2] for p in pairs])
    diff = bl_vals - bj_vals  # BL − BJ
    st = paired_stats([p[1] for p in pairs], [p[2] for p in pairs])

    # Determine if lower diff means BL is "better" for this metric
    # For rmse & cost: negative diff => BL lower => better; For margin: positive diff => BL higher => better
    better_is_negative = mc["direction"] == "↓ better"

    # Violin (single-sided, via sns.kdeplot fill)
    # Use strip + violin
    data_arr = diff
    # KDE plot
    try:
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(data_arr, bw_method=0.4)
        x_range = np.linspace(np.min(data_arr) - abs(np.mean(data_arr))*0.5,
                              np.max(data_arr) + abs(np.mean(data_arr))*0.5, 200)
        density = kde(x_range)
        # Normalize violin width
        max_d = np.max(density) if np.max(density) > 0 else 1.0
        ax.fill_betweenx(x_range, 0, density / max_d * 0.7,
                         color=C_DIFF, alpha=0.25, zorder=1)
        ax.plot(density / max_d * 0.7, x_range, color=C_DIFF, linewidth=1.5, zorder=2)
    except Exception:
        ax.hist(diff, bins=8, orientation="horizontal", color=C_DIFF, alpha=0.3)

    # Strip (individual seed points) with jitter
    jit = np.random.RandomState(99).uniform(-0.05, 0.05, len(diff))
    colors_pts = [C_BL if d < 0 else C_BJ for d in diff] if better_is_negative else \
                  [C_BL if d > 0 else C_BJ for d in diff]
    ax.scatter(np.full(len(diff), 0.0) + jit + 0.42,
               diff, s=40, c=colors_pts, edgecolor="white", linewidth=0.5,
               zorder=4, alpha=0.85)

    # Mean line
    ax.axhline(st["mean_diff"], color="0.2", linewidth=1.8, linestyle="-", zorder=3)
    # CI band
    ax.axhspan(st["mean_diff"] - st["ci95"], st["mean_diff"] + st["ci95"],
               color="0.5", alpha=0.18)

    # Zero reference line
    ax.axhline(0, color="red", linewidth=1.0, linestyle="--", zorder=2)

    # Stats annotation box
    stat_txt = (
        f"Δμ = {st['mean_diff']:.4g}\n"
        f"95% CI = [±{st['ci95']:.4g}]\n"
        f"t = {st['t']:.3f},  {fmt_p(st['p'])}\n"
        f"Cohen's d = {st['d']:.3f}"
    )
    x_box = 0.97 if st["mean_diff"] >= 0 else 0.03
    ha = "right" if st["mean_diff"] >= 0 else "left"
    ax.text(x_box, -0.08, stat_txt, transform=ax.transAxes,
            ha=ha, va="bottom", fontsize=9, fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=C_DIFF, lw=1.2))

    ax.set_xlabel("BL − BJ  (per seed)")
    ax.set_xlim(-0.45, 1.05)
    ax.set_xticks([])
    ax.set_title(f"Paired Difference: {mc['title']}")
    ax.annotate(mc["direction"], xy=(0.02, 0.98), xycoords="axes fraction",
                fontsize=9, fontweight="bold", color="0.3", va="top",
                bbox=dict(boxstyle="round,pad=0.2", fc="lightyellow", ec="0.7", lw=0.5))

# ---- Supertitle + footer ----
fig.suptitle(
    f"Simulation Copilot — 20-Seed Automated Fair Benchmark: BJ vs. BL vs. Frozen\n"
    f"Same-disturbance paired comparison — Frozen n={n_rmse['Frozen']}, BJ n={n_rmse['BJ']}, BL n={n_rmse['BL']}, "
    f"BJ–BL paired n={n_paired['rmse_xyz']}",
    fontsize=15, fontweight="bold", y=0.98
)

footer = (
    f"Scope: {aggregate.get('scope','')}  |  steps={aggregate.get('steps','')}  "
    f"|  config_sha256={aggregate.get('config_sha256','')[:12]}…  "
    f"|  Paired stats: paired-t (two-sided), Cohen's d (paired)"
)
fig.text(0.5, 0.012, footer, ha="center", va="bottom", fontsize=8, color="0.35",
         fontstyle="italic")

# Colors legend key (small inline in top-right via text)
fig.text(0.99, 0.956,
         "Blue dot = BJ    Orange dot = BL    Green = BL−BJ difference density",
         ha="right", va="top", fontsize=8, color="0.4")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight", facecolor="white")
print(f"Saved: {OUT_PNG}")
print(f"File size: {os.path.getsize(OUT_PNG):,} bytes")

# Print key numbers for audit
print("\n=== Paired statistics (BJ vs BL) ===")
for mc in metric_config:
    pairs = mc["paired"]
    st = paired_stats([p[1] for p in pairs], [p[2] for p in pairs])
    print(f"\n[{mc['key']}]  paired n = {st['n']}")
    print(f"  BJ mean = {np.mean([p[1] for p in pairs]):.4g}   "
          f"BL mean = {np.mean([p[2] for p in pairs]):.4g}")
    print(f"  Δμ(BL−BJ) = {st['mean_diff']:.4g}   CI95 ±{st['ci95']:.4g}")
    print(f"  t = {st['t']:.4f},  p = {st['p']:.4f},  Cohen's d = {st['d']:.4f}")
