#!/usr/bin/env python
"""Compare fixed-seed copyBH input-cap sweep artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


STEM = "copyBH_pelican_u_cap_ablation"
CONTROLLERS = {
    "smpc": ("smpc_trajectory_standardized", "smpc_control"),
    "deterministic_mpc": (
        "deterministic_mpc_trajectory_standardized",
        "deterministic_mpc_control",
    ),
}
FAIRNESS_ARRAYS = (
    "reference_standardized",
    "initial_state",
    "initial_control",
    "innovation_latent",
    "wind_velocity_mps",
    "stage_tightening_standardized",
    "stage_mean_bound_standardized",
)


def _array_hash(array):
    values = np.ascontiguousarray(np.asarray(array))
    return hashlib.sha256(values.view(np.uint8)).hexdigest()


def _longest_true_run(mask):
    mask = np.asarray(mask, dtype=bool)
    if not np.any(mask):
        return 0
    padded = np.r_[False, mask, False].astype(np.int8)
    edges = np.flatnonzero(np.diff(padded))
    return int(np.max(edges[1::2] - edges[::2]))


def _load_run(label, directory):
    directory = Path(directory)
    with (directory / f"{STEM}.json").open(encoding="utf-8") as stream:
        summary = json.load(stream)
    archive = np.load(directory / f"{STEM}.npz")
    return {
        "label": str(label),
        "directory": directory,
        "summary": summary,
        "archive": archive,
    }


def _controller_metrics(archive, task_key, input_key, cap, threshold, hard_bound):
    task = np.asarray(archive[task_key], dtype=float)
    command = np.asarray(archive[input_key], dtype=float)
    reference = np.asarray(archive["reference_standardized"], dtype=float)[: len(task)]
    metrics = {
        "max_abs_task": np.max(np.abs(task), axis=0).tolist(),
        "max_u_inf": float(np.max(np.abs(command))),
        "input_saturation_steps": int(
            np.count_nonzero(np.any(np.abs(command) >= cap - 1e-4, axis=1))
        ),
        "input_saturation_fraction": float(
            np.mean(np.any(np.abs(command) >= cap - 1e-4, axis=1))
        ),
        "hard_violation_steps": int(
            np.count_nonzero(np.any(np.abs(task) > hard_bound + 1e-9, axis=1))
        ),
        "axes": {},
    }
    for axis, name in enumerate(("x", "y")):
        target_mask = reference[:, axis] >= threshold
        attained = task[:, axis] >= threshold
        target_count = int(np.count_nonzero(target_mask))
        attained_on_target = attained & target_mask
        metrics["axes"][name] = {
            "target_near_hard_steps": target_count,
            "attained_near_hard_steps": int(np.count_nonzero(attained_on_target)),
            "attained_fraction_during_target": (
                float(np.mean(attained[target_mask])) if target_count else 0.0
            ),
            "longest_attained_seconds": float(
                _longest_true_run(attained_on_target) * 0.01
            ),
            "mean_task_during_target": (
                float(np.mean(task[target_mask, axis])) if target_count else None
            ),
            "rmse_during_target": (
                float(np.sqrt(np.mean(
                    (task[target_mask, axis] - reference[target_mask, axis]) ** 2
                ))) if target_count else None
            ),
        }
    return metrics


def analyze(runs, output_dir):
    loaded = [_load_run(label, directory) for label, directory in runs]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fairness = {}
    for key in FAIRNESS_ARRAYS:
        hashes = {run["label"]: _array_hash(run["archive"][key]) for run in loaded}
        fairness[key] = {"identical": len(set(hashes.values())) == 1, "hashes": hashes}
    if not all(item["identical"] for item in fairness.values()):
        raise RuntimeError("Fairness gate failed: sweep inputs differ")

    first = loaded[0]
    hard_bound = float(first["summary"]["hard_bound_standardized"])
    threshold = hard_bound - 0.4
    terminal_mean_bound = np.asarray(
        first["archive"]["stage_mean_bound_standardized"], dtype=float
    )[-1]
    comparison = {
        "experiment": STEM,
        "scope": "input-cap-only ablation; controller and theory unchanged",
        "hard_bound": hard_bound,
        "near_hard_threshold": threshold,
        "terminal_mean_bound_abs": terminal_mean_bound.tolist(),
        "fairness": fairness,
        "runs": {},
    }
    for run in loaded:
        cap = float(run["summary"]["input_bound_normalized"][1])
        run_metrics = {
            "cap": cap,
            "reference_path": run["summary"]["reference_path"],
            "controllers": {},
        }
        for name, (task_key, input_key) in CONTROLLERS.items():
            controller_metrics = _controller_metrics(
                run["archive"], task_key, input_key, cap, threshold, hard_bound
            )
            controller_summary = run["summary"][name]
            controller_metrics["run_status"] = {
                key: controller_summary[key]
                for key in (
                    "completed_steps",
                    "qp_failure_count",
                    "qp_failure_step",
                    "qp_failure_status",
                    "fallback_count",
                    "active_qp_steps",
                    "maximum_qp_constraint_residual",
                    "minimum_hard_margin_standardized",
                )
            }
            run_metrics["controllers"][name] = controller_metrics
        comparison["runs"][run["label"]] = run_metrics

    json_path = output_dir / "copyBH_u_cap_sweep_comparison.json"
    json_path.write_text(json.dumps(comparison, indent=2), encoding="utf-8")

    figure, axes = plt.subplots(2, 2, figsize=(14, 8.5), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.12, 0.88, len(loaded)))
    reference = np.asarray(first["archive"]["reference_standardized"], dtype=float)
    time = np.arange(len(reference)) * 0.01
    for axis, name in enumerate(("x", "y")):
        panel = axes[axis, 0]
        panel.plot(time, reference[:, axis], "k--", lw=1.5, label="reference")
        for color, run in zip(colors, loaded):
            task = np.asarray(
                run["archive"]["smpc_trajectory_standardized"], dtype=float
            )
            panel.plot(time[:len(task)], task[:, axis], color=color, lw=1.0,
                       label=f"SMPC cap ±{run['label']}")
        panel.axhline(threshold, color="#d62728", ls=":", lw=1.4,
                      label="near-hard threshold" if axis == 0 else None)
        panel.axhline(terminal_mean_bound[axis], color="#ff7f0e", ls="-.", lw=1.2,
                      label="H=18 mean bound" if axis == 0 else None)
        panel.axhline(hard_bound, color="black", ls="-", lw=0.8, alpha=0.6,
                      label="hard bound" if axis == 0 else None)
        panel.set(title=f"{name}-axis sustained-face trial", xlabel="time [s]",
                  ylabel=f"task {name} [standardized]")
        panel.grid(alpha=0.25)
        panel.legend(fontsize=8, loc="best")

    labels = [run["label"] for run in loaded]
    x = np.arange(len(labels), dtype=float)
    width = 0.18
    panel = axes[0, 1]
    series = [
        ("SMPC x", "smpc", "x", -1.5 * width, "#1f77b4"),
        ("SMPC y", "smpc", "y", -0.5 * width, "#2ca02c"),
        ("MPC x", "deterministic_mpc", "x", 0.5 * width, "#9467bd"),
        ("MPC y", "deterministic_mpc", "y", 1.5 * width, "#8c564b"),
    ]
    for title, controller, axis_name, offset, color in series:
        values = [
            100.0 * comparison["runs"][label]["controllers"][controller]
            ["axes"][axis_name]["attained_fraction_during_target"]
            for label in labels
        ]
        panel.bar(x + offset, values, width, label=title, color=color)
    panel.set_xticks(x, [f"±{label}" for label in labels])
    panel.set(title="Near-hard attainment during face target",
              xlabel="input cap", ylabel="target-window samples [%]", ylim=(0, 105))
    panel.grid(axis="y", alpha=0.25)
    panel.legend(fontsize=8, ncol=2)

    panel = axes[1, 1]
    for j, (controller, color) in enumerate((
        ("smpc", "#1f77b4"), ("deterministic_mpc", "#9467bd")
    )):
        values = [
            100.0 * comparison["runs"][label]["controllers"][controller]
            ["input_saturation_fraction"] for label in labels
        ]
        panel.bar(x + (j - 0.5) * 0.32, values, 0.32,
                  label=controller.replace("_", " "), color=color)
    panel.set_xticks(x, [f"±{label}" for label in labels])
    panel.set(title="Input saturation", xlabel="input cap",
              ylabel="steps at selected cap [%]")
    panel.grid(axis="y", alpha=0.25)
    panel.legend(
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=2,
    )
    requested_steps = int(first["summary"]["requested_steps"])
    for index, label in enumerate(labels):
        status = comparison["runs"][label]["controllers"]["deterministic_mpc"][
            "run_status"
        ]
        if status["completed_steps"] < requested_steps:
            panel.text(
                index,
                97,
                f"MPC fail\n@step {status['qp_failure_step']}",
                ha="center",
                va="top",
                fontsize=8,
                color="#8c1d40",
            )
    panel.set_ylim(0, 100)

    figure.suptitle(
        "copyBH: fixed-seed input-cap ablation (no controller/theory change)",
        fontsize=14,
    )
    png_path = output_dir / "copyBH_u_cap_sweep_comparison.png"
    figure.savefig(png_path, dpi=180)
    plt.close(figure)
    print(json.dumps({"json": str(json_path), "png": str(png_path)}, indent=2))
    return comparison


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run", action="append", nargs=2, metavar=("LABEL", "DIRECTORY"),
        required=True,
    )
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    analyze(args.run, args.output_dir)


if __name__ == "__main__":
    _main()
