"""Fair paired 20-seed Frozen/copyBJ/copyBL benchmark.

This is a same-disturbance paired benchmark, not a Monte Carlo probability
calibration. The three existing controllers are imported and called unchanged.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BH_SCRIPT = ROOT / "experiments/copyBH_pelican_u_cap_ablation/run_u_cap_ablation.py"
BJ_SCRIPT = ROOT / "experiments/copyBJ_pelican_online_covariance_smpc/run_online_covariance_smpc.py"
BL_SCRIPT = ROOT / "experiments/copyBL_pelican_dual_gated_covariance_smpc/run_dual_gated_covariance_smpc.py"

def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod

BH = _load("copybm_bh", BH_SCRIPT)
BJ = _load("copybm_bj", BJ_SCRIPT)
BL = _load("copybm_bl", BL_SCRIPT)
BI = BJ.BI

CONTROLLERS = ("Frozen", "BJ", "BL")
ONLINE_WEIGHT = 0.8
HARD_BOUND = float(BI.EXPERIMENT_HARD_BOUND)
DEFAULT_STEPS = int(BI.DEFAULT_STEPS)
CONFIG_VERSION = "copyBM-fair-v1"


def seed_pairs():
    return [(16 + i, 7016 + i) for i in range(20)]


def _sha(*arrays):
    h = hashlib.sha256()
    for array in arrays:
        a = np.ascontiguousarray(np.asarray(array, dtype=np.float64))
        h.update(str(a.shape).encode()); h.update(a.tobytes())
    return h.hexdigest()


def disturbance_sha256(innovation, wind_latent_increment):
    return BH._disturbance_sha256(np.asarray(innovation), np.asarray(wind_latent_increment))


def make_pair_inputs(steps, innovation_seed, wind_seed, *, model=None, E=None, scales=None):
    """Create the pair's arrays once; returned arrays are never mutated by drivers."""
    steps = int(steps)
    if model is None:
        innovation = np.random.default_rng(int(innovation_seed)).normal(size=(steps, 4))
        wind_latent = np.random.default_rng(int(wind_seed)).normal(size=(steps, 4))
        wind_velocity = np.random.default_rng(int(wind_seed)).normal(size=(steps, 3))
    else:
        innovation = BH.innovation_sequence(model, steps, int(innovation_seed))
        wind_velocity = BH.physical_wind_velocity(steps, BH.SIGMA_WIND_MPS, int(wind_seed))
        task_map = E.T @ model["P"][:, :model["R"].shape[1]]
        wind_latent, wind_map = BH.wind_velocity_to_latent_increment(
            wind_velocity, task_map, scales["y_scale"][:3])
        if not np.allclose(wind_map, BH.build_process_covariances(
                model, E, scales["y_scale"][:3], BH.SIGMA_WIND_MPS)["G_w"], atol=1e-14):
            raise RuntimeError("wind map mismatch")
    reference = (BI.hard_edge_face_petal_reference(steps) if steps >= 120
                 else np.zeros((steps, 3), dtype=float))
    return {"reference": reference,
            "innovation": np.asarray(innovation, dtype=float).copy(),
            "wind_velocity": np.asarray(wind_velocity, dtype=float).copy(),
            "wind_latent_increment": np.asarray(wind_latent, dtype=float).copy(),
            "initial_task": np.zeros(3, dtype=float)}


def strict_json_dumps(value):
    def safe(x):
        if isinstance(x, dict): return {str(k): safe(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)): return [safe(v) for v in x]
        if isinstance(x, np.ndarray): return safe(x.tolist())
        if isinstance(x, (np.floating, float)): return float(x) if np.isfinite(x) else None
        if isinstance(x, np.integer): return int(x)
        if isinstance(x, np.bool_): return bool(x)
        return x
    return json.dumps(safe(value), indent=2, allow_nan=False, sort_keys=True)


def _none_if_incomplete(value, complete):
    return value if complete else None


def controller_metrics(result, reference, control_mean, requested_steps):
    completed = int(result.get("completed_steps", 0))
    complete = completed == int(requested_steps)
    S = np.asarray(result.get("S", []), dtype=float).reshape(-1, 3)
    U = np.asarray(result.get("U", []), dtype=float)
    failures = int(result.get("qp_failure_count", 0))
    qp_count = int(result.get("qp_count", 0))
    if "qp_count" not in result:
        # BL exposes every control solve through its multirate loop but does not
        # duplicate BJ/BH's qp_count field; reconstruct the same denominator.
        qp_count = int(np.ceil(completed / BJ.BH.CONTROL_INTERVAL_STEPS)) + failures
    metric = {"completed_steps": completed,
              "complete": complete,
              "rmse_xyz": (np.sqrt(np.mean((S - reference[:completed]) ** 2, axis=0)).tolist()
                           if complete and completed else None),
              "cumulative_realized_stage_cost": None,
              "minimum_hard_margin": (float(result.get("minimum_hard_margin", float("nan"))) if complete else None),
              "qp_count": qp_count,
              "qp_failure_count": failures,
              "qp_failure_rate": (failures / qp_count if qp_count else None)}
    if complete:
        cost = BI.realized_stage_cost(S, reference[:completed], U, control_mean)
        metric["cumulative_realized_stage_cost"] = float(np.sum(cost))
    return metric


def gate_metrics(events):
    events = list(events)
    ftr = sum(bool(e.get("feasibility_triggered", False)) for e in events)
    frej = sum(not bool(e.get("feasibility_accepted", False)) for e in events if e.get("feasibility_triggered", False))
    vtr = sum(bool(e.get("value_triggered", False)) for e in events)
    vacc = sum(bool(e.get("value_accepted", False)) for e in events if e.get("value_triggered", False))
    return {"feasibility_trigger_count": ftr, "feasibility_rejection_count": frej,
            "value_trigger_count": vtr, "value_accept_count": vacc,
            "value_acceptance_rate": vacc / vtr if vtr else None,
            "value_rejection_count": vtr - vacc,
            "update_committed_count": sum(bool(e.get("update_committed", False)) for e in events)}


def _config(steps):
    return {"version": CONFIG_VERSION, "steps": int(steps), "seed_pairs": seed_pairs(),
            "hard_bound": HARD_BOUND, "online_weight": ONLINE_WEIGHT,
            "initial_task": [0.0, 0.0, 0.0], "controllers": list(CONTROLLERS)}


def config_sha256(steps):
    return hashlib.sha256(strict_json_dumps(_config(steps)).encode()).hexdigest()


def load_resume(path, expected_config_sha):
    if not Path(path).exists(): return {"config_sha256": expected_config_sha, "seed_pairs": seed_pairs(), "seeds": []}
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("config_sha256") != expected_config_sha or data.get("seed_pairs") != seed_pairs():
        raise ValueError("resume configuration or seed pair configuration mismatch")
    return data


def _atomic_write(path, value):
    path = Path(path); tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(strict_json_dumps(value), encoding="utf-8"); os.replace(tmp, path)


def _run_one(model, E, covariance, noise, inputs, steps, innovation_seed, wind_seed, initial_sha):
    # Inputs are explicit copies: no controller can mutate the pair's canonical arrays.
    args = (model, E, HARD_BOUND, inputs["reference"].copy(), inputs["innovation"].copy(),
            inputs["wind_latent_increment"].copy(), covariance["Sigma_eps_aug"])
    frozen = BH.simulate_controller(*args, BH.ALPHA_JOINT, initial_task=np.zeros(3))
    bj = BJ.simulate_online_controller(*args, noise["Sigma_obs_proxy"], initial_task=np.zeros(3), online_weight=ONLINE_WEIGHT)
    bl = BL.simulate_dual_gated_controller(*args, initial_task=np.zeros(3))
    sha = disturbance_sha256(inputs["innovation"], inputs["wind_latent_increment"])
    assert frozen["disturbance_sha256"] == bj["disturbance_sha256"] == bl["disturbance_sha256"] == sha
    out = {"seed": {"innovation_seed": innovation_seed, "wind_seed": wind_seed,
                     "disturbance_sha256": sha, "initial_state_sha256": initial_sha}, "controllers": {}}
    for name, result in (("Frozen", frozen), ("BJ", bj), ("BL", bl)):
        m = controller_metrics(result, inputs["reference"], model["u_mean"].ravel(), steps)
        if name == "Frozen": m.update({"gate_acceptance_rate": None})
        elif name == "BJ":
            attempts = int(np.ceil(steps / BJ.BH.CONTROL_INTERVAL_STEPS))
            m.update({"gate_acceptance_rate": None,
                      "online_covariance_update_count": attempts if m["complete"] else None,
                      "online_covariance_update_attempt_count": attempts})
        else:
            gm = gate_metrics(result.get("events", [])); m.update(gm); m["gate_acceptance_rate"] = gm["value_acceptance_rate"]
        m["disturbance_sha256"] = result["disturbance_sha256"]
        out["controllers"][name] = m
    return out


def aggregate(seeds):
    out = {"seed_count": len(seeds), "controllers": {}}
    for name in CONTROLLERS:
        rows = [s["controllers"][name] for s in seeds]
        vals = {}
        for key in ("rmse_xyz", "cumulative_realized_stage_cost", "minimum_hard_margin", "qp_failure_rate", "gate_acceptance_rate"):
            series = [r[key] for r in rows if r.get(key) is not None]
            if key == "rmse_xyz": vals[key] = {axis: _summary([x[i] for x in series]) for i, axis in enumerate("xyz")}
            else: vals[key] = _summary(series)
        out["controllers"][name] = vals
    return out


def _summary(values):
    if not values: return {"n": 0, "mean": None, "std": None, "ci95": None}
    a = np.asarray(values, dtype=float); mean = float(np.mean(a)); std = float(np.std(a, ddof=1)) if len(a) > 1 else 0.0
    return {"n": len(a), "mean": mean, "std": std, "ci95": 1.96 * std / np.sqrt(len(a))}


def save_figure(path, summary):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    metrics = [("cumulative_realized_stage_cost", "Cumulative realized stage cost"), ("minimum_hard_margin", "Minimum hard margin"), ("qp_failure_rate", "QP failure rate")]
    colors = {"Frozen": "#4c78a8", "BJ": "#f58518", "BL": "#54a24b"}
    for ax, (key, title) in zip(axes, metrics):
        for name in CONTROLLERS:
            s = summary["controllers"][name][key]
            if s["mean"] is not None: ax.errorbar(name, s["mean"], yerr=s["ci95"], fmt="o", color=colors[name], capsize=4)
        ax.set_title(title); ax.grid(alpha=.25, axis="y")
    fig.suptitle(f"copyBM same-disturbance paired benchmark (n={summary['seed_count']})\nmean ± 95% CI")
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close(fig)


def run(output_dir=None, steps=DEFAULT_STEPS, max_seeds=20):
    steps = int(steps); max_seeds = min(int(max_seeds), 20)
    results = Path(output_dir) if output_dir else HERE / "results"; results.mkdir(parents=True, exist_ok=True)
    cfg_sha = config_sha256(steps); state = load_resume(results / "per_seed.json", cfg_sha)
    done = {int(s["seed"]["innovation_seed"]): s for s in state.get("seeds", [])}
    model, E, _source, scales, noise = BH.load_identification_and_noise_objects()
    covariance = BH.build_process_covariances(model, E, scales["y_scale"][:3], BH.SIGMA_WIND_MPS)
    all_seeds = []
    for innovation_seed, wind_seed in seed_pairs()[:max_seeds]:
        if innovation_seed in done: all_seeds.append(done[innovation_seed]); continue
        inputs = make_pair_inputs(steps, innovation_seed, wind_seed, model=model, E=E, scales=scales)
        state0, _u0 = BH.equilibrium_initial_condition(model, E, inputs["initial_task"])
        initial_sha = _sha(state0)
        try:
            row = _run_one(model, E, covariance, noise, inputs, steps, innovation_seed, wind_seed, initial_sha)
            row["status"] = "success"
        except Exception as exc:
            # A failed pair is evidence too: retain it rather than dropping it or
            # converting a partial trajectory into a successful comparison.
            row = {"status": "failed", "error": f"{type(exc).__name__}: {exc}",
                   "seed": {"innovation_seed": innovation_seed, "wind_seed": wind_seed,
                            "disturbance_sha256": disturbance_sha256(inputs["innovation"], inputs["wind_latent_increment"]),
                            "initial_state_sha256": initial_sha},
                   "controllers": {name: {"complete": False, "completed_steps": 0,
                                          "rmse_xyz": None, "cumulative_realized_stage_cost": None,
                                          "minimum_hard_margin": None, "qp_count": None,
                                          "qp_failure_count": None, "qp_failure_rate": None,
                                          "gate_acceptance_rate": None}
                                   for name in CONTROLLERS}}
        done[innovation_seed] = row; all_seeds = list(done.values())
        all_seeds.sort(key=lambda r: r["seed"]["innovation_seed"])
        _atomic_write(results / "per_seed.json", {"config_sha256": cfg_sha, "seed_pairs": seed_pairs(), "seeds": all_seeds})
    summary = aggregate(all_seeds)
    summary.update({"config_sha256": cfg_sha, "seed_pairs": seed_pairs(), "steps": steps,
                    "scope": "same-disturbance paired benchmark; exploratory evidence, not Monte Carlo probability calibration"})
    _atomic_write(results / "aggregate.json", summary); save_figure(results / "aggregate_comparison.png", summary)
    with (results / "per_seed.csv").open("w", newline="", encoding="utf-8") as f:
        fields = ["innovation_seed", "wind_seed", "controller", "completed_steps", "rmse_xyz", "cumulative_realized_stage_cost", "minimum_hard_margin", "qp_count", "qp_failure_count", "qp_failure_rate", "gate_acceptance_rate"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for row in all_seeds:
            for name, m in row["controllers"].items():
                w.writerow({"innovation_seed": row["seed"]["innovation_seed"], "wind_seed": row["seed"]["wind_seed"], "controller": name, **{k: json.dumps(m.get(k)) if isinstance(m.get(k), list) else m.get(k) for k in fields[3:]}})
    return summary


def main():
    p = argparse.ArgumentParser(); p.add_argument("--output-dir", type=Path); p.add_argument("--steps", type=int, default=DEFAULT_STEPS); p.add_argument("--max-seeds", type=int, default=20)
    args = p.parse_args()
    print(strict_json_dumps(run(args.output_dir, args.steps, args.max_seeds)))

if __name__ == "__main__": main()
