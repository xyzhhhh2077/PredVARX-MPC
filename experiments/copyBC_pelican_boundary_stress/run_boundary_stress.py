"""Boundary-stress comparison of deterministic MPC and fixed-anchor SMPC.

Both controllers use the same identified Pelican model, cost, input bounds,
initial state, reference, and disturbance realization. They differ only in
whether predicted task constraints are tightened by propagated covariance.
"""

import argparse
import importlib.util
import json
from pathlib import Path

import cvxpy as cp
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
BB_SCRIPT = HERE.parent / "copyBB_pelican_fixed_anchor_smpc" / "run_fixed_anchor_closed_loop.py"


def _load_copybb():
    spec = importlib.util.spec_from_file_location("copybb_boundary_source", BB_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BB = _load_copybb()


def critical_reference(hard_bound, margin=0.02, axis=0):
    reference = np.zeros(2)
    reference[int(axis)] = float(hard_bound) - float(margin)
    return reference


def boundary_scenarios(hard_bound, margin=0.005):
    return {
        "x_critical": critical_reference(hard_bound, margin, axis=0),
        "y_critical": critical_reference(hard_bound, margin, axis=1),
    }


def moving_boundary_reference(hard_bound, steps, margin=0.005,
                              y_amplitude=1.0, period_steps=600):
    k = np.arange(int(steps))
    reference = np.zeros((int(steps), 2))
    reference[:, 0] = float(hard_bound) - float(margin)
    reference[:, 1] = float(y_amplitude) * np.sin(2.0 * np.pi * k / int(period_steps))
    return reference


def sample_disturbances(covariance, steps, seed):
    rng = np.random.default_rng(seed)
    draws = rng.multivariate_normal(np.zeros(covariance.shape[0]), covariance, size=steps)
    return draws[:, :, None]


def boundary_metrics(trajectory, hard_bound):
    trajectory = np.asarray(trajectory)
    margins = float(hard_bound) - np.abs(trajectory)
    return {
        "violation_steps": int(np.sum(np.any(margins < 0.0, axis=1))),
        "minimum_hard_margin": float(np.min(margins)),
        "maximum_abs_position": np.max(np.abs(trajectory), axis=0).tolist(),
    }


def controller_tightening(alpha_joint, axes, horizon, standard_deviation):
    if alpha_joint is None:
        return 0.0
    risk_each = float(alpha_joint) / (2.0 * int(axes) * int(horizon))
    return float(norm.ppf(1.0 - risk_each) * float(standard_deviation))


class BoundaryController(BB.FixedAnchorSMPC):
    def __init__(self, *args, alpha_joint=0.05, **kwargs):
        super().__init__(*args, alpha_joint=0.05 if alpha_joint is None else alpha_joint, **kwargs)
        self.use_tightening = alpha_joint is not None
        self.alpha_joint = alpha_joint

    def step(self, y_std, s_reference):
        z = self.model["R"].T @ (np.asarray(y_std).reshape(-1, 1) - self.y_mean)
        U0 = np.tile(self.u_mean, (self.N, 1))
        U = cp.Variable((self.N * self.nu, 1))
        r_full = self.E @ np.asarray(s_reference).reshape(-1, 1)
        cost = cp.quad_form(U - U0, np.kron(np.eye(self.N), self.Ru))
        constraints = [np.tile(self.u_min, self.N).reshape(-1, 1) <= U,
                       U <= np.tile(self.u_max, self.N).reshape(-1, 1)]
        chance_constraints = []
        risk_each = None if self.alpha_joint is None else (
            self.alpha_joint / (2.0 * self.E.shape[1] * self.N)
        )
        z_quantile = 0.0 if risk_each is None else norm.ppf(1.0 - risk_each)
        predictions = []
        tightened_bounds = []
        for j in range(self.N):
            mu = self.M[j] @ z + self.y_mean + self.G[j] @ (U - U0)
            predictions.append(mu)
            cost += cp.quad_form(mu - r_full, self.Q)
            Sigma_y = self.model["P"] @ self.Sigma_z[j] @ self.model["P"].T
            stage_bounds = []
            for axis in range(self.E.shape[1]):
                hq = self.E[:, axis].reshape(1, -1)
                variance = max(float((hq @ Sigma_y @ hq.T)[0, 0]), 1e-12)
                tight = z_quantile * np.sqrt(variance) if self.use_tightening else 0.0
                bound = self.task_bound - tight
                upper = hq @ mu <= bound
                lower = hq @ mu >= -bound
                constraints += [upper, lower]
                chance_constraints += [upper, lower]
                stage_bounds.append(bound)
            tightened_bounds.append(stage_bounds)
        problem = cp.Problem(cp.Minimize(cost), constraints)
        problem.solve(solver=cp.OSQP, verbose=False, warm_start=True)
        if problem.status not in ("optimal", "optimal_inaccurate"):
            return None, {"status": problem.status}

        stage_slacks = []
        for mu, bounds in zip(predictions, tightened_bounds):
            task_mu = (self.E.T @ mu.value).ravel()
            stage_slacks.extend(np.asarray(bounds) - np.abs(task_mu))
        dual_values = [float(np.max(np.asarray(c.dual_value))) for c in chance_constraints]
        info = {
            "status": problem.status,
            "minimum_predicted_slack": float(np.min(stage_slacks)),
            "maximum_constraint_dual": float(np.max(dual_values)),
        }
        control = np.clip(U.value[:self.nu].ravel(), self.u_min, self.u_max)
        return control, info


def initial_latent_state(model, E, initial_task):
    task_map = E.T @ model["P"]
    task_bias = (E.T @ model["y_mean"]).ravel()
    return np.linalg.pinv(task_map) @ (np.asarray(initial_task) - task_bias)[:, None]


def simulate_one(model, E, hard_bound, disturbances, reference, warmup_steps,
                 d=10, horizon=18, q_weight=80.0, ru=0.18, alpha_joint=0.05):
    steps = disturbances.shape[0]
    controller = BoundaryController(
        model, E, d=d, horizon=horizon, q_weight=q_weight, ru=ru,
        task_bound=hard_bound, alpha_joint=alpha_joint,
    )
    reference = np.asarray(reference)
    reference_trace = (
        np.tile(reference, (steps, 1)) if reference.ndim == 1 else reference.copy()
    )
    if reference_trace.shape != (steps, 2):
        raise ValueError("reference must have shape (2,) or (steps, 2)")
    z = initial_latent_state(model, E, reference_trace[0])
    u = model["u_mean"].ravel()
    trajectory = np.zeros((steps, 2))
    controls = np.zeros((steps, model["B"].shape[1]))
    min_slacks = []
    max_duals = []
    fallback = 0
    for k in range(steps):
        y = model["y_mean"] + model["P"] @ z
        trajectory[k] = (E.T @ y).ravel()
        target = np.zeros(2) if k < warmup_steps else reference_trace[k]
        if k % d == 0:
            candidate, info = controller.step(y.ravel(), target)
            if candidate is None:
                fallback += 1
                u = model["u_mean"].ravel()
            else:
                u = candidate
                min_slacks.append(info["minimum_predicted_slack"])
                max_duals.append(info["maximum_constraint_dual"])
        controls[k] = u
        z = BB.plant_step_with_noise(model, z, u, disturbances[k])
    metrics = boundary_metrics(trajectory, hard_bound)
    metrics.update({
        "fallback_count": fallback,
        "active_qp_steps": int(np.sum(np.asarray(min_slacks) <= 1e-4)),
        "positive_dual_qp_steps": int(np.sum(np.asarray(max_duals) > 1e-5)),
        "minimum_predicted_slack": float(np.min(min_slacks)),
        "maximum_constraint_dual": float(np.max(max_duals)),
        "rmse_to_reference_after_step": np.sqrt(
            np.mean((trajectory[warmup_steps:] - reference_trace[warmup_steps:]) ** 2, axis=0)
        ).tolist(),
    })
    traces = {
        "position": trajectory,
        "reference": reference_trace,
        "control": controls,
        "min_slack": np.asarray(min_slacks),
        "max_dual": np.asarray(max_duals),
    }
    return metrics, traces, controller


def plot_comparison(representatives, hard_bound, tightened_bound, path, seed):
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    colors = {"smpc": "#087e8b", "mpc": "#c44536"}
    for axis, (scenario, label) in enumerate((("x_critical", "x"), ("y_critical", "y"))):
        smpc, mpc = representatives[scenario]
        time = np.arange(len(smpc["position"])) / 100.0
        ax = axes[axis]
        ax.plot(time, smpc["reference"][:, axis], color="#222222", lw=1.2,
                label="critical reference")
        ax.plot(time, smpc["position"][:, axis], color=colors["smpc"], lw=1.2,
                label="SMPC, covariance tightening")
        ax.plot(time, mpc["position"][:, axis], color=colors["mpc"], lw=1.0,
                label="deterministic MPC, same noise")
        ax.axhline(hard_bound, color="#111111", ls="-.", lw=1.0, label="hard bound")
        ax.axhline(-hard_bound, color="#111111", ls="-.", lw=1.0)
        ax.axhline(tightened_bound[axis], color="#7b2cbf", ls="--", lw=1.0,
                   label="terminal tightened mean bound")
        ax.axhline(-tightened_bound[axis], color="#7b2cbf", ls="--", lw=1.0)
        ax.set_ylabel(f"standardized {label} position\n({label}-critical scenario)")
        ax.grid(True, alpha=0.25)
        ax.set_ylim(tightened_bound[axis] - 0.15, hard_bound + 0.08)
    axes[0].legend(loc="lower right", ncol=2, fontsize=8)
    axes[-1].set_xlabel("simulation time (s)")
    fig.suptitle(f"copyBC boundary stress, common disturbances, seed {seed}")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_moving_comparison(smpc, mpc, hard_bound, tightened_bound, path, seed):
    time = np.arange(len(smpc["position"])) / 100.0
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for axis, label in enumerate(("x", "y")):
        ax = axes[axis]
        ax.plot(time, smpc["reference"][:, axis], color="#222222", lw=1.2,
                label="moving reference")
        ax.plot(time, smpc["position"][:, axis], color="#087e8b", lw=1.1,
                label="SMPC")
        ax.plot(time, mpc["position"][:, axis], color="#c44536", lw=0.9,
                label="deterministic MPC")
        ax.axhline(hard_bound, color="#111111", ls="-.", lw=1.0,
                   label="hard bound")
        ax.axhline(tightened_bound[axis], color="#7b2cbf", ls="--", lw=1.0,
                   label="tightened mean bound")
        ax.set_xlabel("simulation time (s)")
        ax.set_ylabel(f"standardized {label} position")
        ax.grid(True, alpha=0.25)
        if axis == 0:
            ax.set_ylim(tightened_bound[axis] - 0.15, hard_bound + 0.08)
    axes[0].legend(loc="lower right", fontsize=8)

    ax = axes[2]
    ax.plot(smpc["reference"][:, 0], smpc["reference"][:, 1], color="#222222",
            lw=2.0, label="moving reference")
    ax.plot(smpc["position"][:, 0], smpc["position"][:, 1], color="#087e8b",
            lw=1.1, label="SMPC path")
    ax.plot(mpc["position"][:, 0], mpc["position"][:, 1], color="#c44536",
            lw=0.9, label="deterministic MPC path")
    ax.axvline(hard_bound, color="#111111", ls="-.", lw=1.0, label="hard x bound")
    ax.axvline(tightened_bound[0], color="#7b2cbf", ls="--", lw=1.0,
               label="tightened x mean bound")
    ax.set_xlim(tightened_bound[0] - 0.15, hard_bound + 0.08)
    ax.set_xlabel("standardized x position")
    ax.set_ylabel("standardized y position")
    ax.set_title("x-y motion along the boundary")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.suptitle(f"copyBC moving boundary stress, common disturbances, seed {seed}")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run_moving_experiment(seeds, steps=1200, margin=0.005, y_amplitude=1.0,
                          period_steps=600, d=10, horizon=18,
                          q_weight=80.0, ru=0.18):
    model, E, _, hard_bound, _ = BB.load_model_and_data()
    reference = moving_boundary_reference(
        hard_bound, steps, margin, y_amplitude, period_steps,
    )
    metrics = {"smpc": {}, "deterministic_mpc": {}}
    representative = None
    controller = None
    for seed in seeds:
        disturbances = sample_disturbances(model["Sigma_eps"], steps, seed)
        smpc_metrics, smpc_traces, controller = simulate_one(
            model, E, hard_bound, disturbances, reference, 0,
            d, horizon, q_weight, ru, alpha_joint=0.05,
        )
        mpc_metrics, mpc_traces, _ = simulate_one(
            model, E, hard_bound, disturbances, reference, 0,
            d, horizon, q_weight, ru, alpha_joint=None,
        )
        metrics["smpc"][str(seed)] = smpc_metrics
        metrics["deterministic_mpc"][str(seed)] = mpc_metrics
        if representative is None:
            representative = (smpc_traces, mpc_traces)
        print(
            f"moving seed {seed}: violations SMPC={smpc_metrics['violation_steps']} "
            f"MPC={mpc_metrics['violation_steps']}"
        )

    risk_each = 0.05 / (2.0 * E.shape[1] * horizon)
    quantile = norm.ppf(1.0 - risk_each)
    Sigma_y_terminal = model["P"] @ controller.Sigma_z[-1] @ model["P"].T
    tightening = []
    for axis in range(E.shape[1]):
        row = E[:, axis].reshape(1, -1)
        variance = max(float((row @ Sigma_y_terminal @ row.T)[0, 0]), 1e-12)
        tightening.append(quantile * np.sqrt(variance))
    tightened_bound = hard_bound - np.asarray(tightening)
    summary = {
        "scope": "moving-reference identified-model-in-the-loop boundary stress test",
        "path": "x fixed near positive hard bound; y sinusoidal",
        "hard_bound": hard_bound,
        "reference_margin_to_hard_bound": margin,
        "y_amplitude": y_amplitude,
        "period_steps": period_steps,
        "terminal_tightening": tightening,
        "steps": steps,
        "seeds": list(seeds),
        "metrics": metrics,
    }
    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    (results / "copyBC_moving_boundary.json").write_text(
        json.dumps(summary, indent=2), encoding="ascii"
    )
    np.savez(results / "copyBC_moving_boundary.npz", summary=summary,
             smpc=representative[0], deterministic_mpc=representative[1])
    plot_moving_comparison(
        representative[0], representative[1], hard_bound, tightened_bound,
        results / "copyBC_moving_boundary.png", seeds[0],
    )
    return summary


def run_experiment(seeds, steps=1200, warmup_steps=0, margin=0.005,
                   d=10, horizon=18, q_weight=80.0, ru=0.18):
    model, E, _, hard_bound, _ = BB.load_model_and_data()
    scenarios = boundary_scenarios(hard_bound, margin)
    all_metrics = {}
    representatives = {}
    controller = None
    for scenario, reference in scenarios.items():
        all_metrics[scenario] = {"smpc": {}, "deterministic_mpc": {}}
        for seed in seeds:
            disturbances = sample_disturbances(model["Sigma_eps"], steps, seed)
            smpc_metrics, smpc_traces, controller = simulate_one(
                model, E, hard_bound, disturbances, reference, warmup_steps,
                d, horizon, q_weight, ru, alpha_joint=0.05,
            )
            mpc_metrics, mpc_traces, _ = simulate_one(
                model, E, hard_bound, disturbances, reference, warmup_steps,
                d, horizon, q_weight, ru, alpha_joint=None,
            )
            all_metrics[scenario]["smpc"][str(seed)] = smpc_metrics
            all_metrics[scenario]["deterministic_mpc"][str(seed)] = mpc_metrics
            if scenario not in representatives:
                representatives[scenario] = (smpc_traces, mpc_traces)
            print(
                f"{scenario} seed {seed}: violations SMPC={smpc_metrics['violation_steps']} "
                f"MPC={mpc_metrics['violation_steps']}; active QPs "
                f"SMPC={smpc_metrics['active_qp_steps']} MPC={mpc_metrics['active_qp_steps']}"
            )

    risk_each = 0.05 / (2.0 * E.shape[1] * horizon)
    quantile = norm.ppf(1.0 - risk_each)
    Sigma_y_terminal = model["P"] @ controller.Sigma_z[-1] @ model["P"].T
    tightening = []
    for axis in range(E.shape[1]):
        row = E[:, axis].reshape(1, -1)
        variance = max(float((row @ Sigma_y_terminal @ row.T)[0, 0]), 1e-12)
        tightening.append(quantile * np.sqrt(variance))
    tightened_bound = hard_bound - np.asarray(tightening)
    summary = {
        "scope": "identified-model-in-the-loop boundary stress test",
        "comparison": "SMPC versus deterministic MPC with common random disturbances",
        "hard_bound": hard_bound,
        "terminal_tightening": tightening,
        "terminal_tightened_mean_bound": tightened_bound.tolist(),
        "critical_references": {key: value.tolist() for key, value in scenarios.items()},
        "reference_margin_to_hard_bound": margin,
        "steps": steps,
        "warmup_steps": warmup_steps,
        "seeds": list(seeds),
        "d": d,
        "horizon": horizon,
        "q_weight": q_weight,
        "ru": ru,
        "metrics": all_metrics,
    }
    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    (results / "copyBC_boundary_stress.json").write_text(
        json.dumps(summary, indent=2), encoding="ascii"
    )
    np.savez(results / "copyBC_boundary_stress.npz", summary=summary,
             representatives=representatives)
    plot_comparison(
        representatives, hard_bound, tightened_bound,
        results / "copyBC_boundary_stress.png", seeds[0],
    )
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", nargs="+", type=int, default=[7, 19, 31])
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--warmup-steps", type=int, default=0)
    parser.add_argument("--margin", type=float, default=0.005)
    args = parser.parse_args()
    run_experiment(args.seeds, args.steps, args.warmup_steps, args.margin)
    run_moving_experiment(args.seeds, args.steps, args.margin)


if __name__ == "__main__":
    main()
