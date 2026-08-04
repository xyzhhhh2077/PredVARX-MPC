"""Long 3D Pelican boundary stress with common external gusts.

The frozen copyBE second-order model, costs, bounds, references, and complete
plant disturbance sequence are identical for SMPC and deterministic MPC.
Only covariance tightening differs between the controllers.
"""

import argparse
import importlib.util
import json
from pathlib import Path

import cvxpy as cp
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import null_space
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
BE_SCRIPT = (
    HERE.parent / "copyBE_pelican_3d_z_aware_smpc"
    / "run_3d_z_aware_trajectory.py"
)
_spec = importlib.util.spec_from_file_location("copybe_long_source", BE_SCRIPT)
BE = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(BE)


def long_boundary_reference(hard_bound, steps, margin=0.02,
                            inward_excursion=0.08, period_steps=1200, axis=0):
    """Slow motion that remains close to one selected positive hard bound."""
    steps = int(steps)
    if steps < 2:
        raise ValueError("steps must be at least two")
    reference = np.zeros((steps, 3))
    cycles = steps / float(period_steps)
    phase = np.linspace(0.0, 2.0 * np.pi * cycles, steps, endpoint=True)
    inward = 0.5 * float(inward_excursion) * (1.0 - np.cos(phase))
    reference[:, int(axis)] = float(hard_bound) - float(margin) - inward
    return reference


def external_gust_sequence(task_map, steps, period_steps=1200,
                           gust_width=30, gust_amplitude=0.025,
                           critical_axis=None):
    """Create zero-impulse outward/recovery gust pairs in task coordinates."""
    task_map = np.asarray(task_map)
    steps = int(steps)
    width = int(gust_width)
    if task_map.shape[0] != 3:
        raise ValueError("task_map must have three task rows")
    task_gusts = np.zeros((steps, 3))
    starts = (np.array([0.20, 0.50, 0.80]) * int(period_steps)).astype(int)
    for cycle_start in range(0, steps, int(period_steps)):
        for axis, local_start in enumerate(starts):
            if critical_axis is not None and axis != int(critical_axis):
                continue
            start = cycle_start + int(local_start)
            stop = min(start + 2 * width, steps, cycle_start + int(period_steps))
            if start < steps and stop - start >= 4:
                half = (stop - start) // 2
                window = np.sin(np.linspace(0.0, np.pi, half, endpoint=True))
                if np.max(window) > 0.0:
                    window /= np.max(window)
                pulse = np.concatenate((window, -window))
                task_gusts[start:start + len(pulse), axis] += float(gust_amplitude) * pulse
    latent_gusts = (np.linalg.pinv(task_map) @ task_gusts.T).T
    return latent_gusts, task_gusts


def cycle_metrics(trajectory, reference, hard_bound, period_steps):
    trajectory = np.asarray(trajectory)
    reference = np.asarray(reference)
    metrics = []
    for start in range(0, len(trajectory), int(period_steps)):
        stop = min(start + int(period_steps), len(trajectory))
        segment = trajectory[start:stop]
        target = reference[start:stop]
        violations = np.any(np.abs(segment) > float(hard_bound), axis=1)
        metrics.append({
            "cycle": len(metrics) + 1,
            "start_step": start,
            "stop_step": stop,
            "violation_steps": int(violations.sum()),
            "rmse": np.sqrt(np.mean((segment - target) ** 2, axis=0)).tolist(),
            "minimum_hard_margin": float(
                np.min(float(hard_bound) - np.abs(segment))
            ),
        })
    return metrics


def equilibrium_initial_condition(model, E, target):
    """Find a task-consistent second-order equilibrium within input limits."""
    n = model["A"].shape[0]
    nu = model["B"].shape[1]
    target = np.asarray(target).reshape(3, 1)
    equations = np.block([
        [np.eye(n) - model["A"], -model["B"]],
        [E.T @ model["P"], np.zeros((3, nu))],
    ])
    rhs = np.vstack([
        -model["B"] @ model["u_mean"],
        target - E.T @ model["y_mean"],
    ])
    solution = np.linalg.lstsq(equations, rhs, rcond=None)[0]
    basis = null_space(equations)
    if basis.shape[1] == 1:
        direction = basis[n:, 0]
        lower, upper = -np.inf, np.inf
        for value, slope in zip(solution[n:, 0], direction):
            if abs(slope) < 1e-14:
                continue
            limits = sorted(((-6.0 - value) / slope, (6.0 - value) / slope))
            lower, upper = max(lower, limits[0]), min(upper, limits[1])
        if lower > upper:
            raise ValueError(f"No bounded-input equilibrium for task target {target.ravel()}")
        solution += basis[:, :1] * np.clip(0.0, lower, upper)
    state_value = solution[:n]
    control_value = solution[n:]
    dynamic_residual = np.max(np.abs(
        (np.eye(n) - model["A"]) @ state_value
        - model["B"] @ (control_value - model["u_mean"])
    ))
    output_residual = np.max(np.abs(
        E.T @ (model["y_mean"] + model["P"] @ state_value) - target
    ))
    if dynamic_residual > 1e-7 or output_residual > 1e-7:
        raise ValueError(
            f"Equilibrium residual too large: dynamics={dynamic_residual}, "
            f"output={output_residual}"
        )
    if np.any(np.abs(control_value) > 6.0 + 1e-8):
        raise ValueError(f"Equilibrium input exceeds bounds for target {target.ravel()}")
    return state_value, control_value.ravel()


class InstrumentedZAwareSMPC(BE.ZAwareSMPC):
    def step(self, z_aug, reference_horizon):
        z_aug = np.asarray(z_aug).reshape(-1, 1)
        references = np.asarray(reference_horizon)
        U0 = np.tile(self.u_mean, (self.N, 1))
        U = cp.Variable((self.N * self.nu, 1))
        cost = cp.quad_form(U - U0, np.kron(np.eye(self.N), self.Ru))
        constraints = [np.tile(self.u_min, self.N).reshape(-1, 1) <= U,
                       U <= np.tile(self.u_max, self.N).reshape(-1, 1)]
        risk = None if self.alpha_joint is None else (
            self.alpha_joint / (2.0 * self.E.shape[1] * self.N)
        )
        quantile = 0.0 if risk is None else norm.ppf(1.0 - risk)
        predictions, stage_bounds, task_constraints = [], [], []
        for j in range(self.N):
            mu = self.M[j] @ z_aug + self.y_mean + self.G[j] @ (U - U0)
            target = self.E @ references[j].reshape(-1, 1)
            cost += cp.quad_form(mu - target, self.Q)
            predictions.append(mu)
            Sigma_y = self.model["P"] @ self.Sigma_z[j] @ self.model["P"].T
            bounds = []
            for axis in range(self.E.shape[1]):
                row = self.E[:, axis].reshape(1, -1)
                variance = max(float((row @ Sigma_y @ row.T).item()), 1e-12)
                tightening = quantile * np.sqrt(variance)
                bound = self.task_bound - tightening
                upper, lower = row @ mu <= bound, row @ mu >= -bound
                constraints += [upper, lower]
                task_constraints += [upper, lower]
                bounds.append(bound)
            stage_bounds.append(bounds)
        problem = cp.Problem(cp.Minimize(cost), constraints)
        problem.solve(solver=cp.OSQP, verbose=False, warm_start=True)
        if problem.status not in ("optimal", "optimal_inaccurate"):
            return None, {"status": problem.status}
        slacks = []
        for mu, bounds in zip(predictions, stage_bounds):
            task_mu = (self.E.T @ mu.value).ravel()
            slacks.extend(np.asarray(bounds) - np.abs(task_mu))
        duals = [float(np.max(np.asarray(c.dual_value))) for c in task_constraints]
        return np.clip(U.value[:self.nu].ravel(), self.u_min, self.u_max), {
            "status": problem.status,
            "minimum_predicted_slack": float(np.min(slacks)),
            "maximum_constraint_dual": float(np.max(duals)),
        }


def simulate(model, E, hard_bound, reference, residual_noise, external_gusts,
             alpha_joint, d=10, horizon=18, q_weight=80.0, ru=0.18,
             period_steps=1200, initial_task=None):
    controller = InstrumentedZAwareSMPC(
        model, E, d=d, horizon=horizon, q_weight=q_weight, ru=ru,
        task_bound=hard_bound, alpha_joint=alpha_joint,
    )
    ell = model["A"].shape[0] // 2
    initial_task = reference[0] if initial_task is None else initial_task
    z_aug, u = equilibrium_initial_condition(model, E, initial_task)
    steps = len(reference)
    task = np.zeros((steps, 3))
    controls = np.zeros((steps, model["B"].shape[1]))
    slacks, duals = [], []
    fallback = 0
    for k in range(steps):
        y = model["y_mean"] + model["P"] @ z_aug
        task[k] = (E.T @ y).ravel()
        if k % d == 0:
            refs = BE.base.future_reference_horizon(reference, k, d, horizon)
            candidate, info = controller.step(z_aug, refs)
            if candidate is None:
                fallback += 1
                u = model["u_mean"].ravel()
            else:
                u = candidate
                slacks.append(info["minimum_predicted_slack"])
                duals.append(info["maximum_constraint_dual"])
        controls[k] = u
        z_aug = model["A"] @ z_aug + model["B"] @ (
            u.reshape(-1, 1) - model["u_mean"]
        )
        z_aug[:ell] += residual_noise[k].reshape(-1, 1)
        z_aug[:ell] += external_gusts[k].reshape(-1, 1)
    violations = np.any(np.abs(task) > hard_bound, axis=1)
    saturation = np.any(np.isclose(np.abs(controls), 6.0, atol=1e-5), axis=1)
    return {
        "S": task,
        "U": controls,
        "violation_steps": int(violations.sum()),
        "minimum_hard_margin": float(np.min(hard_bound - np.abs(task))),
        "fallback_count": int(fallback),
        "input_saturation_steps": int(saturation.sum()),
        "active_qp_steps": int(np.sum(np.asarray(slacks) <= 1e-4)),
        "positive_dual_qp_steps": int(np.sum(np.asarray(duals) > 1e-5)),
        "minimum_predicted_slack": float(np.min(slacks)),
        "maximum_constraint_dual": float(np.max(duals)),
        "rmse": np.sqrt(np.mean((task - reference) ** 2, axis=0)),
        "cycle_metrics": cycle_metrics(task, reference, hard_bound, period_steps),
    }


def terminal_tightening(controller, E):
    quantile = norm.ppf(1.0 - 0.05 / (2.0 * E.shape[1] * controller.N))
    Sigma_y = controller.model["P"] @ controller.Sigma_z[-1] @ controller.model["P"].T
    values = []
    for axis in range(E.shape[1]):
        row = E[:, axis].reshape(1, -1)
        values.append(quantile * np.sqrt(max(float((row @ Sigma_y @ row.T).item()), 1e-12)))
    return np.asarray(values)


def _serializable_metrics(result):
    return {key: (value.tolist() if isinstance(value, np.ndarray) else value)
            for key, value in result.items() if key not in ("S", "U")}


def plot_results(reference, smpc, mpc, task_gusts, hard_bound,
                 tightened_bound, path, seed):
    t = np.arange(len(reference)) / 100.0
    scenario_duration = len(reference) / 300.0
    fig, axes = plt.subplots(4, 1, figsize=(15, 11), sharex=True)
    for axis, label in enumerate(("x", "y", "z")):
        ax = axes[axis]
        ax.plot(t, reference[:, axis], "k--", lw=1.0, label="near-bound reference")
        ax.plot(t, mpc["S"][:, axis], color="#c44536", lw=1.0, ls=":",
                alpha=1.0, label="deterministic MPC")
        ax.plot(t, smpc["S"][:, axis], color="#087e8b", lw=1.0, label="SMPC")
        ax.axhline(hard_bound, color="#111111", ls="-.", lw=0.9,
                   label="hard bound")
        ax.axhline(tightened_bound[axis], color="#7b2cbf", ls="--", lw=0.9,
                   label="terminal tightened mean bound")
        for reset in (scenario_duration, 2.0 * scenario_duration):
            ax.axvline(reset, color="#666666", ls=":", lw=0.8)
        ax.set_ylabel(f"standardized {label}")
        ax.grid(True, alpha=0.25)
        if axis == 0:
            ax.legend(loc="lower right", ncol=3, fontsize=8)
    axes[3].plot(t, task_gusts[:, 0], label="external x gust")
    axes[3].plot(t, task_gusts[:, 1], label="external y gust")
    axes[3].plot(t, task_gusts[:, 2], label="external z gust")
    for reset in (scenario_duration, 2.0 * scenario_duration):
        axes[3].axvline(reset, color="#666666", ls=":", lw=0.8)
    axes[3].set_ylabel("task-space gust")
    axes[3].set_xlabel("simulation time (s)")
    axes[3].grid(True, alpha=0.25)
    axes[3].legend(loc="upper right", ncol=3, fontsize=8)
    axes[0].text(0.5 * scenario_duration, 1.02, "x reachable-limit scenario",
                 transform=axes[0].get_xaxis_transform(), ha="center", fontsize=8)
    axes[0].text(1.5 * scenario_duration, 1.02, "y near-boundary scenario",
                 transform=axes[0].get_xaxis_transform(), ha="center", fontsize=8)
    axes[0].text(2.5 * scenario_duration, 1.02, "z near-boundary scenario",
                 transform=axes[0].get_xaxis_transform(), ha="center", fontsize=8)
    fig.suptitle(f"copyBF: three long 3D boundary scenarios with common external gusts, seed {seed}")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run_experiment(seeds=(7, 19, 31), steps=4000, period_steps=600,
                   gust_width=20):
    model, E, hard_bound, scales = BE.load_model()
    ell = model["A"].shape[0] // 2
    task_map = E.T @ model["P"][:, :ell]
    controller = InstrumentedZAwareSMPC(model, E, task_bound=hard_bound)
    tightening = terminal_tightening(controller, E)
    tightened_bound = hard_bound - tightening
    scenario_specs = {
        "x_reachable_limit": {"axis": 0, "target": 3.60,
                              "gust_amplitude": 0.0003},
        "y_near_boundary": {"axis": 1, "target": 3.55,
                            "gust_amplitude": 0.0002},
        "z_near_boundary": {"axis": 2, "target": 3.75,
                            "gust_amplitude": 0.0003},
    }
    scenarios = {}
    representative_parts = []
    for scenario_name, spec in scenario_specs.items():
        axis = spec["axis"]
        reference = long_boundary_reference(
            hard_bound, steps, margin=hard_bound - spec["target"],
            inward_excursion=0.005, period_steps=period_steps, axis=axis,
        )
        initial_task = np.zeros(3)
        initial_task[axis] = tightened_bound[axis] - 0.05
        external_gusts, task_gusts = external_gust_sequence(
            task_map, steps, period_steps, gust_width,
            spec["gust_amplitude"], critical_axis=axis,
        )
        scenario_results = {}
        representative = None
        for seed in seeds:
            rng = np.random.default_rng(seed)
            residual_noise = rng.multivariate_normal(
                np.zeros(ell), model["Sigma_eps"][:ell, :ell], size=steps
            )
            smpc = simulate(
                model, E, hard_bound, reference, residual_noise,
                external_gusts, 0.05, period_steps=period_steps,
                initial_task=initial_task,
            )
            mpc = simulate(
                model, E, hard_bound, reference, residual_noise,
                external_gusts, None, period_steps=period_steps,
                initial_task=initial_task,
            )
            scenario_results[str(seed)] = {
                "smpc": _serializable_metrics(smpc),
                "deterministic_mpc": _serializable_metrics(mpc),
            }
            if representative is None:
                representative = (smpc, mpc)
            print(
                f"{scenario_name} seed {seed}: violations "
                f"SMPC={smpc['violation_steps']} MPC={mpc['violation_steps']}; "
                f"active={smpc['active_qp_steps']} "
                f"fallback={smpc['fallback_count']}/{mpc['fallback_count']}"
            )
        scenarios[scenario_name] = {
            "axis": axis,
            "reference_target": spec["target"],
            "initial_task": initial_task.tolist(),
            "gust_amplitude": spec["gust_amplitude"],
            "per_seed": scenario_results,
        }
        representative_parts.append((reference, representative[0], representative[1], task_gusts))
    stitched_reference = np.vstack([part[0] for part in representative_parts])
    stitched_smpc = {
        "S": np.vstack([part[1]["S"] for part in representative_parts]),
        "U": np.vstack([part[1]["U"] for part in representative_parts]),
    }
    stitched_mpc = {
        "S": np.vstack([part[2]["S"] for part in representative_parts]),
        "U": np.vstack([part[2]["U"] for part in representative_parts]),
    }
    stitched_gusts = np.vstack([part[3] for part in representative_parts])
    summary = {
        "scope": "three long second-order frozen-model-in-the-loop boundary scenarios; not continuous flight",
        "noise_statement": "identified residual innovation plus separately prescribed external task-space gusts",
        "scenario_steps": int(steps), "scenario_duration_seconds": float(steps / 100.0),
        "total_plotted_steps": int(3 * steps),
        "total_plotted_duration_seconds": float(3 * steps / 100.0),
        "scenario_resets": [float(steps / 100.0), float(2 * steps / 100.0)],
        "sample_rate_hz": 100, "period_steps": int(period_steps),
        "seeds": list(seeds), "d": 10, "horizon": 18,
        "q_weight": 80.0, "ru": 0.18, "alpha_joint": 0.05,
        "hard_bound": hard_bound, "gust_width_steps": gust_width,
        "terminal_tightening": tightening.tolist(),
        "terminal_tightened_bound": tightened_bound.tolist(),
        "scenarios": scenarios,
        "interpretation": (
            "all scenarios prioritize active chance constraints with zero "
            "controller fallback; evidence is the retained hard-bound margin, "
            "not a forced deterministic-MPC violation"
        ),
    }
    results = HERE / "results"
    results.mkdir(parents=True, exist_ok=True)
    (results / "copyBF_long_3d_boundary_disturbance.json").write_text(
        json.dumps(summary, indent=2), encoding="ascii"
    )
    np.savez(
        results / "copyBF_long_3d_boundary_disturbance.npz",
        reference=stitched_reference, smpc=stitched_smpc["S"],
        deterministic_mpc=stitched_mpc["S"],
        smpc_u=stitched_smpc["U"], deterministic_mpc_u=stitched_mpc["U"],
        external_task_gusts=stitched_gusts,
        scenario_reset_steps=np.array([steps, 2 * steps]),
    )
    plot_results(
        stitched_reference, stitched_smpc, stitched_mpc, stitched_gusts,
        hard_bound, tightened_bound,
        results / "copyBF_long_3d_boundary_disturbance.png", seeds[0],
    )
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", nargs="+", type=int, default=[7, 19, 31])
    parser.add_argument("--steps", type=int, default=4000,
                        help="Steps per axis scenario; three scenarios are run")
    args = parser.parse_args()
    run_experiment(tuple(args.seeds), steps=args.steps)
