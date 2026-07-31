# PredVARX-MPC

MATLAB research code for **PredVAR / PredVARX identification, oblique latent-variable models, and stochastic model predictive control (SMPC)**.

This repository studies how a reduced model

$$
z_{k+1}=Az_k+B(u_k-\bar u)+\varepsilon_{k+1},
\qquad
\hat y_k=\bar y+Pz_k,
$$

can preserve control-relevant output directions while satisfying the oblique dual-basis condition

$$
R^\top P=I.
$$

The main questions are:

- how to construct and validate the oblique latent subspace;
- how identification coordinates align with the controlled outputs;
- how prediction covariance enters chance-constrained SMPC;
- which conclusions are mathematical, implementation-level, or empirical.

> Each `copy*` directory is an independent experiment. Results from different copies are comparable only when the plant, data, seed, reference, noise, controller, constraints, and closed-loop horizon are the same.

## Current Reproducible Result

The current paper-aligned experiment is [`copyAR_crte_paper_spectral_validation_unknown_noise`](experiments/copyAR_crte_paper_spectral_validation_unknown_noise/).

It implements the CRTE experimental procedure used for manuscript Sections 3.3–3.4:

- fixed spectral surrogate on the free output complement;
- paper trace normalization
  $N_{tr}(A)=A/\max\{|\operatorname{tr}(A)|/d,10^{-5}\}$;
- task, noise, and reachability gates;
- validation-only candidate selection;
- one-step free-output residual covariance as an unknown-noise proxy;
- no access to the true sensor-noise covariance (`uses_true_Sigma_n = 0`);
- full-data refit followed by a 1200-step SMPC simulation.

The endpoint metrics are excluded. The final search uses

$$
\mu\in\{0.10,0.25,0.50,0.75\},
$$

and selects the best valid interior candidate.

### Latest copyAR run

Plant seed: `20260710`.

| Metric | Result |
|---|---:|
| Selected $\mu$ | **0.10** |
| Validation NRMSE | **0.064557** |
| Valid candidates | 4 / 4 |
| Spectral radius | 0.943072 |
| MAE $[y_1,y_2]$ | [0.06056, 0.05718] |
| RMSE $[y_1,y_2]$ | [0.08335, 0.07923] |
| QP success rate | 100% |
| Fallback / infeasible steps | 0 / 0 |
| Maximum recorded QP residual | $1.53\times10^{-12}$ |
| Constraint-active fraction | 66.29% |

![copyAR closed-loop result](experiments/copyAR_crte_paper_spectral_validation_unknown_noise/results/copyAR_crte_paper_spectral_validation_unknown_noise_fig.png)

The complete auditable outputs are stored beside the experiment:

- [`metrics.txt`](experiments/copyAR_crte_paper_spectral_validation_unknown_noise/results/copyAR_crte_paper_spectral_validation_unknown_noise_metrics.txt)
- [`data.mat`](experiments/copyAR_crte_paper_spectral_validation_unknown_noise/results/copyAR_crte_paper_spectral_validation_unknown_noise_data.mat)
- [`figure.png`](experiments/copyAR_crte_paper_spectral_validation_unknown_noise/results/copyAR_crte_paper_spectral_validation_unknown_noise_fig.png)

## Experiment Lines

### PredVARX geometry and SMPC

| Range | Purpose |
|---|---|
| `copyO`–`copyZ` | Orthogonal/oblique identification, centered SMPC, control-aware subspaces, and fair ablations |
| `copyAA`–`copyAM` | Block dual bases, mathematical corrections, closed-loop diagnostics, covariance variants, and guarantee-layer probes |
| `copyOPQSTR_unified`, `copyALL_unified` | Unified-condition comparisons |

### CRTE algorithm variants

| Path | Role | Status |
|---|---|---|
| `copyAN_crte_fixed_surrogate` → [`copyAR`](experiments/copyAR_crte_paper_spectral_validation_unknown_noise/) | Fixed spectral surrogate plus validation selection | Current manuscript experiment |
| [`copyAO`](experiments/copyAO_crte_teacher_profiled_unknown_noise/) | Profiled/min-teacher selection | Structural research variant, not the manuscript algorithm |
| [`copyAP`](experiments/copyAP_crte_multistep_task_20x4/) | Multi-step task and blocked proxy | Extension |
| [`copyAQ`](experiments/copyAQ_crte_varx_order_20x4/) | Higher-order VARX | Draft; not a completed 20x4 grid |

### Learned output-direction experiments

These copies reuse the same 1500 offline samples and add no training data. They are extensions and diagnostics, not part of the original CRTE algorithm.

| Copy | Question | 1200-step observation |
|---|---|---|
| [`copyAS`](experiments/copyAS_learned_task_anchor_smpc/) | Can the full reference trajectory define a learned task anchor? | The learned anchor did not improve tracking of the original physical outputs |
| [`copyAT`](experiments/copyAT_learned_output_directions/) | Fixed task, supervised recovery, or maximum input authority? | Fixed and supervised directions agree; high input authority alone gives worse task tracking |
| [`copyAU`](experiments/copyAU_soft_preference_output/) | Can soft output preferences define the final controlled outputs? | QP remained feasible, but original-task tracking failed |
| [`copyAV`](experiments/copyAV_hard_preference_output/) | What happens when the highest-priority physical outputs are hard locked? | Good single-seed tracking was recovered; this is not a general optimality result |

The copyAT comparison uses the same old `(u,y)`, seed, closed-loop disturbance, controller, constraints, 1200-step horizon, and five 240-step reference segments. Performance is evaluated on the original physical outputs `y1,y2`.

![copyAT output-direction comparison](experiments/copyAT_learned_output_directions/results/copyAT_learned_output_directions_fig.png)

## Repository Layout

```text
PredVARX-MPC/
├── main/                         # historical baseline
├── experiments/                  # independent copy* experiments
│   ├── copyO_oblique/ ... copyAM_tracked_cov_only/
│   ├── copyAN_crte_fixed_surrogate/
│   ├── copyAO_crte_teacher_profiled_unknown_noise/
│   ├── copyAP_crte_multistep_task_20x4/
│   ├── copyAQ_crte_varx_order_20x4/
│   ├── copyAR_crte_paper_spectral_validation_unknown_noise/
│   ├── copyAS_learned_task_anchor_smpc/
│   ├── copyAT_learned_output_directions/
│   ├── copyAU_soft_preference_output/
│   ├── copyAV_hard_preference_output/
│   ├── copyOPQSTR_unified/
│   └── copyALL_unified/
├── tests/                        # canonical matlab.unittest suite
├── docs/version-summary/         # illustrated version history
├── archive/                      # non-current historical copies
└── README.md
```

`copyAN_closed_loop_audit_grade` is an empty placeholder and does not represent a completed experiment.

## Requirements

- MATLAB R2024a or newer;
- Optimization Toolbox (`quadprog`);
- Git LFS for `.mat`, `.fig`, `.png`, `.jpg`, and `.pdf` artifacts.

```bash
git lfs install
git clone https://github.com/xyzhhhh2077/PredVARX-MPC.git
cd PredVARX-MPC
```

Without Git LFS, large result files are checked out as pointer files.

## Run copyAR

From MATLAB:

```matlab
repo = pwd;
copyAR = fullfile(repo, 'experiments', ...
    'copyAR_crte_paper_spectral_validation_unknown_noise');

cd(copyAR)
addpath(copyAR, fullfile(copyAR, 'tests'))
test_crte_paper_spectral_varx
copyAR_crte_paper_spectral_validation_unknown_noise
```

The full 1200-step SMPC run can take several minutes. Results are written to the experiment's `results/` directory.

## Run Canonical Tests

```matlab
cd('tests')
run_all_tests
```

The root suite covers:

- `tPredvarxIdentifyOblique`: oblique dual-basis identities;
- `tCenteredSmpcStep`: centered prediction and chance-constraint tightening;
- `tControlAwareSubspace`: preservation of controlled axes in the reduced coordinates.

Each experiment may also contain focused contract tests in its own `tests/` directory.

## Evidence Boundaries

1. Mathematical proofs, code contracts, and simulation observations are reported separately.
2. `QP success = 100%` and `fallback = 0` show that one trajectory used the primary SMPC path; they do not prove recursive feasibility or stability.
3. The residual covariance is an engineering proxy under unknown sensor noise, not the true sensor-noise covariance.
4. Oracle quantities may be used only for simulation comparisons, not as inputs to an unknown-noise method.
5. A direction with high input authority is easy to actuate; it is not automatically the output required by the task.
6. A single-seed result does not establish generalization, global optimality, or universal performance improvement.

## Documentation

| Document | Scope |
|---|---|
| [`VERSION_EVOLUTION_AND_METRICS.md`](VERSION_EVOLUTION_AND_METRICS.md) | Evolution and metrics for copies O–Z |
| [`LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md`](LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md) | Legacy letter reuse and current-directory mapping |
| [`PREDVARX_MPC_VERSION_SUMMARY.md`](docs/version-summary/PREDVARX_MPC_VERSION_SUMMARY.md) | Illustrated version summary |
| [`archive/README.md`](archive/README.md) | Archive policy and historical contents |

## Reproducibility Policy

- Large MATLAB data and figures are tracked with Git LFS.
- A new research attempt gets a new `copy*` directory; established experiment identities are not silently overwritten.
- Metrics produced under different protocols are not merged into a single ranking.
- Every reported run should identify its copy, seed, training data, noise convention, closed-loop length, and completion status.
