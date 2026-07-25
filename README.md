# PredVARX-MPC

Public MATLAB experiment suite for **PredVAR / PredVARX + SMPC** copies.

- **GitHub**: https://github.com/xyzhhhh2077/PredVARX-MPC
- **Visibility**: public
- **Default branch**: `main`
- **Large files**: experiment `.mat` / `.png` / `.fig` tracked with **Git LFS**

> Principle: each `copy*` directory is a self-contained version (code, helpers, local tests, results, notes).  
> Do **not** mix different copies into one algorithm. Archive letters may reuse Q/S/Z names with **different** meanings from the active suite.

---

## Clone

```bash
git lfs install
git clone https://github.com/xyzhhhh2077/PredVARX-MPC.git
cd PredVARX-MPC
```

Without LFS you will only get pointer files for mats/figures.

---

## Repository layout

```text
PredVARX-MPC/
├── main/                         # historical baseline (copyS_nosoft)
├── experiments/                  # active copy suite (preferred entry)
│   ├── copyO_oblique/ ...
│   ├── copyAO_crte_teacher_profiled_unknown_noise/
│   ├── copyAP_crte_multistep_task_20x4/
│   └── ...
├── archive/                      # non-active historical MATLAB versions
│   ├── early_snapshot_2026-07-09/
│   ├── organized_legacy_2026-07-10/
│   └── routeC_delta_inc/
├── tests/                        # root matlab.unittest suite
├── docs/version-summary/         # illustrated version report
├── VERSION_EVOLUTION_AND_METRICS.md
├── LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md
└── README.md                     # this file
```

---

## Related docs

| Doc | Content |
|---|---|
| [`VERSION_EVOLUTION_AND_METRICS.md`](VERSION_EVOLUTION_AND_METRICS.md) | Active O–Z evolution, fairness rules, metrics |
| [`LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md`](LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md) | Old workspace letter reuse vs current Git names |
| [`docs/version-summary/PREDVARX_MPC_VERSION_SUMMARY.md`](docs/version-summary/PREDVARX_MPC_VERSION_SUMMARY.md) | Illustrated per-version report with PNGs |
| [`archive/README.md`](archive/README.md) | What is archived and why |

---

## Active experiments (`experiments/`)

Current preferred work lives here. Counts below are from the published tree (code + tracked results).

### A. Identification / geometry / control spine

| Directory | Role |
|---|---|
| `main/` | Historical baseline: `copyS_nosoft` + `predvarx_identify` |
| `copyO_oblique` | Oblique dual-basis identity diagnostics |
| `copyP_centered_smpc` | Centered absolute-prediction SMPC (full-order control check) |
| `copyQ_control_aware` | Low-order control-aware orthogonal subspace |
| `copyR_moqin_oblique` | Mo–Qin Algorithm-1 style whitened realization baseline |
| `copyT_process_lv_smpc` | High-dim process-like plant + online cov scale |
| `copyU_smooth_noise_smpc` | Smooth time-varying noise stress of copyT-style loop |
| `copyV_iterative_ivr` | Tracked-complement iterative IVR (orthogonal) |
| `copyW_fair_identifier_compare` | Fair multi-identifier ablation, fixed plant/controller |
| `copyX_control_aware_oblique` | Control-aware mild oblique reader $R_\alpha$ (no whitening) |
| `copyY_no_whitening_direct_update` | Explicit “copyX fact label”; not a new algorithm |
| `copyZ_strict_whitened_copyX_plant` | Whitening ablation on copyX plant (not full Alg.1 claim) |

### B. Opinion / closure / guarantee probes

| Directory | Role |
|---|---|
| `copyAA_split_control_free_oblique` | Split control/free oblique construction |
| `copyAB_opinion_fixes` | Opinion hard-fix bundle + unit tests |
| `copyAC_open_problems_trial` | Open-problem trial vs AB fair compare |
| `copyAD_closure_ladder` | Closure ladder diagnostics |
| `copyAE_stress_calibration` | Stress / calibration probe |
| `copyAF_op9_and_p1pp` | Opinion-9 / Prop-1 related probe |
| `copyAG_original_alpha_cert` | Original-alpha certificate experiment |
| `copyAH_multistep_alpha_backup` | Multi-step alpha backup after primary fail |
| `copyAI_crosscov_chance` | Cross-covariance chance probe |
| `copyAJ_safety_filter_cert` | Safety-filter certificate probe |
| `copyAK_terminal_set_rf` | Terminal-set / recursive-feasibility probe |
| `copyAL_split_empirical_cov_oblique` | Empirical-cov split oblique variant |
| `copyAM_tracked_cov_only` | Tracked-only covariance chance path |
| `copyAN_crte_fixed_surrogate` | CRTE fixed spectral surrogate in free complement |
| `copyAN_closed_loop_audit_grade` | Placeholder (empty run; branch-name continuity only) |

### C. Unified compare packs

| Directory | Role |
|---|---|
| `copyOPQSTR_unified` | Same plant/noise/SMPC; restore O/P/Q/R/S/T algorithm differences |
| `copyALL_unified` | Broader unified multi-copy compare pack + individual figures |

### D. CRTE profiled-teacher line (current research edge)

| Directory | Role | Status |
|---|---|---|
| `copyAO_crte_teacher_profiled_unknown_noise` | Complete profiled teacher, unknown-noise proxy, strict FWL SVD support, single full run | **Main CRTE structure result** |
| `copyAP_crte_multistep_task_20x4` | Multi-step task stack $t+1:t+H$, $\mu=1$, blocked forward noise proxy | Structure tests + smoke + one full $H=3$ confirm |
| `copyAQ_crte_varx_order_20x4` | VARX order / companion draft for $s>1$ ablation | Draft cores/tests; not a finished 20×4 campaign |

**CRTE reading order**

1. Geometry / FWL / teacher structure → `copyAO`
2. Multi-step task + blocked proxy contract → `copyAP`
3. Higher-order latent AR draft → `copyAQ` (incomplete)

Do not merge AO/AP/AQ numbers into one “best method” claim without a paired protocol.

---

## Archive (`archive/`)

Non-active historical MATLAB. Prefer `experiments/` for new work.

| Path | Content |
|---|---|
| `archive/early_snapshot_2026-07-09/` | Early independent snapshot: main + O/P/Q |
| `archive/organized_legacy_2026-07-10/` | Themed legacy pack (01–10): old Q/S/Z/SM/D/oracle letters |
| `archive/routeC_delta_inc/` | v10 Δ-tracking acceleration line (C/E/F/G/H/I/P1–P4) |

Legacy letter reuse is documented in `LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md`.  
Example: archive `copyQ_*` is **not** the same algorithm as active `copyQ_control_aware`.

---

## How to run (MATLAB R2024a+)

From repo root:

```matlab
% historical baseline
run('main/copyS_nosoft.m')

% active examples
run('experiments/copyP_centered_smpc/copyP_centered_smpc.m')
run('experiments/copyX_control_aware_oblique/copyX_control_aware_oblique.m')

% CRTE main structure run
cd('experiments/copyAO_crte_teacher_profiled_unknown_noise')
test_crte_profiled_teacher_unknown_noise
copyAO_crte_teacher_profiled_unknown_noise

% CRTE multistep task tests
cd('../copyAP_crte_multistep_task_20x4')
runtests('tests/test_copyAP_multistep_task.m')
```

Each active copy writes snapshots under its own `results/` (`.mat` / `.png` / metrics text).

---

## Root tests

```matlab
cd('tests')
run_all_tests
```

| Test class | Coverage |
|---|---|
| `tPredvarxIdentifyOblique` | Oblique dual-basis identities |
| `tCenteredSmpcStep` | Centered prediction + Boole tightening |
| `tControlAwareSubspace` | Tracked-axis coverage under reduction |

Per-copy opinion/diagnostic tests live in `experiments/copy*/tests/`.

---

## Data & figures policy

- **Tracked with Git LFS**: `*.mat`, `*.png`, `*.fig` under experiment/archive results
- **Not tracked**: literature PDF caches, stray `experiments/*.pdf`, local editor helpers
- Result figures are **run evidence**, not automatically “paper-final” claims
- Cross-copy MAE ranking is invalid unless plant, noise, horizon, and controller protocol match

---

## Quick map: what to open first

| Goal | Open |
|---|---|
| Understand O–Z evolution | `VERSION_EVOLUTION_AND_METRICS.md` |
| Fair OPQSTR compare | `experiments/copyOPQSTR_unified/` |
| Current CRTE teacher structure | `experiments/copyAO_crte_teacher_profiled_unknown_noise/` |
| Multi-step task extension | `experiments/copyAP_crte_multistep_task_20x4/` |
| Old letter archaeology | `LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md` + `archive/` |
| Illustrated summary | `docs/version-summary/PREDVARX_MPC_VERSION_SUMMARY.md` |

---

## Notes

- MATLAB path often needs `addpath` to the chosen copy directory (and sometimes `tests/`).
- Optimization Toolbox `quadprog` is required for SMPC loops in most copies.
- Parallel pools are optional; some later scripts support multi-process launchers instead of large parpools.
- When reporting results, always name the **exact copy directory** and whether the figure/mat is smoke, single-seed, or full campaign.
