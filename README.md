# PredVARX-MPC

MATLAB / Python research code for **PredVAR / PredVARX identification**, **oblique latent-variable models**, and **chance-constrained SMPC**.

Reduced model (typical form):

$$
z_{k+1}=A z_k+B(u_k-\bar u)+\varepsilon_{k+1},
\qquad
\hat y_k=\bar y+P z_k,
\qquad
R^\top P=I.
$$

**Each `experiments/copy*` directory is an independent experiment.**  
Compare closed-loop numbers only when plant, data, seed, reference, noise, controller, constraints, and horizon match.

---

## What this repo is (and is not)

| Is | Is not |
|---|---|
| Offline identification + closed-loop SMPC studies | A flight stack or ROS driver |
| Frozen real-data models used as plants (model-in-the-loop) | Automatic claims of outdoor / hardware validation |
| Separate copies for each algorithmic change | One “latest always overwrites” folder |
| Explicit claim boundaries in each copy README | Oracle noise treated as a theorem |

---

## Start here (current endpoints)

### 1. Paper-aligned CRTE — `copyAR`

Manuscript-style spectral validation under **unknown** sensor noise  
(`uses_true_Sigma_n = 0`).

- Path: [`experiments/copyAR_crte_paper_spectral_validation_unknown_noise/`](experiments/copyAR_crte_paper_spectral_validation_unknown_noise/)
- Role: fixed spectral free-complement surrogate, paper \(N_{tr}\), gates, validation-only selection, residual proxy, then SMPC sim
- Not the same as `copyAO` (min-teacher / profiled structural variant)

![copyAR closed-loop](experiments/copyAR_crte_paper_spectral_validation_unknown_noise/results/copyAR_crte_paper_spectral_validation_unknown_noise_fig.png)

### 2. Pelican real-data SMPC — `copyBI` → `copyBJ`

Built on the public [Waterloo WAVELab AscTec Pelican dataset](https://github.com/wavelab/pelican_dataset)  
(54 indoor Vicon flights @ 100 Hz; flights 1–36 ID, 37–54 held out).

| Copy | Role |
|---|---|
| [`copyBI`](experiments/copyBI_pelican_probabilistic_boundary_advantage/) | Frozen model, six-face near-boundary SMPC + visual audit (baseline event) |
| [`copyBJ`](experiments/copyBJ_pelican_online_covariance_smpc/) | Same freeze; **only** online composite-innovation \(\Sigma_\varepsilon\) adaptation |

**copyBI** synchronized run (200-frame GIF):

![copyBI SMPC](experiments/copyBI_pelican_probabilistic_boundary_advantage/results_smpc/copyBI_smpc_summary_200frames.gif)

**copyBJ** online \(\Sigma_\varepsilon\) run (β = 0.8 main trajectory):

![copyBJ online SMPC](experiments/copyBJ_pelican_online_covariance_smpc/results_smpc/copyBJ_smpc_summary_200frames.gif)

Pelican contract (both):

- \(y\in\mathbb{R}^{10}\): position (3), Euler (3), measured motor speeds (4)  
- \(u\in\mathbb{R}^{4}\): commanded motor speeds  
- latent \(z\in\mathbb{R}^{5}\): 3 fixed standardized `xyz` axes + 2 free directions  
- MPC state: second-order \([z_k;z_{k-1}]\in\mathbb{R}^{10}\)

copyBJ one-line claim: **offline ID + online composite \(\Sigma_\varepsilon\)**; not online \(E/A/B/P/R\); not process/measurement/wind split; β-shrink is an experimental regularizer, not a recursive-feasibility theorem. Details: [`copyBJ/README.md`](experiments/copyBJ_pelican_online_covariance_smpc/README.md).

---

## Experiment map (basic)

### Geometry / SMPC core

| Range | Purpose |
|---|---|
| `copyO`–`copyZ` | Orthogonal / oblique ID, centered SMPC, control-aware subspaces, fair ablations |
| `copyAA`–`copyAM` | Dual-basis fixes, diagnostics, covariance variants, guarantee-layer probes |
| `copyOPQSTR_unified`, `copyALL_unified` | Unified-condition comparisons |

### CRTE line

| Copy | Role |
|---|---|
| `copyAN` → **`copyAR`** | Fixed spectral surrogate + validation selection (**paper path**) |
| `copyAO` | Profiled / min-teacher (research variant, not manuscript algorithm) |
| `copyAP` | Multi-step task + blocked residual proxy |
| `copyAQ` | Higher-order VARX (draft) |

### Output preference / direction trials

Same 1500 offline samples where stated; **extensions**, not the original CRTE algorithm.

| Copy | Question (short) |
|---|---|
| `copyAS` | Learned task anchor from full reference? |
| `copyAT` | Fixed vs supervised vs high input-authority directions |
| `copyAU` / `copyAV` | Soft vs hard output preferences |
| `copyAW` | BOPTEST pilot (soft preference) |
| `copyAX` | ControlGym CDR soft preference |
| `copyAY` | Pelican soft preference / real-flight ID side line |

### Pelican position SMPC line

`copyBA` → `copyBB` → … → **`copyBI`** → **`copyBJ`**  
(ID → fixed anchor → boundary stress → 3-D → wind → u-cap → visual audit → online covariance)

---

## Layout

```text
PredVARX-MPC/
├── main/                 # historical baseline
├── experiments/          # independent copy* experiments
├── tests/                # root matlab.unittest suite
├── docs/version-summary/ # illustrated history
├── archive/              # non-current snapshots
└── README.md
```

`copyAN_closed_loop_audit_grade` is an empty placeholder (not a finished run).

---

## Requirements

- MATLAB R2024a+ (Optimization Toolbox / `quadprog` for MATLAB copies)
- Python 3.10+ for Pelican Python copies (`copyBI`, `copyBJ`, …) as documented in each folder
- **Git LFS** for `.mat`, `.fig`, `.png`, `.jpg`, `.pdf`, `.gif`

```bash
git lfs install
git clone https://github.com/xyzhhhh2077/PredVARX-MPC.git
cd PredVARX-MPC
```

Without LFS, large artifacts checkout as pointer files.

---

## Quick runs

### copyAR (MATLAB)

```matlab
repo = pwd;
copyAR = fullfile(repo, 'experiments', ...
    'copyAR_crte_paper_spectral_validation_unknown_noise');
cd(copyAR)
addpath(copyAR, fullfile(copyAR, 'tests'))
test_crte_paper_spectral_varx
copyAR_crte_paper_spectral_validation_unknown_noise
```

### copyBJ (Python)

```bash
python experiments/copyBJ_pelican_online_covariance_smpc/run_online_covariance_smpc.py
python -m pytest tests/test_copybj_pelican_online_covariance_smpc.py -q
# optional BI-style pack from stored beta=0.8 trajectory (no re-solve):
cd experiments/copyBJ_pelican_online_covariance_smpc
python generate_smpc_figures.py
python generate_smpc_summary_gif.py
python verify_smpc_summary_gif.py
```

### Root unit tests (MATLAB)

```matlab
cd('tests')
run_all_tests
```

Covers oblique dual-basis, centered SMPC step, and control-aware subspace contracts. Many copies also ship local tests.

---

## Evidence boundaries (global)

1. **Proof ≠ code contract ≠ single-seed simulation.** Report them separately.  
2. `QP success = 100%` and `fallback = 0` mean that trajectory used the main path — **not** recursive feasibility or stability.  
3. Residual / innovation covariances under unknown sensor noise are **engineering proxies**, not true \(\Sigma_n\) or separated \(Q_w,R_v\).  
4. Oracle quantities are for comparison only; unknown-noise methods must not consume them.  
5. High input authority ≠ the output the task requires.  
6. One seed ≠ generalization or global optimality.  
7. Physical wind injected in the plant is **not** automatically a calibrated term inside chance tightening unless the copy states that explicitly.

---

## More documentation

| Document | Scope |
|---|---|
| [`VERSION_EVOLUTION_AND_METRICS.md`](VERSION_EVOLUTION_AND_METRICS.md) | Copies O–Z evolution / metrics |
| [`LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md`](LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md) | Legacy letter reuse map |
| [`docs/version-summary/`](docs/version-summary/) | Illustrated version summary |
| [`archive/README.md`](archive/README.md) | Archive policy |
| Per-copy `README.md` | Local protocol, seeds, claim boundary |

---

## Reproducibility policy

- New idea → **new `copy*`**. Do not silently overwrite an established identity.  
- Large binaries via **Git LFS**.  
- Do not merge metrics from mismatched protocols into one ranking.  
- Every reported run should name: **copy, seed, data, noise convention, horizon, completion status**.
