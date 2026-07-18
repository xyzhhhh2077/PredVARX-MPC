# copyAI Opinion 6: cross-covariance in chance Sigma_y

Date: 2026-07-19. MATLAB real run. **No commit.**

## What was derived

With centered decomposition `y_c = P z + o`, `z = R' y_c`, dual `R' P = I`:

```text
Sigma_y_full = P Sigma_z P' + P Sigma_zo + Sigma_zo' P' + Sigma_o
Sigma_y_drop = P Sigma_z P' + Sigma_o          % default SMPC path
```

- Samplewise `R' o = 0`  **⇏**  `Sigma_zo := Cov(z,o) = 0`.
- Chance direction: compare `hq' Sigma_y_full hq` vs `hq' Sigma_y_drop hq`.
- For **tracked** rows `hq = e_i` that are columns of `R` (split construction `R = [E, ...]`),
  `Szo * R = 0` sample-algebraically ⇒ `e_i'(P Szo + Szo' P') e_i = 0` when using empirical `Sigma_o`.
  Hence offline `var_ratio = full/drop = 1` on tracked H is **expected**, not a bug.
- Frobenius drop on the full matrix can still be >0 (`drop_cross_rel_err ≈ 3e-2` here) because free / residual directions carry the cross mass.

## What was implemented

| Artifact | Role |
|---|---|
| `lib/cross_cov_diagnostics.m` | z,o,Σz,Σo,Σzo, drop/full recon, optional H variances |
| `lib/sigma_y_chance_blocks.m` | pure full/drop Σy + `hq'Σy hq` |
| `centered_smpc_step.m` | `opt.use_cross_cov` default **false**; ON adds `P Szo + Szo'P'` |
| `lib/build_chance_rows.m` | same flag |
| `run_crosscov_config.m` | closed-loop one config (recovery=none) |
| `copyAI_crosscov_chance.m` | fair multi-seed OFF vs ON |
| `tests/test_opinion06_cross_cov_chance.m` | focused gates |

Default path unchanged when `use_cross_cov` false / missing `Sigma_zo`.

## What was verified (real MATLAB)

### Focused test stdout

```text
PASS opinion06 chance: full_rel=4.393e-16 drop=1.506e+00 ||Szo||=1.056e+00 vr=[1.286 1.006] OFF/ON b_diff=3.039e+00
TEST_DONE
```

Gates: `R'*o~0`, `full_recon_rel_err~0`, `||Szo||_F` can be >0, drop can be >0; ON changes `b_ch` on generic H.

### Offline (identified split dual, same for all seeds)

| metric | value |
|---|---:|
| `\|\|Sigma_zo\|\|_F` | 4.036e-01 |
| `drop_cross_rel_err` | 3.004e-02 |
| `full_recon_rel_err` | 4.23e-16 |
| `cross_block_fro` | 5.75e-01 |
| tracked `var_ratio` (H=e1,e2, So=Σo) | 1.000 / 1.000 |
| dual_err / left_err | ~0 / 0 |

### Closed-loop fair OFF vs ON (y_max=0.55, α=0.10, T_cl=500, N=18, seeds 94001–94003)

| mode | mean cover | MAE1 | MAE2 | RMSE1 | qp_success | active_rate |
|---|---:|---:|---:|---:|---:|---:|
| cross_OFF | 0.9853 | 0.6119 | 0.4449 | 0.7431 | **0.4767** | **0.3857** |
| cross_ON  | 0.9860 | 0.6852 | 0.5814 | 0.7920 | **0.0173** | **0.0216** |

Per-seed CSV: `results/copyAI_crosscov_chance_metrics.csv`.

## Interpretation (honest)

1. **Gap is real in matrix Frobenius**: `||Σzo||_F≈0.40`, `drop_rel≈3%`, full recon machine-zero.
2. **Tracked chance axes do not feel empirical cross quadratic form** under split `R=[E,…]` (var_ratio=1).
3. **Naive ON injection hurts primary QP**: mean `qp_success` 0.48 → 0.017.
   Likely cause: online uses `Sigma_y = PΣzP' + Σ_obs + (PΣzo+·)` where `Σ_obs` is a **proxy**, not empirical `Σo`. The cross block that is PSD-compatible with `Σo` need not keep `Σ_obs+cross` PSD; eig-floor then distorts the tightening landscape and feasibility collapses. MAE worsens because most steps fall to `u_mean` uncertified fallback.
4. Cover stays high (~0.98) mainly because plant/refs often stay under y_max even with open-loopish u_mean — **not** evidence that ON is safer.

## Honest non-claims

- **Does NOT** re-prove Boole / union-bound risk allocation under the full cross-covariance law.
- Does NOT claim rolling online `Σzo` estimation.
- Does NOT claim `Σ_obs ≡ Cov(o)`.
- Does NOT claim joint empirical cover certifies `α_joint`.
- Injecting cross is an **engineering trial**; default remains drop-cross until a PSD-consistent joint law is designed.

## Setup

- plant: n=6, p=30, ell=5, tracked=[1 2], split free dual + declared Σn
- y_max=0.55, alpha_joint=0.10, T_cl=500, N=18, seeds=[94001 94002 94003]
- recovery: none (primary QP only)
- offline data seed rng(20260710); shared y_off/u_off for all jobs
