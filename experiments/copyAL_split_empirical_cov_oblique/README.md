# copyAL_split_empirical_cov_oblique

## Purpose

Correct the **theory/code misalignment** of copyAA's free dual:

- **User design**: fix `P=[E, Nperp*V]`, `R=[E, Nperp*W]`, then choose free dual `W` under the **empirical total free-coordinate covariance** of `Yp = Nperp'*yc`. True sensor noise `Sigma_n` is **unknown** and must **not** enter the identifier.
- **copyAA bug relative to that design**: passes oracle `Sigma_n` and uses `Sigma_perp = Nperp'*Sigma_n*Nperp`.
- The plant generator may use a hidden noise law to create synthetic observations, but
  `copyAL` does not pass the true law or its scale to the identifier/controller. Its
  `Sigma_obs` is initialized and updated from empirical projection residuals only; it is
  an empirical total-residual proxy, not true sensor-noise covariance.

This copy implements the empirical metric. It does **not** claim:

1. sensor-noise optimum theorem, or  
2. PredVAR minimum-innovation-covariance theorem.

## Free dual (alpha = 1)

```text
C_emp = (Yp*Yp')/(T-1) + ridge
W     = C_emp^{-1} V (V' C_emp^{-1} V)^{-1}
```

solves `min tr(W' C_emp W) s.t. W'V = I` under the **empirical-total-covariance metric**.

`oblique_alpha` in `[0,1]`: `W = V + alpha*(Wfull-V)`.  
`alpha=1` is the conditional optimum under that metric; `alpha<1` is only interpolation.

## Files

| File | Role |
|---|---|
| `split_empirical_cov_ivr_varx.m` | identifier (no `Sigma_n` arg) |
| `copyAL_split_empirical_cov_oblique.m` | closed-loop runner |
| `compare_copyAL_vs_copyAA.m` | same-seed plant compare |
| `tests/test_split_empirical_cov_ivr_varx.m` | geometry + metric + no-`Sigma_n` gates |
| `centered_smpc_step.m` / `smooth_noise_profile.m` | local helpers (copied from AA) |

## Run

```matlab
cd('.../experiments/copyAL_split_empirical_cov_oblique')
test_split_empirical_cov_ivr_varx
copyAL_split_empirical_cov_oblique
% or full compare:
compare_copyAL_vs_copyAA
```

On Chinese Windows, prefer ASCII-path runner + `char([20195 30721])` for `代码`.

## Claim boundary

| Object | Meaning |
|---|---|
| free emp-cov objective | `tr(W' C_reg W)` with free-block empirical total cov |
| cover | warm-horizon QP-success **time fraction proxy only** |
| plant `Sigma_n` | used only to simulate sensors; never passed to identifier |

Baseline `main/`, `copyX`, `copyAA` are **not** modified.
