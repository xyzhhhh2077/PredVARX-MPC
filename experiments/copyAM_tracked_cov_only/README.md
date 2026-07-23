# copyAM_tracked_cov_only

## Purpose

Correct the **theory/code misalignment** of copyAA's free dual:

- **User design**: fix `P=[E, Nperp*V]`, `R=[E, Nperp*W]`, then choose free dual `W` under the **empirical total free-coordinate covariance** of `Yp = Nperp'*yc`. True sensor noise `Sigma_n` is **unknown** and must **not** enter the identifier.
- **copyAA bug relative to that design**: passes oracle `Sigma_n` and uses `Sigma_perp = Nperp'*Sigma_n*Nperp`.
- The plant generator may use a hidden noise law to create synthetic observations, but
  `copyAM` does not pass the true law or its scale to the identifier/controller. Its
  `Sigma_obs` is initialized and updated from empirical projection residuals only; it is
  an empirical total-residual proxy, not true sensor-noise covariance.

This copy implements the empirical metric. It does **not** claim:

1. sensor-noise optimum theorem, or  
2. PredVAR minimum-innovation-covariance theorem.

## Tracked-only covariance variant

Unlike AL, this copy does not add `model.Sigma_obs` to the chance-constraint
variance. The chance rows act only on `tracked`, and the split geometry gives
`E_c'*(I-P*R')*y_c = 0` exactly. Therefore the QP uses only the corresponding
block of `P*Sigma_z*P'`. The empirical projection-residual covariance remains
saved as a diagnostic, not as an asserted sensor-noise covariance or calibrated
probability-safe inflation.

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
| `copyAM_tracked_cov_only.m` | closed-loop runner |
| `compare_copyAM_vs_copyAA.m` | same-seed plant compare |
| `tests/test_split_empirical_cov_ivr_varx.m` | geometry + metric + no-`Sigma_n` gates |
| `centered_smpc_step.m` / `smooth_noise_profile.m` | local helpers (copied from AA) |

## Run

```matlab
cd('.../experiments/copyAM_tracked_cov_only')
test_split_empirical_cov_ivr_varx
copyAM_tracked_cov_only
% or full compare:
compare_copyAM_vs_copyAA
```

On Chinese Windows, prefer ASCII-path runner + `char([20195 30721])` for `代码`.

## Claim boundary

| Object | Meaning |
|---|---|
| free emp-cov objective | `tr(W' C_reg W)` with free-block empirical total cov |
| cover | warm-horizon QP-success **time fraction proxy only** |
| plant `Sigma_n` | used only to simulate sensors; never passed to identifier |

Baseline `main/`, `copyX`, `copyAA` are **not** modified.

## Proposition-1 4-piece dual-basis completion (CRTE draft Sec 3.1)

As of 2026-07-24, `split_empirical_cov_ivr_varx.m` additionally exports the
complete dual basis `Pbar, Rbar` and four identity residuals:

| Stat field        | Identity              | Tolerance    |
|---|---|---|
| `dual_error_RP`        | `R' P = I_ell`       | `< 1e-10`    |
| `dual_error_RPbar`     | `R' Pbar = 0`        | `< 1e-10`    |
| `dual_error_RbarP`     | `Rbar' P = 0`        | `< 1e-10`    |
| `dual_error_RbarPbar`  | `Rbar' Pbar = I_{p-ell}` | `< 1e-10` |
| `dual_basis_completion`| all four residuals < 1e-10 | logical 1 |
| `Pi_idempotency_err`   | `||P R' P R' - P R'||_fro` | `< 1e-10` |

Numpy-verified: `R = [E, Nperp*Wfree]` is automatically dual-consistent
(‖R - DualMat[:,1:ell]‖ ≈ 1e-16), so it is not overwritten; only `Pbar` and
`Rbar` are newly constructed via `Pbar = null(R'); DualMat = inv([P Pbar]')`.

Construction: `Pbar = null(R'); DualMat = inv([P Pbar]'); Rbar = DualMat(:, ell+1:end)`.
This is the same path as `copyR_moqin_oblique/predvarx_identify_moqin.m`.
The MATLAB `DualMat = inv(Qfull')` (transpose inside `inv`) is mandatory;
`inv(Qfull)` without the transpose returns the wrong factor.

Live evidence (seed 20260723, p=6, q=2, ell=4, N=600, alpha=1):

```
dual_error_RP    = 2.78e-17
dual_error_RPbar = 1.15e-16
dual_error_RbarP   = 3.40e-17
dual_error_RbarPbar = 1.12e-16
dual_basis_completion = 1
Pi_idempotency_err = 1.76e-16
```

TDD `test_split_empirical_cov_ivr_varx` still PASS:
`dual=3.86e-16; left/right=0; free-oblique=3.09; emp-obj 5.31→4.82; pert dJ=2.96e-08; no Sigma_n`.

## Audit-bug disclosure (CRTE draft Sec 4.3)

The two post-freeze source defects in CRTE draft Sec 4.3 (FWL view copy +
gamma-I shift that breaks generalized eigenvectors) are **not** present in
this MATLAB identifier:

- This experiment does not build the FWL sample-space `B_T, A_T` blocks; it
  only uses the empirical free-coordinate covariance for the free dual, so
  the in-place `current_x = x[start:stop]` view-defect cannot arise here.
- It does not shift the metric by `+ gamma*I`; the free dual solves
  `(C_emp + ridge*I) W = V (V' (C_emp + ridge*I) V)^{-1}` where ridge is tiny
  (`1e-8 * scale + 1e-12`) — a Tikhonov regularization applied to `C_emp`
  directly, not a `gamma*I` change to a generalized eigenproblem.

The CRTE numerical main table (Sec 5) is therefore **not** covered by this
MATLAB copy. Its `CODE_AUDIT_NOTES.md` rerun is a separate Python-side task.
