# Opinion 5 patch report

status: done

## Goal

Expose multiple residual-covariance denominators for `Sigma_eps` while using the conditional zero-mean Gaussian ML denominator `N_res` as the primary return value. The OLS residual-DOF scale remains diagnostic only.

## Math

After the coupled VARX fit, residual matrix `Eps` has `N = T-1` columns (one lag lost). With Gram `G = Eps*Eps'`:

| Name | Formula | Denominator | Role |
|---|---|---|---|
| `Sigma_eps` (primary) | `G / max(N,1)` | `N = T-1` | Conditional zero-mean Gaussian ML / population scaling |
| `Sigma_eps_ml` | `G / max(N,1)` | `N` | Gaussian ML / population MLE |
| `Sigma_eps_ols` | `G / max(N-(ell+m),1)` | residual DOF after regressing `ell+m` params | OLS residual covariance |

All three are symmetrized `(S+S')/2`.

## Code change

File: `split_control_free_ivr_varx.m` only (Opinion 5 section).

- Compute `Nres`, `Gram_eps`, then three denoms and matrices.
- **Return value** `Sigma_eps` uses the primary zero-mean Gaussian ML scale.
- `stats` records:
  - `Sigma_eps`, `Sigma_eps_ml`, `Sigma_eps_ols`
  - `N_residual`, `Sigma_eps_denom_primary`, `Sigma_eps_denom_ml`, `Sigma_eps_denom_ols`

## Downstream safety

- `copyAB_opinion_fixes.m` and `centered_smpc_step.m` consume the 5th return `Sigma_eps` → conditional zero-mean Gaussian ML scale.
- ML/OLS variants are opt-in via `stats` only; no caller switch was added.

## Test

`tests/test_opinion05_sigma_denoms.m`:

1. Required multi-denom fields exist.
2. Denoms match `N-1`, `N`, `max(N-(ell+m),1)` with `N=T-1`.
3. Return `Sigma_eps == stats.Sigma_eps`.
4. Gram consistency across the three denoms.
5. All three matrices PSD (symmetric + `eig >= -1e-8`).

### Real MATLAB run (2026-07-18)

```
PASS opinion05 sigma denoms: N=399; denoms t2/ml/ols=398/399/392; eigmin t2/ml/ols=3.349e-03/3.341e-03/3.401e-03
PASS split extractor: iso ||R-P||=9.065e-16; ... (core geometry still green)
ALL_PASS
```

## Boundary

Reporting multiple denoms does **not** make the OLS scale valid under the present ridge/data-dependent-subspace pipeline; the primary return uses the conditional zero-mean Gaussian ML scale, while OLS DOF is retained only as a clearly labeled diagnostic.
