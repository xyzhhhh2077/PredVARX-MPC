# Opinion 7 patch report

status: done

## Opinion (source)

PDF §7「残差自由度与协方差支撑空间」:

- Residual scale uses denominator `(p-ell)`, so residual is understood as living in a `(p-ell)`-dimensional unresolved subspace.
- Writing `Sigma_obs = sigma^2 I_p` (or any full-rank full-space matrix) as if it were `Cov(o)` is inconsistent with that support.
- Canonical residual form: `Sigma_o = Q_o Sigma_xi Q_o'` (isotropic: `sigma^2 Q_o Q_o'`).
- Full-space declared sensor noise `Sigma_n` is a **different** object and must not be silently identified with `Cov((I-PR')y)`.

## Design implemented

Diagnostic only (does **not** change the online SMPC covariance law).

1. New library `lib/sigma_obs_support_diag.m`
   - Inputs: `o_window` (p×L projected residuals), `Sigma_n_declared` (p×p), `p`, `ell`
   - Outputs: `rank_o`, `rank_n`, `trace_o`, `trace_on`, `note`
   - `rank_o` / `trace_o` from empirical `Cov(o)`
   - `rank_n` from declared `Sigma_n`
   - `trace_on` from the **same online proxy** as the runner:
     `Sigma_on = max(sigma_hat^2,1e-8) * shape(Sigma_n)`
     with `sigma_hat` using residual DOF `(p-ell)`
   - `note` explicitly states: online agent is **not** `Cov(o)`

2. Runner `copyAB_opinion_fixes.m` (minimal intrusion)
   - `addpath(.../lib)`
   - After each online `Sigma_obs` update, optionally call the diagnostic
   - Store last result in `last_sigma_obs_diag`
   - Save as `out.sigma_obs_support_diag`
   - Wrapped in `try/catch` so missing/failing diagnostics never break the loop
   - **Does not alter** `model.Sigma_obs` (still scaled declared shape)

3. Test `tests/test_opinion07_sigma_obs.m`
   - Synthetic `o = Q_o xi` ⇒ `rank_o = p-ell`
   - Heteroscedastic PD `Sigma_n` ⇒ `rank_n = p`
   - Checks `trace_on` matches manual proxy
   - Asserts note contains NOT-identity language
   - Guards bad input dimensions

## Explicit claim discipline

| Object | Support / rank | Meaning |
|--------|----------------|---------|
| `Sigma_n` declared | typically rank `p` | calibrated full-space sensor noise |
| empirical `Cov(o)` | ≈ `p-ell` | projected residual (noise + leak) |
| online `Sigma_obs` | rank `p` | sensor-floor **proxy** = scale × shape(`Sigma_n`) |

**The online agent is not the identity `Cov(o)`.**  
It is a conservative full-space floor used inside `Sigma_y = P Sigma_z P' + Sigma_obs`. Rank/trace diagnostics expose the support mismatch; they do not force `Sigma_obs := Q_o ... Q_o'`.

## Files

| Path | Action |
|------|--------|
| `lib/sigma_obs_support_diag.m` | created |
| `tests/test_opinion07_sigma_obs.m` | created |
| `copyAB_opinion_fixes.m` | path + optional last-diag record |
| `patches/opinion_7.md` | this report |

Baseline `main/`, `copyX`, `copyAA` untouched.

## Verification

```matlab
cd('.../experiments/copyAB_opinion_fixes');
addpath(pwd); addpath(fullfile(pwd,'lib')); addpath(fullfile(pwd,'tests'));
test_opinion07_sigma_obs;
```

Expected: `PASS opinion07 Sigma_obs support: rank_o=p-ell, rank_n=p, ... note: ... NOT Cov(o) ...`
