# copyAC_open_problems_trial

## Purpose

New experiment copy for **open pure-math opinions that were explained but not algorithmically closed** in copyAA/copyAB:

| Opinion | Trial method | Default in this runner |
|---|---|---|
| 5 Sigma_eps denom | switch `t2` / `ml` / `ols` | **ols** |
| 6 Cross cov | inject `Sigma_zo` into chance `Sigma_y` | **ON** |
| 7 Sigma_obs support | `declared_shape` / `residual_support` / `additive` | **additive** |
| 8 Fallback guarantee | soft recovery ladder before uncertified hold | **ON** |
| 9 Recursive feas/stability | optional terminal cost | **ON** (soft only) |
| 10 Input-conditional IVR | residualize-on-U free IVR | **ON** (candidate) |

Does **not** modify `main/`, `copyX`, `copyAA`, `copyAB`.

## Honest non-claims

1. Soft recovery uses **inflated risk / shorter N / bound-only QP** — not the original chance certificate.
2. Cross-term `Sigma_y` is an experimental inclusion; Boole allocation is not re-proved under cross blocks.
3. `residual_support` / `additive` Sigma_obs are engineering trial objects, not theorems equating them to `Cov(o)`.
4. Terminal cost is a soft regularizer only — **not** recursive feasibility or closed-loop stability.
5. Input residualize IVR is a **candidate** free-direction scheme, not proved input-conditional PredVARX.

## Files

- `copyAC_open_problems_trial.m` — full closed-loop trial runner
- `centered_smpc_step.m` — optional cross cov + terminal cost
- `split_control_free_ivr_varx.m` — multi-denom + residualize (from copyAB)
- `lib/soft_recovery_smpc.m` — L2 recovery ladder
- `lib/build_sigma_obs_trial.m` — Sigma_obs trial modes
- `lib/cross_cov_diagnostics.m`, `fallback_certify_step.m`, `sigma_obs_support_diag.m`
- `tests/test_copyAC_open_problem_trials.m`

## Run

```matlab
cd('.../experiments/copyAC_open_problems_trial');
addpath(pwd); addpath('lib'); addpath('tests');
test_copyAC_open_problem_trials;
run('copyAC_open_problems_trial.m');
```


## Fair comparison vs copyAB

See `results/copyAB_vs_copyAC_fair_compare.md` (and `.txt`) for a same-seed re-run of copyAB defaults vs this trial stack.
