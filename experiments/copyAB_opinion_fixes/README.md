# copyAB_opinion_fixes

## Purpose

New experiment copy derived from `copyAA_split_control_free_oblique`.  
**Does not modify** `main/`, `copyX`, or `copyAA`.

Goal: try a concrete code/doc fix for each of the 10 pure-math review opinions, with real MATLAB tests where possible.

## Design logic (before edits)

| Opinion | Fix type | Owned artifacts |
|---|---|---|
| 1 Geometry | assert left/right/column | `split_*` + `tests/test_opinion01_*` |
| 2 True noise metric | declared Σ_n free-block optimum flags | `split_*` + `tests/test_opinion02_*` |
| 3 M vs τ | comments + small numeric check | `split_*` comments + `tests/test_opinion03_*` |
| 4 Stop ≠ converge | export `ivr_trace` fields | `split_*` + `tests/test_opinion04_*` |
| 5 Sigma denom | primary zero-mean Gaussian ML scale; OLS DOF diagnostic | `split_*` + `tests/test_opinion05_*` |
| 6 Cross cov | diagnostic library | `lib/cross_cov_diagnostics.m` |
| 7 Sigma_obs support | diagnostic library | `lib/sigma_obs_support_diag.m` |
| 8 Fallback cert | runner catch path | `copyAB_opinion_fixes.m` |
| 9 Terminal cost | optional flag default off | `centered_smpc_step.m` |
| 10 Input residual IVR | optional flag default off | `split_*` |

## Hard boundaries

- Geometry guarantees are **reconstruction-layer** only, not closed-loop tracking RMSE=0.
- Free dual is optimal only for **declared** `Sigma_n` with fixed free loading `V`.
- IVR remains output-history predictable unless optional residualize-on-`U` is enabled; even then it is a candidate, not a proved optimum.
- Fallback without re-checking chance rows is **uncertified**.
- **Strict convexity ≠ recursive feasibility ≠ stability.** Enabling
  `opt.use_terminal_cost=true` (default **false**) only adds an optional
  discrete-Lyapunov quadratic terminal weight on \(z_N\); it is **not** a
  terminal-set / Lyapunov certificate and does **not** prove recursive
  feasibility or closed-loop stability.

## Run

```matlab
cd('.../experiments/copyAB_opinion_fixes');
addpath(pwd); addpath(fullfile(pwd,'lib')); addpath(fullfile(pwd,'tests'));
test_copyAB_core_geometry;
% then per-opinion tests if present
run('copyAB_opinion_fixes.m');
```

## Patch reports

See `patches/opinion_1.md` … `opinion_10.md` written by subagents.
