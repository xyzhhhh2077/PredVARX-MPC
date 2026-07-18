# copyAD_closure_ladder

## Purpose

Sequential **empirical** closure attempts for remaining open pure-math opinions, in priority order:

| Priority | Opinion | What this copy does |
|---|---|---|
| P0 | 8 soft path | Tight `y_max` stress: soft OFF vs ON (and +cross) |
| P1 | 5 | `Sigma_eps` t2/ml/ols coverage proxy comparison |
| P1 | 7 | `Sigma_obs` declared_shape / residual_support / additive |
| P2 | 6 | cross-cov ON vs OFF |
| P2 | 10 | input residualize ON vs OFF |

Does **not** modify `main/`, `copyX`, `copyAA`, `copyAB`, `copyAC`.

## Non-claims

- Soft recovery ≠ original chance certificate.
- Empirical joint cover ≠ calibrated residual DOF theorem.
- Sigma_obs modes ≠ proved Cov(o).
- Cross ON ≠ re-proved Boole allocation.
- Residualize ON ≠ input-conditional PredVARX optimum.
- Opinion 9 stability/recursive feasibility **not** closed here.

## Run

```matlab
cd('.../experiments/copyAD_closure_ladder');
addpath(pwd); addpath('lib');
run('copyAD_closure_ladder.m');
```

Outputs:

- `results/copyAD_closure_ladder_metrics.csv`
- `results/copyAD_closure_ladder_report.md`
