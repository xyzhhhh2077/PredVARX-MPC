# copyAE stress calibration (P0+ / P1+)

Multi-seed soft-stage breakdown and tight-limit coverage probes. **Not theorem closures.**

`T_cl=600`. Offline data `rng(20260710)`. Soft stages: risk_inflate / short_horizon / bound_only.

## P0+ Soft multi-seed (y_max=0.55)

| seed | soft | primary_fail | soft_ok | uncert | succ|fail | cover | soft_step_viol | risk | shortN | bound |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 71001 | 0 | 351 | 0 | 351 | 0.000 | 0.9850 | NaN | 0 | 0 | 0 |
| 71001 | 1 | 352 | 352 | 0 | 1.000 | 0.9967 | 0.0057 | 63 | 289 | 0 |
| 71002 | 0 | 365 | 0 | 365 | 0.000 | 0.9833 | NaN | 0 | 0 | 0 |
| 71002 | 1 | 366 | 366 | 0 | 1.000 | 0.9950 | 0.0055 | 83 | 283 | 0 |
| 71003 | 0 | 330 | 0 | 330 | 0.000 | 0.9617 | NaN | 0 | 0 | 0 |
| 71003 | 1 | 330 | 330 | 0 | 1.000 | 0.9983 | 0.0000 | 37 | 293 | 0 |
| 71004 | 0 | 356 | 0 | 356 | 0.000 | 0.9383 | NaN | 0 | 0 | 0 |
| 71004 | 1 | 356 | 356 | 0 | 1.000 | 0.9983 | 0.0028 | 73 | 283 | 0 |
| 71005 | 0 | 358 | 0 | 358 | 0.000 | 0.9417 | NaN | 0 | 0 | 0 |
| 71005 | 1 | 358 | 358 | 0 | 1.000 | 0.9917 | 0.0140 | 44 | 314 | 0 |

**Aggregate soft-ON:** mean succ|fail = 1.000; fraction seeds with uncert=0: 1.00 (5/5).

Non-claim: stage mix shows *how* soft recovers, not original-alpha validity.

## P1+ Sigma_eps under tight y_max (soft ON)

| config | y_max | cover | MAE1 | MAE2 | soft_rate | uncert_rate |
|---|---:|---:|---:|---:|---:|---:|
| P1_eps_t2_y0.80 | 0.80 | 1.0000 | 0.2777 | 0.1908 | 0.337 | 0.000 |
| P1_eps_ml_y0.80 | 0.80 | 1.0000 | 0.2776 | 0.1901 | 0.325 | 0.000 |
| P1_eps_ols_y0.80 | 0.80 | 1.0000 | 0.2794 | 0.1953 | 0.438 | 0.000 |
| P1_eps_t2_y0.65 | 0.65 | 0.9967 | 0.3462 | 0.2722 | 0.515 | 0.000 |
| P1_eps_ml_y0.65 | 0.65 | 0.9967 | 0.3453 | 0.2714 | 0.510 | 0.000 |
| P1_eps_ols_y0.65 | 0.65 | 0.9967 | 0.3529 | 0.2786 | 0.560 | 0.000 |
| P1_eps_t2_y0.55 | 0.55 | 0.9983 | 0.4177 | 0.3418 | 0.555 | 0.000 |
| P1_eps_ml_y0.55 | 0.55 | 0.9983 | 0.4164 | 0.3408 | 0.552 | 0.000 |
| P1_eps_ols_y0.55 | 0.55 | 1.0000 | 0.4280 | 0.3522 | 0.613 | 0.000 |

## P1+ Sigma_obs at y_max=0.65

| config | cover | MAE1 | MAE2 | soft_rate |
|---|---:|---:|---:|---:|
| P1_obs_declared_shape_y0.65 | 0.9983 | 0.3422 | 0.2635 | 0.513 |
| P1_obs_residual_support_y0.65 | 0.9933 | 0.3217 | 0.2398 | 0.498 |
| P1_obs_additive_y0.65 | 0.9950 | 0.3346 | 0.2559 | 0.505 |

## Takeaways

1. P0+: multi-seed confirms soft path reliability under stress when enabled.
2. Stage histogram (risk/short/bound) documents recovery mix.
3. P1+: if cover separates across denoms/modes at tight limits, use as empirical ranking only.
4. Still open: original-alpha certificate; DOF theorem; Cov(o) identity; opinion 9 stability.
