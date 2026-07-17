# Fair comparison: copyAB vs copyAC

Same plant/seed/parameter family. Real MATLAB re-runs in one batch.

| metric | copyAB (defaults) | copyAC (trials ON) |
|---|---:|---:|
| `tracked_left_error` | 0.000000000000e+00 | 0.000000000000e+00 |
| `tracked_right_error` | 0.000000000000e+00 | 0.000000000000e+00 |
| `tracked_column_error` | 0.000000000000e+00 | 0.000000000000e+00 |
| `dual_error` | 5.744789948873e-17 | 3.311397765341e-16 |
| `free_oblique_norm` | 1.042752953743e+00 | 1.042265973139e+00 |
| `free_noise_improvement` | 3.093576438203e-03 | 3.093466702187e-03 |
| `MAE` | 0.065618457626 0.061325142093 | 0.058029832742 0.056244207183 |
| `RMSE` | 0.092136024871 0.086041727267 | 0.078527016026 0.079262049091 |
| `upper_violation_rate` | 0.000000000000 0.000000000000 | 0.000000000000 0.000000000000 |
| `qp_certified_rate` | 1.000000000000 | 1.000000000000 |
| `qp_success_rate` | 1.000000000000 | NA |
| `uncertified_fallback_rate` | 0.000000000000 | 0.000000000000 |
| `soft_recovery_rate` | NA | 0.000000000000 |
| `max_recorded_qp_constraint` | -4.440892098501e-16 | -1.274536032270e-13 |
| `constraint_active_rate` | 0.664761904762 | 0.045714285714 |
| `avg_costJ` | 61.604110048318 | 64.946948151589 |
| `cross_drop_rel_err` | NA | 2.937249345567e-02 |
| `cross_Sigma_zo_fro` | NA | 3.952071804491e-01 |
| `sigma_eps_mode` | NA | ols |
| `use_cross_cov` | NA | 1 |
| `sigma_obs_mode` | NA | additive |
| `input_residualize` | NA | 1 |
| `use_terminal_cost` | NA | 1 |

## Setup

- Shared: `rng(20260710)`, `n=6,m=3,p=30,ell=5`, `T_cl=1200`, `N=18`, `y_max=2`, heteroscedastic/correlated declared `Sigma_n`.
- **copyAB**: split free dual on declared `Sigma_n`; IVR residualize **OFF**; no cross term in chance `Sigma_y`; online `Sigma_obs` = residual-scale x declared shape; terminal **OFF**; soft recovery **OFF**.
- **copyAC**: `Sigma_eps` **ols**; residualize **ON**; cross-cov **ON**; `Sigma_obs` **additive**; terminal **ON**; soft recovery **ON**.

## How to read

1. Geometry errors near 0 on both: reconstruction-layer left/right/column/dual.
2. RMSE/MAE differences = closed-loop empirical effect of the trial stack, not a stability theorem.
3. If soft_recovery_rate=0 and uncertified=0, primary QP never failed; soft path was not stress-tested.
4. Active-rate change can come from different Sigma_y / Sigma_obs / Sigma_eps, not only from better control.

## Non-claims

- Does not prove recursive feasibility, closed-loop stability, or input-conditional PredVARX optimality.
- Does not select a statistically optimal Sigma_eps denominator by coverage calibration.

## Numeric deltas (AC − AB)

| quantity | copyAB | copyAC | delta (AC-AB) |
|---|---:|---:|---:|
| MAE y1 | 0.06562 | 0.05803 | **-0.00759** |
| MAE y2 | 0.06133 | 0.05624 | **-0.00508** |
| RMSE y1 | 0.09214 | 0.07853 | **-0.01361** |
| RMSE y2 | 0.08604 | 0.07926 | **-0.00678** |
| avg cost J | 61.604 | 64.947 | **+3.343** |
| constraint active rate | 0.6648 | 0.0457 | **-0.619** |
| qp certified rate | 1.000 | 1.000 | 0 |
| upper violation rate | 0 / 0 | 0 / 0 | 0 |
| free noise improvement | 3.0936e-3 | 3.0935e-3 | ~0 |
| dual error | 5.7e-17 | 3.3e-16 | ~0 |

### Interpretation (empirical, not theorems)

1. **Tracking improved** on both quality outputs under the AC trial stack (lower MAE/RMSE).
2. **Cost rose slightly** (`avg J` +3.3): more conservative / different Sigma path and terminal regularizer can trade cost for tracking.
3. **Constraint activity collapsed** (0.665 → 0.046): AC uses larger/additive observation noise and cross terms, so Boole tightening changes; fewer near-active chance rows does **not** by itself prove better safety theory.
4. **Geometry unchanged**: both keep reconstruction-layer left/right/column ~ 0.
5. **Soft recovery unused** (rate 0): this fair seed never stressed QP infeasibility; soft path remains unvalidated under pressure.
6. **Do not claim** recursive feasibility, stability, or proved input-conditional PredVARX from these deltas.
