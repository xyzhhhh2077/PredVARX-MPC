# copyAF: Opinion 9 probes + P1++ coverage

**Not theorem closures.**

## Opinion 9 numerical probe

| quantity | value |
|---|---:|
| rho(Ahat) | 0.942272 |
| is_schur | 1 |
| dlyap_ok | 1 |
| lyap residual F | 4.743e-14 |
| V-decrease OK ratio | 1.000 |
| ctrb rank | 5 |
| left geometry err | 0.000e+00 |

CL terminal OFF vs ON (y_max=0.70, soft ON):

| | cover | MAE | soft_rate |
|---|---:|---:|---:|
| term OFF | 1.0000 | 0.3204/0.2415 | 0.334 |
| term ON | 1.0000 | 0.3210/0.2424 | 0.334 |

Non-claim: Schur + Lyapunov decrease on **unconstrained free dynamics** does **not** prove recursive feasibility or chance-constrained closed-loop stability.

## P1++ Sigma_eps multi-seed tight limits

| y_max | eps | n | mean cover | std cover | mean MAE1 | mean soft_rate |
|---:|---|---:|---:|---:|---:|---:|
| 0.50 | t2 | 5 | 0.9956 | 0.0017 | 0.4528 | 0.576 |
| 0.50 | ml | 5 | 0.9952 | 0.0018 | 0.4516 | 0.567 |
| 0.50 | ols | 5 | 0.9968 | 0.0027 | 0.4636 | 0.657 |
| 0.45 | t2 | 5 | 0.9952 | 0.0023 | 0.4879 | 0.648 |
| 0.45 | ml | 5 | 0.9948 | 0.0023 | 0.4866 | 0.641 |
| 0.45 | ols | 5 | 0.9880 | 0.0125 | 0.4939 | 0.707 |

## P1++ Sigma_obs multi-seed y_max=0.50

| obs | n | mean cover | std cover | mean MAE1 | mean soft_step_viol |
|---|---:|---:|---:|---:|---:|
| declared_shape | 5 | 0.9956 | 0.0017 | 0.4528 | 0.0076 |
| residual_support | 5 | 0.9860 | 0.0063 | 0.4319 | 0.0202 |
| additive | 5 | 0.9908 | 0.0033 | 0.4436 | 0.0161 |

## Takeaways

1. Opinion 9: if Ahat is Schur and dlyap residual tiny, free dynamics admit a Lyapunov function — still no RF under constraints/chance/noise.
2. P1++: multi-seed means/std of cover at y_max=0.50/0.45; use only if cover separates.
3. Soft ON keeps uncert near zero; rankings remain empirical.
4. Still open: true RF/stability theorems; original-alpha soft cert; DOF optimality; Cov(o) identity; input-conditional PredVARX.
