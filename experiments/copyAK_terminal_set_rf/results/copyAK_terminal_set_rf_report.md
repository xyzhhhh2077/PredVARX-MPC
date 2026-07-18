# copyAK: Opinion 9 terminal set / RF numerical probes

**Not theorem closures. No chance-constraint RF. No closed-loop stability theorem.**

## Minimal RF condition chain

| ID | Condition | Status in this probe |
|---|---|---|
| C1 | rho(Ahat)<1 (Schur latent model) | PROBED-PASS (rho=0.942272) |
| C2 | Pterm=dlyap(A',Qf), residual small | PROBED-PASS (res=4.743e-14) |
| C3 | V-decrease on free z+=A z | PROBED-PASS (ratio=1.000) |
| C4 | terminal law keeps Xf +invariant under constraints | **UNPROVED** |
| C5 | stage feasible => z_N in exact Xf | **UNPROVED** (optional box rows only) |
| C6 | chance-constraint RF under noise | **UNPROVED** |

### What is implemented (engineering probe)

1. Terminal cost (opt, default OFF): `Pterm≈dlyap(A',Qf)`, add `z_N'Pterm z_N`.
2. Terminal set (opt, default OFF): spectral-box outer approx of
   `z_N' Pterm z_N ≤ α_term` as linear QP rows (not exact ellipsoid SOCP).
3. `α_term` calibrated from free unit-ball V quantile + offline latent energy.
4. Soft recovery remains engineering ladder (not original-α cert).

### alpha_term used

`alpha_term = 447.969`

## Offline Lyapunov probe

| quantity | value |
|---|---:|
| rho(Ahat) | 0.942272 |
| is_schur | 1 |
| dlyap_ok | 1 |
| lyap residual F | 4.743e-14 |
| V-decrease OK ratio | 1.000 |
| ctrb rank | 5 |
| left geometry err | 0.000e+00 |

## Closed-loop OFF vs ON (multi-seed)

Modes: OFF = no terminal; COST = terminal cost only; SET = cost + box terminal set.

| family | n | mean cover | mean MAE1 | mean MAE2 | mean qp | mean soft | mean unc | mean V_term_med |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| MED_OFF | 2 | 0.9975 | 0.3393 | 0.2587 | 0.555 | 0.445 | 0.000 | NaN |
| MED_COST | 2 | 0.9975 | 0.3394 | 0.2597 | 0.556 | 0.444 | 0.000 | 5.262 |
| MED_SET | 2 | 0.9975 | 0.3394 | 0.2597 | 0.556 | 0.444 | 0.000 | 5.262 |
| TGT_OFF | 1 | 0.9950 | 0.4599 | 0.3761 | 0.435 | 0.565 | 0.000 | NaN |
| TGT_COST | 1 | 0.9950 | 0.4605 | 0.3761 | 0.432 | 0.568 | 0.000 | 4.546 |
| TGT_SET | 1 | 0.9950 | 0.4605 | 0.3761 | 0.432 | 0.568 | 0.000 | 4.546 |

### Per-run detail

| name | cover | MAE1 | MAE2 | qp | soft | unc | lyap_res | Vdec | V_term_med |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| MED_OFF_y0.70_s91001 | 0.9950 | 0.3413 | 0.2602 | 0.595 | 0.405 | 0.000 | 4.74e-14 | 1.000 | NaN |
| MED_COST_y0.70_s91001 | 0.9950 | 0.3417 | 0.2616 | 0.595 | 0.405 | 0.000 | 4.74e-14 | 1.000 | 8.456 |
| MED_SET_y0.70_s91001 | 0.9950 | 0.3417 | 0.2616 | 0.595 | 0.405 | 0.000 | 4.74e-14 | 1.000 | 8.456 |
| MED_OFF_y0.70_s91002 | 1.0000 | 0.3373 | 0.2572 | 0.515 | 0.485 | 0.000 | 4.74e-14 | 1.000 | NaN |
| MED_COST_y0.70_s91002 | 1.0000 | 0.3371 | 0.2579 | 0.517 | 0.482 | 0.000 | 4.74e-14 | 1.000 | 2.068 |
| MED_SET_y0.70_s91002 | 1.0000 | 0.3371 | 0.2579 | 0.517 | 0.482 | 0.000 | 4.74e-14 | 1.000 | 2.068 |
| TGT_OFF_y0.50_s92001 | 0.9950 | 0.4599 | 0.3761 | 0.435 | 0.565 | 0.000 | 4.74e-14 | 1.000 | NaN |
| TGT_COST_y0.50_s92001 | 0.9950 | 0.4605 | 0.3761 | 0.432 | 0.568 | 0.000 | 4.74e-14 | 1.000 | 4.546 |
| TGT_SET_y0.50_s92001 | 0.9950 | 0.4605 | 0.3761 | 0.432 | 0.568 | 0.000 | 4.74e-14 | 1.000 | 4.546 |

## Non-claims (mandatory)

1. No recursive feasibility theorem under chance constraints / noise.
2. No closed-loop stability theorem (with or without terminal cost/set).
3. Free-dynamics Schur + V-decrease ≠ constrained Xf positive invariance.
4. Spectral-box rows are a **coarse outer approximation**, not exact ellipsoid Xf.
5. Soft recovery ≠ original joint-risk certificate.
6. Empirical cover/MAE shifts under terminal ON are probes, not RF proofs.

## Takeaways

1. If C1-C3 pass, free latent dynamics admit a discrete Lyapunov function (numeric).
2. Terminal cost/set are optional regularizers/constraints default OFF.
3. Compare OFF/COST/SET at medium and tight y_max for empirical effect size.
4. Still open: true Xf construction under chance rows; kappa terminal law; CL stability.
