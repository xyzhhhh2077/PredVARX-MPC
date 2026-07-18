# Opinion 6 patch (copyAI): cross-cov into chance Sigma_y

**Status: PASS (derivation + test + fair CL run). No Boole re-proof.**

## Scope

Only `experiments/copyAI_crosscov_chance/`. Default SMPC path unchanged
(`opt.use_cross_cov` default false).

## Derivation checked

```text
y_c = P z + o,   z = R' y_c,   R' P = I
Sigma_y_full = P Sz P' + P Szo + Szo' P' + So
Sigma_y_drop = P Sz P' + So
R' o = 0  =/=>  Szo = 0
```

Tracked H=e1,e2 under split R=[E,...]: `Szo R = 0` ⇒ tracked
`hq'(P Szo+Szo'P') hq = 0` with empirical So ⇒ var_ratio=1 expected.

## Code

- `lib/cross_cov_diagnostics.m` — full/drop + optional H variances
- `lib/sigma_y_chance_blocks.m` — block builder
- `centered_smpc_step.m` / `build_chance_rows.m` — `use_cross_cov`
- `run_crosscov_config.m`, `copyAI_crosscov_chance.m`
- `tests/test_opinion06_cross_cov_chance.m`
- Fix in runner: `Ipr = eye(p) - Phat*Rhat'` (was missing transpose)

## Test stdout (MATLAB R2024a)

```text
PASS opinion06 chance: full_rel=4.393e-16 drop=1.506e+00 ||Szo||=1.056e+00 vr=[1.286 1.006] OFF/ON b_diff=3.039e+00
TEST_DONE
```

## Closed-loop stdout (excerpt)

```text
smoke OFF cover=1.000 MAE1=0.053 ||Szo||=4.036e-01 drop=3.004e-02 full=4.230e-16
smoke ON  cover=1.000 MAE1=0.053 var_ratio=[1.000 1.000] qp=1.000
[1/6] AI_cross_OFF_s94001
  cover=0.9940 MAE=0.586/0.447 ... qp=0.494 fail=253 act=0.392 ||Szo||=4.036e-01 drop=3.004e-02 vr=[1.000 1.000]
[2/6] AI_cross_OFF_s94002
  cover=0.9860 MAE=0.660/0.453 ... qp=0.452 fail=274 act=0.377
[3/6] AI_cross_OFF_s94003
  cover=0.9760 MAE=0.589/0.435 ... qp=0.484 fail=258 act=0.389
[4/6] AI_cross_ON_s94001
  cover=0.9960 MAE=0.638/0.558 ... qp=0.008 fail=496 act=0.010
[5/6] AI_cross_ON_s94002
  cover=0.9860 MAE=0.733/0.609 ... qp=0.000 fail=500 act=0.000
[6/6] AI_cross_ON_s94003
  cover=0.9760 MAE=0.684/0.577 ... qp=0.044 fail=478 act=0.055
COPYAI_DONE
```

## Aggregate

| mode | cover | MAE1 | qp_success | active | \|\|Szo\|\|_F | drop_rel | vr_tracked |
|---|---:|---:|---:|---:|---:|---:|---:|
| OFF | 0.9853 | 0.6119 | 0.4767 | 0.3857 | 0.404 | 0.030 | 1.0 |
| ON  | 0.9860 | 0.6852 | 0.0173 | 0.0216 | 0.404 | 0.030 | 1.0 |

## Non-claims

- **Cannot** claim Boole was re-proved under cross terms.
- Naive ON + offline Szo + online Sigma_obs proxy **hurts** QP feasibility;
  need PSD-consistent joint law before recommending ON as default.
- Artifacts: `results/copyAI_crosscov_chance_metrics.csv`, `results/REPORT.md`,
  `results/copyAI_crosscov_chance_data.mat`.
