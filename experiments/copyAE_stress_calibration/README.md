# copyAE_stress_calibration

## Purpose

Recommended follow-ups after copyAD:

1. **P0+** multi-seed soft recovery with stage breakdown (`risk_inflate` / `short_horizon` / `bound_only`) and soft-step violation rate.
2. **P1+** tighter `y_max` so empirical joint cover can drop below 1, enabling denom/obs comparisons.

Does not modify earlier copies.

## Run

```matlab
cd('.../experiments/copyAE_stress_calibration');
addpath(pwd); addpath('lib');
run('copyAE_stress_calibration.m');
```

## Outputs

- `results/copyAE_stress_calibration_metrics.csv`
- `results/copyAE_stress_calibration_report.md`
