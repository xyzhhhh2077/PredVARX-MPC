# copyV_iterative_ivr

This experiment is copied from the committed `copyU_smooth_noise_smpc` baseline and changes only the reduced-subspace estimator.

- copyU: one SVD in the non-tracked complement + one VARX OLS fit.
- copyV: iterative IVR in the non-tracked complement, followed by the same VARX and centered SMPC/QP pipeline.
- tracked quality axes remain exactly retained.
- copyU is not modified.

Run from the repository root:

```matlab
run('experiments/copyV_iterative_ivr/copyV_iterative_ivr.m')
```

The run writes `results/copyV_iterative_ivr_data.mat`, `results/copyV_iterative_ivr_metrics.txt`, and `results/copyV_iterative_ivr_fig.png`.

The estimator records `stats.ivr_iter`, `stats.ivr_trace`, and `stats.ivr_subspace_delta`. Compare copyV and copyU using identical random seed and controller settings. Do not claim improvement from IVR based on QP success alone; compare prediction error, orthogonal-noise leakage, chance-constraint violations, tracking error, and control effort.
