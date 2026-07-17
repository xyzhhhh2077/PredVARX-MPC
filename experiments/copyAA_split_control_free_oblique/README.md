# copyAA_split_control_free_oblique

## Purpose

This is a new self-contained successor experiment to `copyX_control_aware_oblique`. It does not modify `main/` or copyX.

The algorithm separates two roles that copyX mixed together:

1. **Tracked/control directions** are hard-fixed in both loading and extraction matrices.
2. **Free latent directions** are learned by IVR and receive a true sensor-noise-covariance-weighted dual extractor.
3. The two blocks are merged into one latent coordinate, and the complete coupled VARX model is refitted before SMPC.

**Boundary (opinion 10):** IVR on free directions still uses only projected output lags (`Ylag/Ycur` in `null(E')`). This yields a **predictable-in-output-history** free subspace, **not** a proved input-conditional optimal PredVARX subspace. Inputs enter only at the subsequent VARX step `Phi=[zc;ur]`.

## Mathematics

Let

```text
P = [E, Nperp*V]
R = [E, Nperp*W]
```

where `E` contains the tracked standard axes and `V` is learned by IVR in `null(E')`.

For a declared positive-definite sensor-noise covariance `Sigma_n`, define

```text
Sigma_perp = Nperp' * Sigma_n * Nperp
W = Sigma_perp^{-1} V (V' Sigma_perp^{-1} V)^{-1}
```

This minimizes the true free-subspace read-noise objective under `W'*V=I`.

The construction guarantees:

```text
R'*P = I
P*R'*E = E
E'*P*R' = E'
```

Thus tracked directions have both right-invariant coverage and value preservation. Only the free block can be oblique.

## Isotropic boundary

When `Sigma_n = sigma^2*I`, the strict optimum satisfies `W=V`, hence `R=P`. This is a required boundary, not a failure. A nontrivial statistically meaningful oblique free extractor requires heteroscedastic or correlated sensor noise.

## Files

- `split_control_free_ivr_varx.m` — split identifier and free-noise-optimal extractor.
- `copyAA_split_control_free_oblique.m` — full process-like closed-loop experiment.
- `centered_smpc_step.m` — centered chance-constrained SMPC copied from copyX.
- `smooth_noise_profile.m` — smooth nonstationary noise envelope.
- `tests/test_split_control_free_ivr_varx.m` — TDD regression for isotropic and heteroscedastic cases.
- `results/` — generated MAT/PNG remain untracked; the small metrics text is committed as evidence.

## TDD evidence

RED, before production implementation:

```text
'split_control_free_ivr_varx' could not be recognized
```

GREEN focused test:

```text
PASS split extractor: iso ||R-P||=9.065e-16;
het free-oblique=3.020e+00;
noise objective 0.429293 -> 0.146887
(improvement 2.824e-01)
```

## Run

```matlab
addpath('experiments/copyAA_split_control_free_oblique/tests');
test_split_control_free_ivr_varx;
run('experiments/copyAA_split_control_free_oblique/copyAA_split_control_free_oblique.m');
```

## Interpretation boundary

The full experiment deliberately declares a heteroscedastic/correlated `Sigma_n` so the free-space oblique dual has a genuine minimum-read-noise interpretation. The isotropic case remains covered by the regression test and correctly reduces to `R=P`.

## Full closed-loop verification

```text
tracked right/left/column errors = 0 / 0 / 0
dual error = 5.745e-17
free oblique norm = 1.04275
true free-noise objective = 1.15842e-2 -> 8.49066e-3
MAE = [0.0656, 0.0613]
RMSE = [0.0921, 0.0860]
QP success = 100%, fallbacks = 0
max chance-row residual = -4.441e-16
upper/absolute violations = 0 / 0
```

The PNG was generated and decoded structurally. The local vision tool did not receive the file correctly, so visual inspection is recorded as inconclusive rather than claimed.
