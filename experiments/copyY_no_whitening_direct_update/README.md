# copyY_no_whitening_direct_update

## Purpose

This is a new experiment derived from `copyX_control_aware_oblique`. The existing copyX and `main/` baseline remain unchanged.

The version removes the normalization operation shown in the referenced Algorithm-1 line,

$$
Y^*=YUD^{-1/2},
$$

and performs every iterative predictable-direction update directly in the raw centered complement coordinates

$$
Y_\perp=N_\perp^T(Y-\bar Y).
$$

If the current direction matrix is $V_i$, the new dominant eigenspace directly replaces it:

$$
V_{i+1}\leftarrow V_{\mathrm{new}}.
$$

There is no whitening, variance normalization, normalized-coordinate iteration, or de-normalization step.

## What remains identical to copyX

- Fixed tracked axes and exact tracked-output coverage.
- Complement-space predictable-direction iteration.
- Covariance-weighted oblique extractor and $\alpha=0.02$ regularization.
- Re-extraction of latent coordinates followed by centered VARX regression.
- Centered chance-constrained SMPC and online covariance-window update.
- Plant, random seed, reference, constraints, and plotting layout.

## Theoretical boundary

This is an engineering direct-update variant, not a faithful implementation of the whitening/de-normalization stages of Mo--Qin Algorithm 1. It should be compared with strict normalized PredVAR and with copyX. Geometric dual-basis checks alone do not establish statistical equivalence.

## Verification

Run the focused test and then the experiment:

```matlab
addpath('experiments/copyY_no_whitening_direct_update/tests');
test_control_aware_direct_update_varx;
run('experiments/copyY_no_whitening_direct_update/copyY_no_whitening_direct_update.m');
```

The persisted metrics explicitly report:

- `normalization_applied 0`;
- `update_rule direct_eigenspace_replacement`;
- dual-basis and tracked-coverage errors;
- closed-loop tracking, chance-constraint, feasibility, and control metrics.
