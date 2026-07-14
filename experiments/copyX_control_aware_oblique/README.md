# copyX_control_aware_oblique

## Purpose

This experiment tests the user's question: should the control-aware `copyV` model use a general oblique projector `P*R'` instead of the orthogonal special case `P*P'`?

It is a separate copy. `copyV`, `copyR`, and `main/` are not modified.

## Design

1. Keep exactly the same control-aware loading subspace as `copyV`:
   `P = [E_tracked, Nperp*V]`.
2. Replace `R=P` by a covariance-weighted dual extractor
   `Rfull = Sigma_y^{-1} P (P' Sigma_y^{-1} P)^{-1}`.
3. Regularize the obliqueness with
   `R_alpha = P + alpha*(Rfull-P)`.
   This preserves `R_alpha'*P=I` for every alpha.
4. Recompute `z=R_alpha'*(y-y_mean)` and refit the final `A,B` in that coordinate.
5. Use the existing centered controller, which already applies `R'` for extraction and `P` for output prediction.

## Verification

Run:

```matlab
addpath('experiments/copyX_control_aware_oblique/tests');
test_control_aware_oblique_ivr_varx;
run('experiments/copyX_control_aware_oblique/screen_oblique_alpha.m');
run('experiments/copyX_control_aware_oblique/copyX_control_aware_oblique.m');
```

The focused test verifies `R'*P=I`, `P~=R`, non-symmetric `P*R'`, exact tracked-output coverage, dimensions, and positive-semidefinite innovation covariance.

## Held-out alpha screen

On the same synthetic process and seed, held-out one-step output RMSE was minimized at `alpha=0` (`R=P`). Increasing obliqueness monotonically worsened reconstruction and prediction:

| alpha | output RMSE | tracked RMSE | reconstruction residual |
|---:|---:|---:|---:|
| 0.00 | 0.1176 | 0.0803 | 0.117 |
| 0.02 | 0.1176 | 0.0806 | 0.117 |
| 0.10 | 0.1202 | 0.0882 | 0.120 |
| 0.40 | 0.1586 | 0.1692 | 0.155 |
| 1.00 | 0.3011 | 0.3925 | 0.279 |

The full covariance-weighted dual (`alpha=1`) also failed closed-loop: QP success 17.9%, 985 fallbacks, and projector asymmetry 7.98.

## Closed-loop result for alpha=0.02

The smallest genuinely oblique setting remained safe and feasible:

- dual error: `1.04e-14`
- tracked oblique coverage error: `1.03e-14`
- QP success: 100%
- fallbacks: 0
- upper/absolute violations: 0

Compared with `copyV` (`alpha=0`):

- y1 MAE worsened 0.98%; y2 MAE improved 0.27%
- y1 RMSE worsened 1.04%; y2 RMSE improved 0.20%
- average cost worsened 1.16%
- reconstruction residual worsened 0.09%
- input RMS increased slightly

## Conclusion

For the current equal-variance independent sensor-noise experiment, a general oblique extractor is not automatically better. Full obliqueness strongly amplifies errors and breaks chance-constrained feasibility; mild obliqueness is usable but provides no overall performance gain over `R=P`.

This does **not** prove oblique projection is useless. A separate heteroscedastic/correlated sensor-noise experiment is required to test the setting where read and write directions should genuinely differ. The strict free-subspace Mo--Qin `copyR` is a different question: its QP fails at step 1 because it does not guarantee tracked-output coverage, despite dual-basis errors near `1e-14`.
