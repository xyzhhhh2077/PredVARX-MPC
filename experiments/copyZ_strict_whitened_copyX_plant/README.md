# copyZ_strict_whitened_copyX_plant

## Question

What happens if the copyX plant and SMPC are kept fixed, but the identifier is replaced by a rank-full, first-order Mo--Qin-style normalized-space IVR and direct de-normalization ablation?

> **Correction after source-to-paper audit:** copyZ performs the whitening and direct de-normalization equations, but it is not a complete line-by-line implementation of Algorithm 1. It fixes the VAR order to $s=1$, uses a tolerance-based stopping rule rather than accepting updates only while the predicted trace increases, relies on the full-rank noisy case instead of implementing the paper's rank-null complement, and reports separately fitted VARX covariance quantities rather than all prescribed Algorithm-1 covariance outputs.

## Fair comparison with copyX

The following are kept identical to copyX:

- random seed and synthetic plant;
- dimensions, offline and closed-loop sample counts;
- offline excitation sequence and noise realization;
- reference trajectory and smooth noise envelopes;
- prediction horizon, tracking/input weights, input bounds;
- joint chance-constraint allocation and online covariance-window update.

Only the identification geometry changes.

## Implemented whitening path

For centered output data $Y^c$, compute

$$
\Sigma_y=UDU^T,
\qquad
Y^*=D^{-1/2}U^TY^c.
$$

Iterative variance-ratio estimation is performed in $Y^*$ and returns an orthonormal $P^*$:

$$
P^{*T}P^*=I.
$$

The realization is then de-normalized directly:

$$
P=UD^{1/2}P^*,
\qquad
R=UD^{-1/2}P^*.
$$

The complementary pair uses the same construction:

$$
\bar P=UD^{1/2}\bar P^*,
\qquad
\bar R=UD^{-1/2}\bar P^*.
$$

No QR orthogonalization, post-hoc SVD alignment, generic dual completion, tracked-axis insertion, or $R_\alpha$ interpolation is applied to this rank-full ablation.

## PredVARX boundary

The original paper identifies a VAR model without an explicit control input. After the implemented whitening/de-normalization realization is obtained, this experiment adds the centered VARX extension

$$
z_k=R^T(y_k-\bar y),
\qquad
z_{k+1}=Az_k+B(u_k-\bar u)+\varepsilon_k.
$$

This input channel and the subsequent SMPC are project extensions, not claims about the original paper.

## Control-coverage boundary

The strict predictive subspace is left free. It is not forced to satisfy

$$
PR^TE_c=E_c.
$$

Therefore a loss of QP feasibility or tracking is retained as experimental evidence that predictive-subspace identification alone does not guarantee control coverage.

## Verification

The focused test checks normalized orthogonality, all four full dual-basis identities, direct raw de-normalization, and PSD covariance estimates. The full runner persists whitening, dual-basis, tracked-coverage, tracking, feasibility, chance-residual, noise, and control metrics.

## Real MATLAB result

The implemented whitening/de-normalization identities are algebraically correct:

- normalized orthogonality: `1.88e-15`;
- four full dual-basis errors: `[8.01e-15, 1.88e-15, 3.69e-14, 1.04e-14]`.

It does not cover the two tracked output axes on this control task:

- tracked coverage error: `1.6327`;
- QP success: `39.67%`;
- fallbacks: `724/1200`;
- MAE: `[0.5172, 0.7246]`;
- average full cost: `786.36`.

For the same plant and controller, copyX gives exact tracked coverage, 100% QP success, no fallbacks, MAE `[0.0973, 0.0950]`, and average cost `109.16`.

Thus this rank-full, $s=1$ whitening/de-normalization ablation fails as a direct replacement for the control-aware subspace in this experiment. This supports the narrower conclusion that the implemented free predictive subspace does not guarantee control-output coverage. It must not be cited as a complete reproduction or definitive falsification of Mo--Qin Algorithm 1.
