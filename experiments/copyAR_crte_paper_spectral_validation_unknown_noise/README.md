# copyAR_crte_paper_spectral_validation_unknown_noise

## What this copy is

Paper draft **experimental algorithm** (Sec. 3.3–3.4) with **unknown-noise residual proxy**:

1. Free-complement metric \(C_\mu\).
2. Fixed free-source \(S_{yu}\), \(A_T\).
3. Noise term = **one-step free-output residual covariance** (never true \(\Sigma_n\)).
4. Paper \(N_{tr}(A)=A/\max\{|\mathrm{tr}(A)|/d,\varepsilon\}\), with \(\varepsilon=10e-6\).
5. \(A_{CRTE}=N_{tr}(S)+N_{tr}(A_T)-N_{tr}(C_n)\) with weights 1.
6. Dual basis + VARX per \(\mu\); task/noise/reach gates.
7. **Validation-only** selection; rebuild on all offline data; closed-loop SMPC.

The final runner searches the interior grid \(\mu\in\{0.10,0.25,0.50,0.75\}\), excluding the two endpoint metrics. For plant seed `20260710`, the best valid interior candidate is \(\mu=0.10\) with validation NRMSE approximately `0.064557` before the all-data refit.

Not profiled teacher. Not min-teacher copyAO. Not known-\(\Sigma_n\) Oracle.

## Diff

| Copy | Surrogate | Noise | Selection |
|---|---|---|---|
| copyAN | fixed spectral + Frobenius Ntr | residual proxy | validation |
| **copyAR** | fixed spectral + **paper Ntr** | residual proxy | validation |
| copyAO | profiled teacher | residual proxy | min teacher |

## Run

```matlab
cd('.../experiments/copyAR_crte_paper_spectral_validation_unknown_noise')
test_crte_paper_spectral_varx
copyAR_crte_paper_spectral_validation_unknown_noise
```

Plant seed `20260710`.
