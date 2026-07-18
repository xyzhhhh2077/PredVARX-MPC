# copyAG_original_alpha_cert

## Purpose

Opinion **8 guarantee layer** experiment after soft recovery was shown to work as a *relaxed* certificate.

Compare three recovery modes under tight `y_max=0.55`:

| mode | meaning | original-alpha? |
|---|---|---|
| `none` | `sat(u_mean)` on primary fail | no |
| `soft_relaxed` | old risk-inflate / short-N / bound-only | **no** |
| `original_alpha` | only accept inputs that re-pass original `alpha_joint` rows | **yes (model-based)** |

## Certificate definition

- **Primary success**: the full planned `U` from the N-step original-alpha QP satisfies `max(A_ch*U-b_ch)<=tol` (certified by construction).
- **Recovery candidate `uk`**: only accepted if the constant-hold extension of `uk` over the original N passes original-alpha rows.

`soft_relaxed` may restore feasibility under a *different* risk/horizon; we also measure how often those recovered `uk` still pass original-alpha constant-hold (`soft_also_orig_rate`).
## Non-claims

- Not recursive feasibility
- Not closed-loop stability
- Not infinite-horizon joint chance constraint
- `soft_relaxed` success ≠ original-alpha certificate

## Run

```matlab
cd('.../experiments/copyAG_original_alpha_cert');
addpath(pwd); addpath('lib');
run('copyAG_original_alpha_cert.m');
```

## Outputs

- `results/copyAG_original_alpha_cert_metrics.csv`
- `results/copyAG_original_alpha_cert_report.md`
