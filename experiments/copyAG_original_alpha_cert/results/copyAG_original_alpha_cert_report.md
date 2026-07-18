# copyAG original-alpha certificate experiment

Opinion 8 guarantee layer. **Not recursive feasibility / stability.**

Tight stress: `y_max=0.55`, `T_cl=500`, `alpha_joint=0.10`, 5 seeds.

Recovery modes:

1. `none` — uncertified sat(u_mean)
2. `soft_relaxed` — old risk-inflate/short-N/bound-only (NOT original alpha)
3. `original_alpha` — only accept recovery inputs that re-pass original-alpha constant-hold rows

## Aggregate by mode

| mode | n | mean cover | mean MAE1 | mean pf | mean orig_cert_rate | mean uncert_rate | mean orig_succ|fail | mean soft_also_orig |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| none | 5 | 0.9636 | 0.5586 | 240.2 | 0.5196 | 0.4804 | 0.0000 | 0.0000 |
| soft_relaxed | 5 | 0.9960 | 0.4166 | 240.2 | 0.5196 | 0.0000 | 0.0000 | **0.0000** |
| original_alpha | 5 | 0.9636 | 0.5586 | 240.2 | 0.5196 | 0.4804 | **0.0000** | 0.0000 |

## Per-seed original_alpha stage mix

| seed | pf | orig_rec | unc | shortN | one-step | prev | umean | cover | orig_cert_rate |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 93001 | 219 | 0 | 219 | 0 | 0 | 0 | 0 | 0.9820 | 0.5620 |
| 93002 | 246 | 0 | 246 | 0 | 0 | 0 | 0 | 0.9860 | 0.5080 |
| 93003 | 256 | 0 | 256 | 0 | 0 | 0 | 0 | 0.9160 | 0.4880 |
| 93004 | 257 | 0 | 257 | 0 | 0 | 0 | 0 | 0.9460 | 0.4860 |
| 93005 | 223 | 0 | 223 | 0 | 0 | 0 | 0 | 0.9880 | 0.5540 |

## Certificate definition

- **Primary success**: full planned `U` from N-step original-alpha QP has `max(A_ch*U-b_ch)<=tol` (certified by construction).
- **Recovery `uk`**: accepted only if constant-hold extension of `uk` over original N passes original-alpha rows.

## Empirical takeaways (honest)

1. **Primary original-alpha certification works** on successful QP steps (~52% of horizon under this stress).
2. **soft_relaxed restores engineering feasibility** (uncert→0, cover→0.996, better MAE) but **`soft_also_orig_rate = 0`**: none of the relaxed recoveries pass original-alpha constant-hold. Soft remains a **different certificate**.
3. **original_alpha recovery ladder failed completely** on all 5 seeds after primary fail (`orig_success_given_fail=0`). Short-horizon / one-step / u_prev / u_mean never produced an original-alpha-certified single-step hold under this stress.
4. Therefore, under tight `y_max=0.55`, the guarantee layer is **incomplete**: we can label primary successes, but **cannot currently recover under the original α after primary infeasibility**.
5. Practical implication: either keep soft as an explicitly *relaxed* certificate, or design a stronger original-α backup (terminal set / safety filter / multi-step certified backup plan). Constant-hold recheck alone is too strict here.

## Non-claims

1. Not recursive feasibility.
2. Not infinite-horizon joint chance guarantee.
3. soft_relaxed success does **not** imply original-alpha certification (empirically confirmed: 0%).
4. original_alpha recovery incompleteness is a real remaining open item, not a documentation gap only.
