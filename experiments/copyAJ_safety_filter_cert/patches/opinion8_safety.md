# Opinion 8 — safety filter / backup certificate (copyAJ)

**Status:** IMPLEMENTED + MATLAB verified  
**Copy:** `experiments/copyAJ_safety_filter_cert/` only  
**Date:** 2026-07-19

## Goal

After primary N-step original-α chance QP fails, **do not** stack more same-α online heuristics
(copyAG constant-hold / copyAH multi-step both had `succ|fail=0`).

Instead **change the certified object** to a one-step deterministic mean hard-bound
safety filter and label it honestly as L2-prime:

| level | meaning |
|---|---|
| `qp_primary` | primary chance QP OK |
| `soft_relaxed` | old soft recovery |
| `mean_safety_filter` | L2' mean hard-bound filter |
| `uncertified_fallback` | sat(u_mean) |

Hard naming rule: recovery is **never** called `qp_original` / `original_alpha`.

## Implementation

- `lib/mean_safety_filter.m` — one-step QP:
  - `H * y_mean_pred(uk) <= h` (+ lower mirror), `uk` in input bounds
  - fallback candidates: screened `u_prev`, `u_mean`
- `run_safety_config.m` — modes `none | soft_relaxed | mean_safety_filter`
- `copyAJ_safety_filter_cert.m` — stress `y_max=0.55`, `T_cl=500`, seeds `95001..95003`

## MATLAB evidence (real run)

Stdout marker: `COPYAJ_DONE`

### Aggregate (3 seeds)

| mode | mean cover | mean MAE1 | mean pf | mean uncert | mean new-level rate | msf\|fail | msf_also_orig |
|---|---:|---:|---:|---:|---:|---:|---:|
| none | 0.9960 | 0.6006 | 262.0 | 0.5240 | 0 | 0 | — |
| soft_relaxed | **0.9973** | 0.4129 | 262.0 | **0.0000** | soft 0.524 | 0 | — |
| mean_safety_filter | 0.7473 | **0.3088** | 262.0 | **0.0000** | **msf 0.524** | **1.000** | **0.000** |

### MSF stage mix

All recovered steps used `one_step_mean_qp` (274/258/254); prev/umean backups unused.
`uncert_rate=0` under MSF on all 3 seeds.

### Diagnostics

- `msf_also_orig_rate = 0` on all seeds: MSF inputs **never** pass original-α constant-hold.
- `soft_also_mean_rate = 1` on soft mode: soft recoveries happen to also satisfy mean hard bounds
  (different certificate; does not make soft = original-α).

## Honest takeaways

1. **Object change works as an L2' layer:** after every primary fail, a mean-hard-bound
   `uk` was found (`msf_success_given_fail=1`). This is a real recovery rate improvement
   vs copyAG/AH original-α ladders (`succ|fail=0`).
2. **This is NOT original-α recovery.** `msf_also_orig=0` confirms the certificate is different.
3. **This is NOT recursive feasibility / multi-step chance guarantee.**
4. **Empirical cover drops** (0.75 vs ~0.996 for none/soft) because the mean hard bound
   drops Boole/Gaussian tightening — stochastic overshoot of `y_max` is expected.
5. **Soft remains best engineering cover/MAE tradeoff** under this stress; MSF trades
   cover for a *named* deterministic one-step mean certificate and lower MAE.
6. Soft is still L2 `soft_relaxed`; MSF is L2' `mean_safety_filter` — do not collapse labels.

## Non-claims (mandatory)

- Not recursive feasibility
- Not original-α chance-constraint recovery success
- Not infinite-horizon joint chance guarantee
- Mean-feasible open-loop one-step ≠ closed-loop safety invariant set

## Files

```
experiments/copyAJ_safety_filter_cert/
  copyAJ_safety_filter_cert.m
  run_safety_config.m
  lib/mean_safety_filter.m
  results/REPORT.md
  results/copyAJ_safety_filter_cert_metrics.csv
  results/copyAJ_safety_filter_cert_data.mat
  patches/opinion8_safety.md
```

## PASS/FAIL vs task gates

| gate | result |
|---|---|
| L2' not named qp_original/original_alpha | PASS |
| vs none/soft, y_max=0.55, ≥3 seeds | PASS (3 seeds) |
| metrics: cover, MAE, uncert, new-level rate, msf\|fail | PASS |
| MATLAB real run | PASS (`COPYAJ_DONE`) |
| REPORT.md + patches/opinion8_safety.md | PASS |
| honest non-claims | PASS |
| no commit | PASS |
| only copyAJ edited | PASS |
