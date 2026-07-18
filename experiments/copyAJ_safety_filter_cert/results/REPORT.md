# copyAJ safety-filter certificate (Opinion 8 L2-prime)

Guarantee-layer upgrade after primary alpha fail: **change the certified object**.

Stress: `y_max=0.55`, `T_cl=500`, `alpha_joint=0.10`, seeds=[95001 95002 95003].

## Certificate levels (naming hard rule)

| level | meaning |
|---|---|
| `qp_primary` | primary N-step chance QP success (original alpha by construction) |
| `soft_relaxed` | old soft recovery (risk inflate / short N / bound-only) |
| `mean_safety_filter` | **L2-prime** one-step deterministic mean hard-bound filter |
| `uncertified_fallback` | sat(u_mean), no certificate |

**Never** labeled as `qp_original` / `original_alpha` for recovery.

## Aggregate by mode

| mode | n | mean cover | mean MAE1 | mean pf | mean uncert | mean msf_rate | mean soft_rate | mean msf_given_fail | mean msf_also_orig |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| none | 3 | 0.9960 | 0.6006 | 262.0 | 0.5240 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| soft_relaxed | 3 | 0.9973 | 0.4129 | 262.0 | 0.0000 | 0.0000 | 0.5240 | 0.0000 | 0.0000 |
| mean_safety_filter | 3 | 0.7473 | 0.3088 | 262.0 | 0.0000 | 0.5240 | 0.0000 | 1.0000 | 0.0000 |

## mean_safety_filter stage mix (per seed)

| seed | pf | n_msf | n_unc | qp | prev | umean | cover | MAE1 | msf_given_fail | msf_also_orig |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 95001 | 274 | 274 | 0 | 274 | 0 | 0 | 0.7360 | 0.3123 | 1.0000 | 0.0000 |
| 95002 | 258 | 258 | 0 | 258 | 0 | 0 | 0.7480 | 0.3090 | 1.0000 | 0.0000 |
| 95003 | 254 | 254 | 0 | 254 | 0 | 0 | 0.7580 | 0.3051 | 1.0000 | 0.0000 |

## L2-prime object definition

After primary fail, solve for first-step `uk` in input bounds s.t.

```text
H * y_mean_pred(uk)  <=  h     (and lower-side mirror)
y_mean_pred(uk) = y_bar + P A z + P B (uk - u_bar)
```

Optional backup: accept `u_prev` / `u_mean` only if they pass the same mean rows.

## Honest non-claims

1. **Not recursive feasibility.**
2. **Not original-alpha chance-constraint recovery.** `msf_also_orig` is diagnostic only.
3. **Not** a multi-step joint chance guarantee.
4. Mean hard bound is a **different certified object** (deterministic one-step mean).
5. Soft remains an engineering relaxed certificate; MSF is L2-prime mean-safety.

## Empirical takeaways (filled after real MATLAB run)

| mode | mean cover | mean MAE1 | mean pf | mean uncert | new-level rate | msf_given_fail | msf_also_orig |
|---|---:|---:|---:|---:|---:|---:|---:|
| none | 0.9960 | 0.6006 | 262.0 | 0.5240 | — | 0 | — |
| soft_relaxed | **0.9973** | 0.4129 | 262.0 | **0** | soft 0.524 | 0 | — |
| mean_safety_filter | 0.7473 | **0.3088** | 262.0 | **0** | **msf 0.524** | **1.000** | **0** |

1. **L2-prime recovery rate:** `msf_success_given_fail = 1.0` on 3/3 seeds (all primary fails got a mean-feasible `uk` via `one_step_mean_qp`).
2. **Not original-alpha:** `msf_also_orig = 0` everywhere — safety filter never restores original-alpha constant-hold certificate (consistent with AG/AH).
3. **Cover tradeoff:** MSF empirical joint cover falls to ~0.75 because the certified object dropped chance tightening; none/soft stay ~0.996.
4. **Tracking:** MSF has best MAE1 (~0.31) under this stress; soft is middle (~0.41); none worst (~0.60).
5. **Soft still best engineering cover** with zero uncertified steps; MSF gives a *named different* certificate with full primary-fail coverage but weaker stochastic cover.
6. **Relation to AG/AH:** those proved same-alpha backup ladders yield `succ|fail=0`. copyAJ does not fight that; it certifies a weaker object on purpose.

## Metrics definition

- `cover` = 1 - fraction of steps with any tracked y > y_max
- `MAE` = warm-start MAE on tracked channels
- `uncert_rate` = fraction of steps labeled `uncertified_fallback`
- `msf_rate` = fraction labeled `mean_safety_filter`
- `msf_given_fail` = n_msf / primary_fail  (safety-filter success after primary fail)
- `msf_also_orig` = diagnostic only (constant-hold original-alpha recheck)

## Artifact paths

- `results/copyAJ_safety_filter_cert_metrics.csv`
- `results/copyAJ_safety_filter_cert_data.mat`
- `patches/opinion8_safety.md`
