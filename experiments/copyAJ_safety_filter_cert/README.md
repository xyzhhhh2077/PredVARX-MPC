# copyAJ_safety_filter_cert

Opinion 8 guarantee-layer upgrade (L2'): **change the certified object** after primary chance-QP fail.

## Design

- Primary path unchanged: N-step Boole chance QP at original `alpha_joint`.
- On primary fail, do **not** re-stack original-alpha heuristics (copyAG/AH already show `succ|fail=0`).
- Instead run a **one-step deterministic mean hard-bound safety filter**:
  - certificate level: `mean_safety_filter` (never `qp_original` / `original_alpha`)
  - object: `H * y_mean_pred(uk) <= h` (+ lower mirror), `uk` in input bounds
  - ladder: mean-constrained one-step QP → screened `u_prev` → screened `u_mean`
  - else: `uncertified_fallback` = sat(`u_mean`)

## Modes compared (tight `y_max=0.55`, ≥3 seeds)

| mode | recovery after primary fail |
|---|---|
| `none` | uncertified sat(u_mean) |
| `soft_relaxed` | old risk-inflate / short-N / bound-only |
| `mean_safety_filter` | L2' mean hard-bound filter |

## Non-claims

- Not recursive feasibility / stability
- Not original-α chance recovery success
- Soft success ≠ original-α; MSF success ≠ original-α

## Entry

```matlab
cd(.../experiments/copyAJ_safety_filter_cert)
copyAJ_safety_filter_cert
```

Outputs: `results/REPORT.md`, `results/copyAJ_safety_filter_cert_metrics.csv`, `patches/opinion8_safety.md`.
