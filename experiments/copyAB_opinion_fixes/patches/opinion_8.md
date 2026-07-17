# Opinion 8 patch report

status: done

## Claim / design

**Uncertified fallback does not inherit the QP chance-constraint certificate.**

When `centered_smpc_step` throws (infeasible / exitflag ≤ 0), the runner applies

$$
u_k = \mathrm{sat}(\bar u;\, u_{\min}, u_{\max}).
$$

This guarantees **input bounds only**. It does **not** re-solve $A_{ch}U\le b_{ch}$ and does **not** recompute $\mu+\kappa\sigma$ Boole/Gaussian tightening. Therefore the step is labeled

| `cc_cert_level` | Meaning |
|---|---|
| `'qp'` | Successful chance-constrained QP; residual `max(A_ch U - b_ch)` is meaningful |
| `'uncertified_fallback'` | Heuristic sat(u_mean); **no chance-constraint guarantee** |

A one-step **deterministic mean** prediction is computed as a simplified diagnostic hard-bound residual against `opt.h` (`y_max`). Even if the mean lies under `y_max`, the step remains `uncertified_fallback` — mean safety ≠ chance-constraint certification.

## Files changed

| File | Role |
|---|---|
| `lib/fallback_certify_step.m` | **new** helper: mean check + always-uncertified cert grade |
| `copyAB_opinion_fixes.m` | catch path, metrics split, save/print |
| `tests/test_opinion08_fallback_flag.m` | **new** logic unit test (no real QP) |

## Catch-path behavior (runner)

1. `uk = sat(model.u_mean)`
2. `fb = fallback_certify_step(yk, model, opt, uk)`
3. Record:
   - `cc_cert_level{k} = 'uncertified_fallback'`
   - `exitflag(k) = -1`
   - `max_cc_violation(k) = NaN` (no QP residual; not treated as feasible)
   - `fallback_det_mean_violation(k) = max(H * mean_pred - h)` diagnostic
4. First 3 fallbacks print: `QP uncertified fallback at k=... | ... NOT chance-constraint certified`
5. Success path: `cc_cert_level{k} = 'qp'`

## Metrics

- `qp_certified_count` / `qp_certified_rate` — only `'qp'` steps
- `uncertified_fallback_count` / `uncertified_fallback_rate` — only fallbacks
- `max_recorded_qp_constraint` and `constraint_active_rate` use **QP steps only**
- Fallbacks are **not** counted as chance-constraint-guaranteed successes

## What we deliberately do **not** claim

- Fallback is **not** `backup-certified` (no robust invariant set / safety filter proof).
- A safe one-step mean does **not** restore Boole chance guarantees.
- Fallback rate is **not** the same as chance-constraint violation rate.

## Verification

```matlab
cd('.../experiments/copyAB_opinion_fixes');
addpath(pwd); addpath(fullfile(pwd,'lib')); addpath(fullfile(pwd,'tests'));
test_opinion08_fallback_flag
```

Expected: `PASS opinion08: ...`
