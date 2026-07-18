# copyAH multi-step original-alpha backup

Opinion 8 guarantee layer v2. **Not RF/stability.**

Stress: `y_max=0.55`, `T_cl=500`, `alpha_joint=0.10`, 5 seeds.

## Aggregate

| mode | mean cover | mean MAE1 | mean pf | mean orig_cert | mean uncert | mean succ\|fail | mean soft_also_orig |
|---|---:|---:|---:|---:|---:|---:|---:|
| none | 0.9848 | 0.6034 | 261.8 | 0.4764 | 0.5236 | 0.0000 | 0.0000 |
| soft_relaxed | **0.9976** | **0.4170** | 261.4 | 0.4772 | **0.0000** | 0.0000 | **0.0000** |
| original_alpha_multistep | 0.9848 | 0.6034 | 261.8 | 0.4764 | 0.5236 | **0.0000** | 0.0000 |

## original_alpha_multistep stage mix

| seed | pf | orig_rec | unc | short | one | feas | redQ | prev | umean | cover | succ\|fail |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 94001 | 253 | 0 | 253 | 0 | 0 | 0 | 0 | 0 | 0 | 0.9940 | 0.0000 |
| 94002 | 274 | 0 | 274 | 0 | 0 | 0 | 0 | 0 | 0 | 0.9860 | 0.0000 |
| 94003 | 258 | 0 | 258 | 0 | 0 | 0 | 0 | 0 | 0 | 0.9760 | 0.0000 |
| 94004 | 271 | 0 | 271 | 0 | 0 | 0 | 0 | 0 | 0 | 0.9960 | 0.0000 |
| 94005 | 253 | 0 | 253 | 0 | 0 | 0 | 0 | 0 | 0 | 0.9720 | 0.0000 |

## Ladder attempted (all failed after primary infeasibility)

1. short-horizon same-α QP → embed to full N → recheck full U  
2. one-step same-α → embed hold → recheck  
3. feasibility QP `min ||U-U0||^2` s.t. original rows + bounds  
4. reduced-Q tracking at original α  
5. certified constant-hold `u_prev` / `u_mean`

## Honest takeaways

1. **Primary original-α certification still works** on successful QP steps (~48% of horizon under this stress).
2. **soft_relaxed remains the only engineering recovery** that clears uncertified steps and improves MAE/cover.
3. **`soft_also_orig_rate = 0` again**: relaxed recovery never yields an original-α constant-hold certificate.
4. **Multi-step original-α backup also fails completely** (`succ|fail=0` on 5/5 seeds), including pure feasibility QP at original α.
5. Interpretation: under this tight `y_max`, once the primary N-step original-α QP is infeasible, the **original-α chance-row set is effectively empty (or numerically unreachable)** for the tested backup generators. This is stronger evidence than copyAG that the guarantee gap is not just “constant-hold too strict”.
6. Closing the gap likely requires changing the certified object itself:
   - certified **smaller risk allocation only with explicit relabeling** (not claiming original α), or
   - a **true safety filter / invariant set** designed offline for original α, or
   - softening constraints inside a *declared* degraded certificate class.

## Non-claims

- Not recursive feasibility
- Not closed-loop stability
- Not a proof that original-α recovery is impossible in all regimes—only that this ladder fails on this stress suite
