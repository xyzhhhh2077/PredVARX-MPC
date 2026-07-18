# copyAH_multistep_alpha_backup

## Purpose

Follow-up to copyAG: after primary original-alpha QP fails, try **multi-step**
backup plans that still satisfy **original** `alpha_joint` chance rows on the
**full planned U**, not only constant-hold of first-step `uk`.

Ladder:

1. short-horizon same-alpha QP → embed/hold to full N → recheck full U
2. one-step same-alpha QP → embed hold → recheck
3. feasibility QP at original alpha (`min ||U-U0||^2` s.t. rows)
4. reduced-Q tracking at original alpha
5. certified constant-hold `u_prev` / `u_mean`

Compare vs `none` and `soft_relaxed`.

## Non-claims

- Not recursive feasibility / stability
- soft_relaxed still not original-alpha unless measured otherwise

## Run

```matlab
cd('.../experiments/copyAH_multistep_alpha_backup');
addpath(pwd); addpath('lib');
run('copyAH_multistep_alpha_backup.m');
```
