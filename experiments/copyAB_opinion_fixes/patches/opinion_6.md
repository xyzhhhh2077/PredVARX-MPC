# Opinion 6 patch report

status: done

## Gap

$R^\top o = 0$ samplewise (when $R^\top P = I$) does **not** imply $\mathrm{Cov}(z,o)=0$.
The SMPC construction in `centered_smpc_step`

```matlab
Sigma_y = P*Sigma_z*P' + model.Sigma_obs
```

silently drops the cross blocks $P\Sigma_{zo} + \Sigma_{zo}^\top P'$.

## Algebra

With centered samples $y_c$, $z=R^\top y_c$, $o=(I-PR^\top)y_c$:

$$
y_c = Pz + o
$$

$$
\mathrm{Cov}(y) = P\Sigma_z P^\top + P\Sigma_{zo} + \Sigma_{zo}^\top P^\top + \Sigma_o
$$

Full block reconstruction relative error is $\sim 0$ by construction.
Drop-cross relative error can be strictly positive for generic / oblique data.

## Owned artifacts

| File | Role |
|---|---|
| `lib/cross_cov_diagnostics.m` | diagnostic: `y,P,R → z,o,Σz,Σo,Σzo,drop_cross_rel_err` |
| `tests/test_opinion06_cross_cov.m` | random-data gates: full~0, drop>0, split oblique>0 |
| `copyAB_opinion_fixes.m` | offline one-shot `fprintf` after identification (main loop untouched) |

## What this does **not** fix

- Does not change the online SMPC $\Sigma_y$ formula (still drop-cross).
- Does not estimate or inject $\Sigma_{zo}$ into chance constraints.
- Diagnostic only: documents the missing cross term magnitude on offline data.

## Verification

```matlab
cd('.../experiments/copyAB_opinion_fixes');
addpath(pwd); addpath(fullfile(pwd,'lib')); addpath(fullfile(pwd,'tests'));
test_opinion06_cross_cov
```

Expected: `PASS opinion06 cross-cov: full_rel≈0, drop_generic>0, ...`
