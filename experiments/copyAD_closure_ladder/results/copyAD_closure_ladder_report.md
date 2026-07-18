# copyAD closure ladder report

Sequential empirical trials for remaining open opinions. **Not theorem closures.**

Shared plant/seed family with copyAB/AC offline data `rng(20260710)`; each config has its own `cl_seed`.
Horizon here `T_cl=800` (shorter ladder runs).

## P0 Soft recovery under tight y_max=0.55

| config | primary_fail | soft_ok | uncert | soft_success|fail | joint_cover | MAE |
|---|---:|---:|---:|---:|---:|---:|
| P0_stress_no_soft | 476 | 0 | 476 | 0.000 | 0.9625 | 0.5383/0.4981 |
| P0_stress_with_soft | 478 | 478 | 0 | 1.000 | 0.9975 | 0.4110/0.3234 |
| P0_stress_soft_plus_cross | 749 | 749 | 0 | 1.000 | 1.0000 | 0.5949/0.4259 |

Non-claim: soft success is a **different** certificate (risk inflate / short N / bounds-only), not original alpha.

## P1 Sigma_eps denominators (coverage proxy)

| config | joint_cover | target | MAE1 | MAE2 | active |
|---|---:|---:|---:|---:|---:|
| P1_eps_t2 | 1.0000 | 0.90 | 0.0593 | 0.0576 | 0.146 |
| P1_eps_ml | 1.0000 | 0.90 | 0.0593 | 0.0576 | 0.134 |
| P1_eps_ols | 1.0000 | 0.90 | 0.0593 | 0.0576 | 0.240 |

Non-claim: joint empirical cover is only a proxy; not a proof of correct residual DOF.

## P1 Sigma_obs modes

| config | joint_cover | MAE1 | MAE2 | active |
|---|---:|---:|---:|---:|
| P1_obs_declared_shape | 1.0000 | 0.0554 | 0.0554 | 0.191 |
| P1_obs_residual_support | 1.0000 | 0.0554 | 0.0554 | 0.182 |
| P1_obs_additive | 1.0000 | 0.0554 | 0.0554 | 0.186 |

Non-claim: modes are engineering objects; none identified with Cov(o) theorem.

## P2 Cross-cov and residualize

| config | joint_cover | MAE1 | MAE2 | qp | active |
|---|---:|---:|---:|---:|---:|
| P2_cross_off | 1.0000 | 0.0610 | 0.0570 | 1.000 | 0.168 |
| P2_cross_on | 1.0000 | 0.0610 | 0.0570 | 1.000 | 0.232 |
| P2_resid_off | 1.0000 | 0.0599 | 0.0585 | 1.000 | 0.218 |
| P2_resid_on | 1.0000 | 0.0599 | 0.0585 | 1.000 | 0.218 |
| REF_AB_like | 1.0000 | 0.0605 | 0.0595 | 1.000 | 0.168 |
| REF_AC_like | 1.0000 | 0.0605 | 0.0595 | 1.000 | 0.337 |

## Still open (theorem level)

1. Recursive feasibility / stability (opinion 9) — not attempted beyond soft terminal cost elsewhere.
2. Soft recovery original-risk certificate.
3. Statistically optimal Sigma_eps denom; Sigma_obs = Cov(o); Boole with cross terms; input-conditional PredVARX optimum.

## Empirical takeaways (this ladder)

### P0 Soft recovery — **practical path works under stress**

| | no soft | with soft | soft+cross |
|---|---:|---:|---:|
| primary QP fail | 476 | 478 | 749 |
| soft recoveries | 0 | **478** | **749** |
| uncertified | **476** | **0** | **0** |
| soft success \| fail | 0 | **1.000** | **1.000** |
| joint cover | 0.9625 | **0.9975** | **1.0000** |
| MAE | 0.54/0.50 | **0.41/0.32** | 0.59/0.43 |

- Enabling soft recovery **eliminated uncertified fallback** on this tight-\(y_{\max}\) seed and improved cover/MAE vs hard uncertified hold.
- Soft+cross made primary QP fail more often (tighter \(\Sigma_y\)) but soft still absorbed **all** failures.
- **Still not** the original \(\alpha\) certificate.

### P1 denoms / obs modes — **coverage saturated; active rate differs**

- At \(y_{\max}=2\), joint empirical cover = 1.000 for all t2/ml/ols and all three \(\Sigma_{obs}\) modes → **cannot rank statistical correctness** by violation rate here.
- Useful signal is **chance-row activity**, not cover:
  - eps: t2 0.146, ml 0.134, **ols 0.240** (ols tightens more on this seed)
  - obs modes: activity ~0.18–0.19, MAE identical to 1e-6
- Need rarer-event / multi-seed stress to calibrate cover; this ladder only shows **local indifference of MAE** under loose limits.

### P2 cross / residualize — **small local effect**

- Cross ON vs OFF: MAE ~identical; active 0.168 → **0.232** (cross changes tightening, little tracking change).
- Residualize ON vs OFF: MAE differs only at ~1e-5; activity same on this seed.
- Longer \(T_{cl}=1200\) AB/AC fair compare previously showed larger stack effect; single-factor flips here are weak.

### Opinion 9

Not closed. No terminal-set / Lyapunov recursive feasibility experiment in copyAD.
