# copyAO_crte_teacher_profiled_unknown_noise

Complete candidate-profiled CRTE teacher evaluation under unknown true sensor noise.

## Evaluated objective

For every complete free subspace `X in St(d, ell-q)` and metric `mu`:

1. build metric-dual `(P,R)`;
2. extract candidate DLVs and refit VARX;
3. evaluate multi-step free-DLV prediction residual;
4. evaluate exact OLS-FWL task trace on `V'*B_T*V > 0`;
5. evaluate two-fold cross-fitted free-output residual-noise proxy;
6. evaluate candidate-specific `Ru^{-1}` finite-horizon authority;
7. reject invalid/rank-deficient/unstable candidates;
8. select minimum profiled teacher objective.

No true `Sigma_n` enters the function signature or scoring path. This is an
unknown-noise proxy variant, not the known-`Sigma_n` formula.

## Search claim

Finite deterministic candidate verification only:

- five `mu` values;
- nine spectral initializers per `mu`;
- thirty seeded random Stiefel candidates per `mu` in the full run.

No global Stiefel/Grassmann optimum is claimed.

## Run

```matlab
cd('.../copyAO_crte_teacher_profiled_unknown_noise')
test_crte_profiled_teacher_unknown_noise
copyAO_crte_teacher_profiled_unknown_noise
```

## Last focused test

```text
PASS profiled teacher unknown-noise: candidates=39 feasible=39 selected=12
J=0.117243 [pred=0.0746292 task=0.000193774 noise=0.0854213]
reach=2.900e-02 dual=2.54e-15 rho=0.9201 val=0.1081
```

## Last full run

```text
selected mu=0.500 source=spectral
teacher J=0.0332907 [pred=0.0237261 task=0.00105902 noise=0.0106236]
reach=0.3943 val NRMSE=0.1205
195/195 candidates feasible; uses true Sigma_n=0
dual=1.38e-15; rho(A)=0.9356
MAE=[0.0595 0.0611]; RMSE=[0.0805 0.0829]
QP success=1; fallback=0; actual upper violations=[1 0]
```

## Important result boundaries

- Teacher optimum validation rank was 88/195; validation optimum teacher rank
  was 135/195. Teacher and prediction validation are not equivalent.
- Task term was 1-3 orders smaller than prediction/noise over the pool; term
  scaling or preregistered alpha/beta calibration is still required.
- Unknown-noise residual proxy is not true sensor-noise covariance.
- One observed plant violation remains despite all hard QPs succeeding.
- No recursive-feasibility, stability, calibration, or global-optimum claim.
