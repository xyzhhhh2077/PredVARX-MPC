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

## Strict Sec. 5.3 FWL support

Before candidate generation, the implementation computes
`Z0=W*M0=Us*Ss*Vs'`, fixes a numerical rank, and generates every candidate
inside `range(B_T)` as `V0=Us*Ss^{-1}*Zeta`. It then applies paired metric
normalization `V=V0*(V0'*C_mu*V0)^{-1/2}`. Rank smaller than `ell_f` raises
`InsufficientFWLRank`; denominator ridge is not used.

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
PASS profiled teacher unknown-noise + Sec5.3 SVD support:
rank=8 support-I=3.58e-15 support-res=5.72e-15
candidates=39 feasible=39 selected=12
J=0.117243 [pred=0.0746292 task=0.000193774 noise=0.0854213]
reach=2.900e-02 dual=1.15e-15 rho=0.9201 val=0.1081
```

## Last full run

```text
selected mu=1.000 source=random
teacher J=0.0297283 [pred=0.0174975 task=0.000147351 noise=0.0123782]
reach=0.0100 val NRMSE=0.1419
195/195 candidates feasible; uses true Sigma_n=0
support rank=28; support-I=8.69e-14; max support residual=2.70e-14
dual=4.50e-16; rho(A)=0.9289
MAE=[0.0630 0.0626]; RMSE=[0.0847 0.0845]
QP success=1; fallback=0; actual upper violations=[1 0]
```

## Important result boundaries

- Teacher optimum validation rank was 169/195; validation optimum teacher rank
  was 164/195. Teacher and prediction validation are not equivalent.
- Task term was 1-3 orders smaller than prediction/noise over the pool; term
  scaling or preregistered alpha/beta calibration is still required.
- Unknown-noise residual proxy is not true sensor-noise covariance.
- One observed plant violation remains despite all hard QPs succeeding.
- No recursive-feasibility, stability, calibration, or global-optimum claim.
