# copyAX — ControlGym CDR (Convection-Diffusion-Reaction) soft-preference validation

Apply the copyAU anchored-CRTE soft-preference construction to a **linear,
fully-actuated, disturbance-free benchmark** from
[ControlGym](https://github.com/xiangyuan-zhang/controlgym) (Zhang et al.,
L4DC 2024, [arXiv:2311.18736](https://arxiv.org/abs/2311.18736)).

This is the replacement validation object for the BOPTEST HVAC pilot
(`copyAW_boptest_soft_preference`), which failed because the dominant
exogenous weather/occupancy disturbance cannot be modeled in the current
PredVARX form (no `d_k` channel).

## Why CDR fits the acceptance criteria

| Criterion | CDR status |
|---|---|
| `u` actively manipulated | 8 localized actuator supports, PRBS-excited |
| `y` high-dimensional | 30 sensors over the 200-grid field (C has one 1 per row) |
| public data/simulator | ControlGym open-source, seeded, deterministic |
| original paper + standard evaluation | L4DC paper; one-step prediction vs persistence |
| MPC-connectable | first-order linear `x_{k+1}=A x_k+B2 u_k`, LQ benchmark built-in |
| explicit `d_k` channel | absent (process noise disabled) — matches "d_k may be absent" |
| model form | CDR is *exactly* a first-order VARX in the state; the latent model must recover it from 30 sensors |

## Model being verified

Environment (analytical, linear):
```
x_{k+1} = A x_k + B2 u_k + w_k,      w_k = 0 (process noise off)
y_k     = C x_k + v_k,               v_k ~ N(0, 1e-4 I)  (std 0.01)
```
Identified model (copyAU construction, segmented):
```
z_k = R' (y_k - ybar),  z_k in R^5
z_{k+1} = Ahat z_k + Bhat u_k + eps
```
with the task anchor `E` (q=2) learned from a soft spatial preference
(central sensors 1.0 / edge sensors 0.85) blended with input reachability.

**Time alignment (standard causality)**: the ControlGym `step()` returns the
observation computed from the state *before* the update, so `u_k` drives
`y_{k+1}`. All three MATLAB functions use `ur = uc(:, valid)`. The contract
test `tests/cdrSegmentContractTest.m` verifies this by recovering
`y(k)=0.4 y(k-1)+1.7 u(k-1)` from the fit and by rejecting the BOPTEST-style
(mis)alignment.

## Layout

```
collect_cdr_data.py                  seeded PRBS collection -> data/*.mat
copyAX_cdr_soft_preference.m         main runner (learn E -> fit -> evaluate)
fit_segmented_anchored_varx.m        copyAU anchored CRTE fit (aligned u_k)
learn_segmented_preferred_output_directions.m
evaluate_segmented_prediction.m
tests/test_cdr_dataset.py            Python dataset contract (6 tests)
tests/cdrSegmentContractTest.m       MATLAB contract (3 tests)
data/copyAX_cdr_dataset.mat          y 30x8008, u 8x8008, segments 1..8 (LFS)
data/copyAX_cdr_dataset.json          full generation record
results/                             metrics txt + fig png + full state .mat
```

## Data

- 8 segments x 1001 columns; segments 1–6 train (6000 transitions),
  7–8 validation (2000 transitions); contiguous, never shuffled.
- Two-level PRBS, amplitude ±1, held 4 samples, per-channel seeded:
  `seed = 20260801 + segment*101`.
- `sample_time = 0.01` (10 Hz). Default `reaction_constant = 0.1` gives an
  open-loop growing DC mode (`exp(0.001)` per step); 1000-step segments keep
  the field in O(1) range, well inside the unconstrained limits.
- Sensor noise std 0.01 (SNR ~40 dB at O(1) field scale); no process noise.

## Acceptance criteria (this copy)

1. `dual_error < 1e-7` and anchor preservation ~0 (geometry intact).
2. Identified `Ahat` stable (`spectral_radius < 1`).
3. Validation `task_R2 > 0` and `task_RMSE < persistence_RMSE` — the
   criterion BOPTEST failed (`R2 = [-1.136, 0.3257]`).
4. Latent model improves on persistence in standardized units.
5. Learned contribution spatially interpretable (concentrated near the
   central actuated region, consistent with the soft preference).

If CDR fails structurally, the fallback order is ControlGym **Wave** (linear,
oscillatory) then **Kuramoto–Sivashinsky** (nonlinear stress test).

## Verification summary (2026-08-01) — CDR is suitable

Followed the model → parameters → results loop; CDR passed on the first run,
so no fallback object was needed.

**Model adequacy**
- Identified latent spectrum stable: `spectral_radius = 0.9981` (the true DC
  mode grows `exp(+0.001)` per step; the latent model lands just below 1).
- Geometry exact: `dual_error = 2.0e-15`, anchor preservation `3.2e-16`.
- Task-axis residuals are tiny: latent RMSE `[0.00193 0.00277 ...]` vs
  persistence `[0.00207 0.00297 ...]` — the two controlled axes are the best
  modeled directions.
- Residual autocorrelation is high (lag-1 `max|acf| ≈ 0.91`) on the free-block
  dims — expected for a 5-dim truncation of a 200-dim field, not a failure
  criterion; the task axes (what the MPC chance constraints use) are clean.

**Parameter adequacy**
- Regression Gram condition number `6.55e4` (BOPTEST was `6.5e10`; PRBS
  excitation makes the identification well-posed).
- Sample/parameter ratio `92.3` (BOPTEST `1.35`), 6000 training transitions.
- `S_yu` eigenvalue decay `[14.5 6.77 5.32 1.15 0.81 0.61 ...]` supports
  `ell = 5`; free-block CRTE eigenvalues `[25.2 2.2 1.0]` show one dominant
  residual direction.

**Result validity**
- Both held-out segments beat persistence: task `R2 = [0.80 0.80]` (segment 7,
  n=1000) and `[0.82 0.78]` (segment 8, n=1000); pooled `task_R2 = [0.814 0.802]`,
  `task_RMSE < task_persistence_RMSE` on both axes (BOPTEST: `R2 = [-1.136 0.326]`).
- Full 30-output standardized one-step RMSE `0.491`, every latent dim improves
  on its persistence baseline.
- Contribution is spatially interpretable: top sensors 16/17/15/20 (grids
  90–120, the central actuated region around control supports 3–4), minimum
  contribution 1e-4 at domain edges; consistent with the central soft
  preference (1.0) vs edges (0.85).

Verdict: **suitable**. The linear, disturbance-free CDR benchmark isolates the
CRTE construction (learned output directions + anchored free block) from the
confound that sank BOPTEST (dominant exogenous `d_k`), and the construction
clearly wins on held-out segments.

## Scope

Offline, one-step prediction validation on held-out segments only. No MPC
loop is closed here; the linear CDR benchmark is MPC-ready for a follow-up.
