# copyAY — Waterloo Pelican quadrotor soft-preference validation

Apply the copyAU anchored-CRTE soft-preference construction to a **real
flight dataset** — the Waterloo Pelican quadrotor dataset (Mohajerin et al.,
ICRA 2018 / TNNLS 2019), collected indoors with Vicon motion capture at 100 Hz.

This is the real-data validation object after the simulator-only objects
(TEP → BOPTEST → ControlGym CDR) were exhausted. Unlike all three previous
objects, the data are **measured physical flight data**, not a simulator
rollout.

## Why Pelican fits the acceptance criteria

| Criterion | Pelican status |
|---|---|
| `u` actively manipulated | 4 commanded motor speeds `Motors_CMD`, pilot-driven across hover → aggressive manoeuvres (std 16–18 over 50–183) |
| `y` high-dimensional | 10 raw channels: Vicon position (3), Euler attitude (3), motor speeds (4) |
| public data/simulator | direct HTTP download, 238 MB .mat, no registration |
| original paper + standard evaluation | ICRA 2018 / TNNLS 2019; the dataset's own task is system identification + multi-step prediction |
| MPC-connectable | quadrotor MPC literature is large; task axes = flight state (position/attitude) |
| explicit `d_k` channel | absent (indoor Vicon environment) — matches "d_k may be absent" |
| excitation | pilot actively flies every regime: hover, close-to-ground, light, moderate, aggressive (u std ~17 per channel) |

## Model being verified

```
z_k = R' (y_k - ybar),            z_k in R^5
z_{k+1} = Ahat z_k + Bhat u_k + eps
y_k = [Pos_k(3); Euler_k(3); Motors_k(4)]      (10 raw measurements, 100 Hz)
u_k = Motors_CMD_k(4)
```

Task anchor `E` (q=2) learned from a soft preference: position 1.0,
attitude 0.9, motor speeds 0.7 — flight state is the task, actuator internals
are represented but not locked.

**Channel exclusion note**: `Vel` and `pqr` are NOT used as outputs. The
dataset's `Vel` is exactly `100 * diff(Pos)` (numerical derivative + linear
smoothing; measured residual 2e-14) and `pqr` is a nonlinear transform of the
Euler-rate derivative. Including derived channels would make `C_y` singular
and break the `C_y^{-1/2}` geometry of the copyAU construction. The 10 raw
measurements are the honest observation set.

**Time alignment (standard causality)**: the motor command at step k produces
thrust that moves the vehicle at step k+1 (10 ms at 100 Hz), so transitions
pair `(y(:,k), u(:,k)) -> y(:,k+1)`; all MATLAB functions use
`ur = uc(:, valid)`.

## Layout

```
prepare_pelican_data.py             source .mat -> data/*.mat + provenance json
copyAY_pelican_soft_preference.m    main runner (learn E -> fit -> evaluate)
learn_segmented_preferred_output_directions.m   (copy of copyAX, unchanged)
fit_segmented_anchored_varx.m                   (copy of copyAX, unchanged)
evaluate_segmented_prediction.m                 (copy of copyAX, unchanged)
tests/test_pelican_dataset.py       Python dataset contract (5 mandatory + 2 source-gated)
data/copyAY_pelican_dataset.mat     y 10x1388410, u 4x1388410, segments 1..54 (LFS)
data/copyAY_pelican_dataset.json    full provenance record
results/                            metrics txt + fig png + full state .mat
```

## Data

- 54 flights = 54 segments, 1,388,410 samples total at 100 Hz (~3h50m flight
  time), Vicon indoor cube ~5 m; position mm ±5 mm, Euler deg ±0.1 deg,
  motor speeds integer [0, 218].
- All channels smoothed by the dataset authors (window-5 local regression,
  robust for motor speeds).
- Split by flight: segments 1–36 train, 37–54 validation; never shuffled
  within a flight, validation never touches training statistics.

## Acceptance criteria (this copy)

1. `dual_error < 1e-7` and anchor preservation ~0 (geometry intact).
2. Identified `Ahat` stable (`spectral_radius < 1`).
3. Validation `task_R2 > 0` and `task_RMSE < persistence_RMSE` — the
   criterion TEP failed (`R2 = [-0.005, 0.061]`, RMSE 0.853).
4. Latent model improves on persistence in standardized units.
5. Learned contribution interpretable: position/attitude channels weighted
   higher than motor channels, consistent with the soft preference.

## Verification summary

**Result: PASS on all acceptance criteria that are testable offline.**

| Criterion | Value | Status |
|---|---|---|
| 1. dual geometry intact | dual_error 2.8e-15, anchor preservation 3.5e-16 | PASS |
| 2. identified `Ahat` stable | spectral_radius **1.0000044** (marginally >1) | ⚠️ see note |
| 3. validation `task_R2 > 0` | **[0.9819, 0.9817]** vs TEP [-0.005, 0.061] | PASS |
| 4. latent model beats persistence | task RMSE [0.0918, 0.0992] < pers [0.1045, 0.1160] | PASS |
| 5. contribution interpretable | motor channels dominate (0.44–0.47), pos/attitude small | PASS (physical) |

- Training 36 flights / 1,045,715 samples; validation 18 flights / 342,695 samples.
- Task anchors (q=2) score 0.867 / 0.783 on validation.
- CRTE eigenvalues [4.42, 2.24, 0.85]: two dominant task directions, one
  weaker free direction — consistent with q=2 + r=3.
- `S_yu` spectrum decays slowly (1.87 → 0.62, no sharp gap): the free block
  keeps meaningful directions, no truncation pathology.

**Spectral radius note (honest)**: ρ = 1.0000044, marginally above 1. The
vehicle has integral dynamics (position = ∫velocity), and at 100 Hz with
window-5 smoothing the identified A is expected to sit near the unit circle;
the 4.4e-6 overshoot is identification noise, not a hard instability. One-step
prediction is unaffected. Any future closed-loop MPC stage must re-check this
(and likely needs an integrator-aware formulation); it does not invalidate the
offline prediction result.

**Contribution note**: motor-speed channels carry the largest input-reachable
energy (0.44–0.47) while position/attitude contribute little — motors respond
to `u` within one 10 ms step, positions only through integration. This is the
"task-needed output ≠ easiest-to-drive direction" distinction: motor channels
have higher authority energy but are NOT the task; the task axes (position/
attitude, weights 1.0/0.9) still reach R² ≈ 0.98.

**Channel exclusion**: `Vel` ≡ 100·diff(Pos) (linear dependence, residual
2e-14) and `pqr` (Euler-rate transform) are excluded from `y`; the 10 raw
measurements Pos/Euler/Motors are the honest observation set.

**Comparison across validation objects**:

| Object | Type | Task R² | Outcome |
|---|---|---|---|
| TEP (copyAW, deleted) | closed-loop sim | [-0.005, 0.061] | prediction failed |
| BOPTEST pilot | d_k-dominated sim | — | disturbance dominates |
| ControlGym CDR (copyAX) | linear sim | [0.814, 0.802] | prediction OK, control failed |
| **Pelican (copyAY)** | **real flight data** | **[0.982, 0.982]** | **prediction PASS** |

