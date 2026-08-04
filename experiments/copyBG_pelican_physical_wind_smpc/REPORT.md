# copyBG verification report

This is a **model-in-the-loop** comparison built on the model identified from
the public Waterloo Pelican flight dataset. It is not a real-flight or
hardware-in-the-loop experiment.

## Run configuration (default: rising circular helix)

- Closed-loop length: 12000 samples at 100 Hz (120 s)
- Control interval: 1 sample
- Prediction horizon: N = 18 (0.18 s)
- Control weights: Q = 80 on the tracked xyz error, Ru = 0.18
- Joint chance budget: 0.05 (SMPC only)
- Normalized input bound: [-6, 6] on each of the four inputs
- Hard bound on tracked position (standardized): 3.4
- Reference mode: **`circle-helix`** (default)
  - Geometry: bottom dwell 24 s at z = −3.0, three rising circular turns
    `(rx cos θ, ry sin θ)` with `rx = ry = 2.0` and z linear −3.0 → +3.0,
    then top dwell 24 s. Horizontal section satisfies
    `x²/rx² + y²/ry² = 1` (smooth cylinder wall; no 90° corners).
  - Pressures: `(2.0, 2.0, 3.0)` inside hard bound 3.4
- Optional modes still available: `--reference square-spiral` (L∞ edge stress
  at 3.2) and `--reference tour` (C1 axis tour). Prior square-spiral
  calibrated artifacts are kept under `results_square_spiral/`.
- Closed-loop physical wind: independent Gaussian with std
  `sigma_wind_mps = 0.02 m/s`, fixed seed
- Latent innovation: drawn from the identified `Sigma_eps`, fixed seed
- Both controllers receive the same innovation array and the same wind
  realization, verified by shared SHA-256
  `375c6ac93ac6b89b43013b3bf23748534f4e037a3a76dc80d0a76f7040512583`
- QP solver: OSQP with `eps_abs = eps_rel = 1e-8` and polishing

## Why the default changed

History lookup (see skill
`predvarx-mpc/references/pelican-reference-shape-history-lookup.md`) showed:

| remembered “good cylindrical” run | actual shape | RMSE |
|---|---|---|
| copyBE | yz ring, x fixed | ~0.05, 0 viol |
| copyBF tour | cylinder-surface C1 tour | ~0.10, 0 viol |
| copyBF tightened spiral | rising elliptical helix | **failed** (unit bug, not controller) |
| copyBG square spiral | L∞ edge corners | x/y RMSE ~2.2; MPC dies @10848 |

User request: restore rising **circular** helix under the copyBG physical-wind
contract; keep square spiral as optional pressure mode.

## Initialization diagnostics

| Check                                | Value           |
|--------------------------------------|-----------------|
| Dynamic residual norm                | 2.66e-15        |
| Tracked output residual norm         | 9.77e-15        |
| Maximum absolute steady input        | 3.553           |

The spiral starts at `(2.0, 0, −3.0)`. Equilibrium input is well inside ±6.

## Hard-bound contact and SMPC tightening

| Quantity                                                     | Value              |
|--------------------------------------------------------------|--------------------|
| Hard bound (standardized)                                    | 3.4                |
| y-scale used for unstandardization (m), per axis              | `[0.909, 0.919, 0.887]` |
| Hard position upper bound (m) = 3.4 · y_scale                | `[3.090, 3.124, 3.016]` |
| SMPC terminal tightening (standardized), per axis            | `[0.0376, 0.0362, 0.0354]` |
| SMPC terminal mean-bound (standardized), per axis            | `[3.362, 3.364, 3.365]` |

## Closed-loop results (σ_w = 0.02, circle helix)

| Metric                                                 | SMPC             | Deterministic MPC |
|--------------------------------------------------------|------------------|-------------------|
| Completed steps                                        | **12000**        | **12000**         |
| QP failure step / status                               | none             | none              |
| QP failure count                                       | 0                | 0                 |
| Fallback count                                         | 0                | 0                 |
| Active predicted-bound steps                           | 36               | 29                |
| Positive-dual predicted-bound steps                    | 36               | 29                |
| Hard violation rate                                    | **0.0**          | **0.0**           |
| Minimum hard margin (std), per axis                    | `[0.335, 0.475, 0.0317]` | `[0.335, 0.475, 0.00124]` |
| Peak \|x\|, \|y\|, \|z\| realized (std)                | `[3.065, 2.925, 3.368]` | `[3.065, 2.925, 3.399]` |
| RMSE standardized `[x, y, z]`                          | `[1.529, 1.732, 0.300]` | `[1.529, 1.732, 0.300]` |
| Input saturation steps                                 | 0                | 0                 |
| Max QP constraint residual                             | 1.4e-14          | 2.6e-14           |
| Mean QP solve time                                     | 1.28e-3 s        | 1.29e-3 s         |

Interpretation:

- **Both controllers finish the full 120 s** with zero hard violations and zero
  fallback. On the square spiral (same Q/Ru/N/dt/wind), MPC became infeasible
  at step 10848; the circular path removes the 90° corners that forced that
  failure under identical disturbance seeds.
- SMPC keeps a larger z hard margin (0.032 vs 0.001) and a lower z peak
  (3.368 vs 3.399). Chance tightening is active for 36 steps; the mean-bound
  guide is visible on the z channel near the dwells.
- x/y RMSE ~1.5–1.7 is still large relative to the reference radius 2.0: the
  identified slow x/y channels (τ ≈ 20–25 s) lag a multi-turn helix under
  fixed Q/Ru/N. This is **not** a wind artifact and not a unit bug — both
  controllers share the lag almost identically. z tracking (RMSE 0.30) is
  the tightest axis and the one that presses the hard bound.
- Do **not** claim RMSE comparable to copyBE’s 0.05 ring: that run used a
  much smaller yz circle (r≈0.55) over 12 s. This run is a 120 s rising
  helix at radius 2.0 with independent physical wind.

## Comparison to square-spiral pressure mode (archived)

Artifacts: `results_square_spiral/` (σ_w = 0.02, pressures 3.2).

| Metric | circle-helix default | square-spiral archive |
|---|---|---|
| SMPC completed | 12000 | 12000 |
| MPC completed | **12000** | 10848 (infeasible) |
| SMPC hard viol | 0 | 0 |
| SMPC RMSE xyz | 1.53 / 1.73 / 0.30 | 2.20 / 2.38 / 0.35 |
| SMPC active tight steps | 36 | 199 |

Square spiral remains the correct **stress** tool when the goal is to drive
deterministic MPC infeasible while SMPC survives. Circle helix is the default
when the goal is long cylindrical tracking under physical wind.

## Channel contract (unchanged)

- `Sigma_eps`: offline latent innovation from flights 1–36; enters QP via
  total process covariance.
- `Sigma_obs_proxy`: diagnostic only; never enters the QP
  (`uses_true_Sigma_n = 0`).
- Physical wind: independent Gaussian velocity (m/s) → `Δp = v·dt` →
  standardized → `G_w` into latent increment. Separate RNG stream from
  innovation.
- Controller state: reconstructed from outputs only
  (`R.T @ (y − y_mean)` + previous latent block).

## Verification

```text
python experiments/copyBG_pelican_physical_wind_smpc/verify_results.py
→ status PASS
  smpc_completed_steps 12000
  deterministic_mpc_completed_steps 12000
  smpc_active_qp_steps 36
  smpc_hard_violation_rate 0.0
```

Unit tests: `pytest tests/test_copybg_physical_wind_smpc.py` → 17 passed.

## Honest limits

1. Model-in-the-loop on a frozen identified model + synthetic physical wind —
   not real-flight validation.
2. x/y lag under multi-turn helix with fixed Q/Ru/N is structural; reducing it
   would require retuning the frozen control contract (not done here).
3. At pressures (2,2,3) both controllers stay feasible; the SMPC-vs-MPC
   **safety separation** is milder than on the square spiral (margin gap on z,
   not MPC death). For a hard infeasibility contrast, use
   `--reference square-spiral`.
4. No-wind ablation (same circle helix, `results_nowind/`) isolates the wind
   channel; numbers below.

## No-wind ablation (σ_w = 0.0, same circle helix)

Shared disturbance SHA-256 prefix: `08f167a39874c447…` (innovation only; wind array is zero).

| Metric | SMPC | Deterministic MPC |
|---|---|---|
| Completed steps | **12000** | **12000** |
| QP failure / fallback | 0 / 0 | 0 / 0 |
| Active predicted-bound steps | **0** | **0** |
| Hard violation rate | **0.0** | **0.0** |
| Peak \|z\| (std) | 3.192 | 3.192 |
| Min z hard margin (std) | 0.208 | 0.208 |
| RMSE standardized `[x, y, z]` | `[1.408, 1.608, 0.264]` | `[1.408, 1.608, 0.264]` |

No-wind takeaways:

- Both controllers are **identical** to machine precision on the tracked path
  when wind is off (same peaks, margins, RMSE). Chance tightening never
  activates (`active_qp_steps = 0`) because the mean trajectory stays well
  inside the tightened mean-bound.
- x/y lag (RMSE ~1.4–1.6) remains with **no wind**, confirming it is the slow
  identified channels vs multi-turn helix under fixed Q/Ru/N — not a wind
  artifact.
- Wind σ_w = 0.02 raises z peak 3.192 → 3.368 (SMPC) / 3.399 (MPC), activates
  SMPC tightening (36 steps), and opens a z-margin gap (0.032 vs 0.001). That
  is the physical-wind channel doing its job.
