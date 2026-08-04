# copyBJ: Pelican online covariance-adaptive SMPC

copyBJ keeps the **copyBI** Pelican model, task geometry, reference, input box,
horizon, objective, and disturbance seeds fixed. The only algorithmic change is a
**copyAU-style online update of the composite latent-innovation covariance**
`Σ_ε`. Geometry `E, A, B, P, R` stays frozen.

This is closed-loop simulation on a frozen real-data model
(**model-in-the-loop**). It is not online re-identification, not
process/measurement-noise separation, not probability calibration, and not
hardware or flight validation.

## Best synchronized GIF (β = 0.8 online main trajectory)

![Synchronized copyBJ Pelican online SMPC run](results_smpc/copyBJ_smpc_summary_200frames.gif)

The animation is **200 frames / 18 s** (90 ms per frame) and compresses the
180-second experiment by ×10. Layout (same evidence pack as copyBI):

1. synchronized `xyz` tracking + 3-D reference/trajectory (hard box ±3.8)
2. four standardized control commands and absolute bounds
3. composite latent innovation and injected three-axis physical wind
4. realized one-step stage cost

Keyframe contact sheet (frames 1 / 67 / 134 / 200):

![copyBJ summary keyframes](results_smpc/copyBJ_smpc_summary_keyframes.png)

## What is updated

At every 100 Hz plant sample, copyBJ forms the closed-loop latent innovation

```text
e_k = z_k - A z_aug,k-1 - B (u_{k-1} - u_mean)
```

and keeps the newest **40** samples. After ≥5 samples it replaces the SMPC
innovation covariance with the centered sample covariance plus a `1e-8` PSD floor.
At every 10-sample control decision, the full 18-stage chance tightening is
recomputed from that covariance and fed to the already compiled parameterized QP
(no rebuild).

The projection-residual proxy `Σ_obs` is updated with the same 40-sample window
for diagnostics only. It does **not** enter the tracked-axis main QP.

Physical wind is injected into the plant but **not** given to the controller as a
calibrated risk term. Online `e_k` is therefore a **composite residual**
(identified innovation + wind effects + closed-loop mismatch), not a separated
`Q_w` / `R_v` / wind covariance.

## Two candidates (must be reported separately)

### A. Raw copyAU-style (`online_weight = 1`)

| Quantity | Value |
|---|---:|
| Completed samples | **8,490 / 18,000** (stops at 84.9 s) |
| Stop reason | 18-step chance tightening empties the mean interval (fail-closed) |
| Hard-bound violations | 0 |
| Fallback calls | 0 |
| Min hard margin before stop (std) | 0.2855 |

The short-window online covariance can push tightenings past the box width on this
long boundary trajectory. The failure is retained; the controller is **not**
silently continued.

### B. Shrinkage candidate (experimental regularizer, **not** a theorem)

```text
Σ_eff = 0.2 Σ_offline + 0.8 Σ_online,40
```

| Quantity | Frozen copyBI | Online copyBJ β=0.8 |
|---|---:|---:|
| Completed samples | 18,000 | 18,000 |
| Hard-bound violations | 0 | 0 |
| Fallback calls | 0 | 0 |
| Min hard margin (std) | **0.007150** | **0.155262** |
| RMSE xyz (std) | [0.4617, 0.3215, 0.0916] | [0.5158, 0.3824, 0.1117] |
| Chance-active QP steps | 564 | 1,403 |
| Input-cap samples | 7,400 | 7,130 |
| OSQP same-QP cold-start retries | 0 | 14 (all succeeded; no substitute control) |

On this single shared seed the shrinkage run stays **farther from the box**
(more conservative) and tracks **slightly worse**. The `0.8/0.2` mix is an
**experimental feasibility regularizer** — not part of copyAU theory and **not** a
recursive-feasibility / stability / probability-calibration proof.

Online `Σ_ε` RMS std range ≈ `0.01227`–`0.02406` (start ≈ `0.01886`, end ≈ `0.01917`).

## Frozen vs online comparison figure

![Frozen copyBI vs online copyBJ](results/copyBJ_online_covariance_comparison.png)

## BI-style four static panels (β = 0.8)

### xyz + reference + hard bounds

![xyz tracking](results_smpc/copyBJ_smpc_xyz_reference.png)

### 3-D trajectory + hard box

![3d trajectory](results_smpc/copyBJ_smpc_3d_trajectory.png)

### Composite innovation + physical wind

![innovation and wind](results_smpc/copyBJ_smpc_noise_timeseries.png)

### Realized stage cost

![stage cost](results_smpc/copyBJ_smpc_stage_cost.png)

Realized stage cost is `Q‖s−r‖² + R‖u−ū‖²` (one-step diagnostic), **not** the full
finite-horizon QP objective.

## Fixed vs updated objects

| Object | Treatment |
|---|---|
| `E, A, B, P, R` | Frozen offline Pelican identification |
| `Σ_ε` | Online 40-sample covariance |
| `Σ_obs` proxy | Online 40-sample **diagnostic only** |
| Physical wind covariance | **Not** supplied to the controller |
| Reference, seeds, `Q_y`, `H`, horizon | Identical to copyBI |

Shared disturbance seed SHA256:
`7671c9ab8a456370278ed9052267a126890699a23d8731eeb42e38fe263b3375`

## Reproduce

From the repository root:

```bash
python experiments/copyBJ_pelican_online_covariance_smpc/run_online_covariance_smpc.py
python -m pytest tests/test_copybj_pelican_online_covariance_smpc.py -q
```

Controller artifacts → `results/`:

- `copyBJ_online_covariance_summary.json`
- `copyBJ_online_covariance_comparison.npz`
- `copyBJ_online_covariance_comparison.png`

BI-style pack from the stored **β=0.8** trajectory (**no controller re-run**):

```bash
cd experiments/copyBJ_pelican_online_covariance_smpc
python generate_smpc_figures.py
python generate_smpc_summary_gif.py
python verify_smpc_summary_gif.py
```

Pack artifacts → `results_smpc/`:

- `copyBJ_smpc_xyz_reference.png`
- `copyBJ_smpc_3d_trajectory.png`
- `copyBJ_smpc_noise_timeseries.png`
- `copyBJ_smpc_stage_cost.png`
- `copyBJ_smpc_figure_data.npz` / `copyBJ_smpc_figure_summary.json`
- `copyBJ_smpc_summary_200frames.gif`
- `copyBJ_smpc_summary_keyframes.png`
- `copyBJ_smpc_summary_gif_audit.json`

## Claim boundary

| Claim | Status |
|---|---|
| Online composite-innovation covariance adaptation on frozen Pelican model | ✅ |
| Raw 40-sample window can fail closed mid-trajectory (negative evidence) | ✅ |
| β=0.8 shrinkage completes this seed with 0 violation / 0 fallback | ✅ |
| BI-style four panels + 200-frame GIF for the β=0.8 main run | ✅ |
| Online learning of `E/A/B/P/R` | ❌ |
| Process / measurement / wind noise separation | ❌ |
| Probability calibration or recursive feasibility theorem | ❌ |
| Outdoor / hardware flight validation | ❌ |

## Related

- Baseline event: `experiments/copyBI_pelican_probabilistic_boundary_advantage/`
- Unit tests: `tests/test_copybj_pelican_online_covariance_smpc.py` (5 passed)
