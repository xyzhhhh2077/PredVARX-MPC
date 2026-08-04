# copyBJ: copyBI with online covariance learning

This experiment keeps the copyBI Pelican model, task geometry, reference, input
box, horizon, objective, and disturbance seeds fixed. Its only algorithmic
change is the copyAU-style online update of uncertainty statistics.

## Update loop

At every 100 Hz plant sample, copyBJ computes the closed-loop latent innovation

```text
e_k = z_k - A z_aug,k-1 - B (u_k-1 - u_mean)
```

and retains the newest 40 samples. After at least five samples, it replaces the
SMPC innovation covariance with the centered sample covariance plus a `1e-8`
PSD floor. At every 10-sample control decision, the complete 18-stage chance
tightening is recomputed from that covariance and supplied to the already
compiled parameterized QP.

The projection-residual proxy is also updated with a 40-sample window for the
copyAU diagnostic contract. It does not enter the tracked-axis main QP.

## Full-length result

The direct copyAU-style update (`online_weight = 1`) is not feasible for the
complete copyBI event. At sample 8,490 (84.9 s), the 40-sample covariance gives
an 18-step tightening larger than the `3.8` task half-width, so the mean chance
interval is empty. The run stops without fallback and preserves this failure.

An explicitly separate regularized candidate uses

```text
Sigma_effective = 0.2 Sigma_offline + 0.8 Sigma_online,40
```

This candidate completes 18,000 steps with zero hard-bound violations and zero
fallbacks. Its minimum hard margin is `0.155262`; axis RMSE is
`[0.515811, 0.382387, 0.111666]`. The frozen copyBI baseline has minimum margin
`0.007150` and RMSE `[0.461697, 0.321477, 0.091572]`. Thus the candidate backs
off farther from the boundary but tracks less accurately. Fourteen OSQP solves
reached the ordinary iteration limit and succeeded after a cold-start resolve
of the same QP with a higher iteration limit; no substitute control was used.

## Fixed and updated objects

| Object | Treatment |
|---|---|
| `E, A, B, P, R` | Frozen from the offline Pelican identification |
| `Sigma_eps` | Online 40-sample covariance update |
| `Sigma_obs_proxy` | Online 40-sample diagnostic update |
| Physical wind covariance | Not supplied to the controller |
| Reference and seeds | Identical to copyBI |

Because physical wind is injected into the plant but omitted from the prediction
model, the observed online innovation is a composite residual. It absorbs the
sampled identified innovation, injected wind effects, and any closed-loop model
mismatch. It is not a separated estimate of physical process noise or sensor
measurement noise.

## Run

```bash
python experiments/copyBJ_pelican_online_covariance_smpc/run_online_covariance_smpc.py
python -m pytest tests/test_copybj_pelican_online_covariance_smpc.py -q
```

Artifacts are written to `results/`:

- `copyBJ_online_covariance_summary.json`
- `copyBJ_online_covariance_comparison.npz`
- `copyBJ_online_covariance_comparison.png`

## BI-style four figures + 200-frame GIF

After the verified online run exists, build the same evidence pack as copyBI from
the stored beta=0.8 trajectory (no controller re-run):

```bash
cd experiments/copyBJ_pelican_online_covariance_smpc
python generate_smpc_figures.py
python generate_smpc_summary_gif.py
python verify_smpc_summary_gif.py
```

Artifacts in `results_smpc/`:

- `copyBJ_smpc_xyz_reference.png`
- `copyBJ_smpc_3d_trajectory.png`
- `copyBJ_smpc_noise_timeseries.png`
- `copyBJ_smpc_stage_cost.png`
- `copyBJ_smpc_figure_data.npz`
- `copyBJ_smpc_figure_summary.json`
- `copyBJ_smpc_summary_200frames.gif`
- `copyBJ_smpc_summary_keyframes.png`
- `copyBJ_smpc_summary_gif_audit.json`

## Claim boundary

This is online uncertainty adaptation on a frozen identified Pelican model. It
is not online identification of `E/A/B/P/R`, not process/measurement-noise
separation, not probability calibration, and not real-flight validation. The
`0.8/0.2` shrinkage is an experimental feasibility regularization, not part of
the copyAU theory and not a general recursive-feasibility result.
