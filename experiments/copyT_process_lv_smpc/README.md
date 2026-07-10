# copyT_process_lv_smpc

Purpose: synthetic industrial-process / batch-process validation for the recommended application direction: high-dimensional process measurements, two tracked quality variables, latent-variable PredVARX identification, and chance-constrained MPC.

This experiment does **not** modify `main/`. It is a synthetic process-like validation of the geometry and controller plumbing implied by the following LV-MPC / dynamic latent-variable papers, not a claim of real industrial-data validation.

## Literature anchors synchronized to Zotero

Target Zotero collection: `应用方向-LV-MPC-PredVARX参数系统验证` under parent collection `H34H6SQR`.

| Paper | Zotero key | What it contributes to this experiment | Parameter/system setting mirrored here |
|---|---:|---|---|
| Golshan 2010 — LV-MPC batch trajectory tracking | `R44SE7UE` | Batch-wise high-dimensional trajectories, PCA/LV control, quality-trajectory tracking | high-dimensional process outputs, phase-like reference changes, tracked quality variables |
| Laurí 2010 — continuous LV-MPC | `27YF4QRZ` | Continuous-process LV-MPC, PLS dynamic matrix, long-horizon MPC, input/input-rate constraints | continuous closed-loop plant, constrained manipulated variables, reduced latent predictor |
| Laurí 2013 — validity-constrained LV-MPC | `B5XK9AB4` | Keep future predictions inside the valid latent/data domain via hard validity constraints | logged QP residual `max(AU-b)` and latent/control-aware coverage checks |
| Zhu 2025 — QbD-LV-MPC | `FXCJXTKE` | Batch quality target, design-space/update logic, product-quality emphasis | two tracked quality variables and chance constraints on quality outputs |
| Ying 2024/2025 — AR-DLV | `X5HRFPD9` | High-order dynamic latent variables, TE/CSTR/semiconductor-style process monitoring | dynamic latent model with high-dimensional measurements and stochastic residual covariance |
| Yu 2024 — Kernel LaVAR / Kernel PredVAR | `ZQ87ITDG` | Nonlinear latent VAR lineage; predictable DLVs for monitoring/control bridge | linear synthetic variant here, with README noting kernelization as next extension |

## Synthetic process-like settings

The synthetic plant is deliberately shaped like process-control data:

- high-dimensional observations: `p = 30` sensor/quality variables;
- low-dimensional true latent state: `n = 6`;
- manipulated variables: `m = 3`;
- tracked quality variables: `y_1, y_2`;
- reduced control-aware latent dimension: `ell = 5`;
- offline excitation data: `T_off = 1500`;
- closed-loop run: `T_cl = 1200`;
- SMPC horizon: `N = 18`;
- input limits: `[-3, 3]`;
- quality upper chance constraint: `y_1,y_2 <= 2.6` with joint risk `alpha_joint = 0.10`.

## Main runner

```matlab
copyT_process_lv_smpc
```

Outputs are written to `results/`:

- `copyT_process_lv_smpc_data.mat`
- `copyT_process_lv_smpc_fig.png`
- `copyT_process_lv_smpc_metrics.txt`

The generated figure is aligned with the updated baseline diagnostic style and contains five rows:

1. tracked quality outputs `y_1,y_2` versus references and quality upper bound;
2. four-line noise diagnosis: estimated latent process noise, estimated observation noise, true process noise `sigma_w`, true sensor noise `sigma_e`;
3. manipulated inputs and input bounds;
4. full MPC cost `J`;
5. QP chance-constraint residual `max(AU-b)`.

## Interpretation gate

A successful synthetic validation should show:

- `tracked_projection_error` near machine precision, proving tracked quality axes are retained by the reduced subspace;
- moderate reconstruction residual despite high-dimensional observations;
- MAE/RMSE after warm-up substantially below the reference step magnitudes;
- empirical upper-constraint violation rates near zero or compatible with the risk design;
- QP success rate near 1 and fallback count 0 or justified;
- `max_recorded_qp_constraint <= 0`;
- no actuator saturation unless intentionally stressed.

## Limitation and next real-data layer

This validates the application geometry and controller plumbing on synthetic process-like data. It does **not** prove performance on real industrial data. The next validation layer should use one of:

1. Tennessee Eastman Process: high-dimensional process benchmark, close to Ying/Yu dynamic latent-variable papers;
2. CSTR: low-to-medium dimensional sanity benchmark;
3. IndPensim / penicillin fermentation: closest to Zhu 2025 QbD-LV-MPC batch-quality narrative;
4. Data-center thermal process: closest to high-dimensional engineering MIMO deployment.
