# copyBL: Pelican dual-gated covariance SMPC

copyBL is an independent copy of the copyBJ Pelican experiment. It keeps the
identified geometry and dynamics `E, A, B, P, R`, objective, horizon, reference,
and disturbance seeds fixed. Only the composite closed-loop latent-innovation
covariance `Sigma_eps` is eligible for an online update.

Each covariance candidate is transactional:

1. **Feasibility gate (D3RMPC-inspired engineering proxy):** test the 18-stage
   chance interval and the current parameterized QP. If needed, backtrack along
   `Sigma(eta) = (1-eta) Sigma_acc + eta Sigma_cand`.
2. **Control-value gate (Lin-inspired paired closed-loop realized-cost proxy):** every 400 plant
   samples, replay the next 180 stored simulation samples from the same state
   and previous control, with identical future innovation, wind, and the true
   global future reference horizon, under the accepted and feasible candidate
   covariances. Commit only if both replays complete, neither violates the hard
   bound, and the candidate reduces paired realized cost by more than the 1%
   deadband.

The feasibility gate runs at every 10-sample control decision. A feasible
candidate remains pending until the scheduled value gate accepts it. Rejection
always preserves the last accepted online covariance; the frozen covariance is
used only for cold-start initialization.

## Artifacts

The completed run writes controller data and an explicit event log to `results/`:

- `copyBL_dual_gated_summary.json`
- `copyBL_dual_gated_run.npz`
- `copyBL_dual_gate_events.json`

The evidence pack in `results_media/` contains exactly four static figures:

1. `copyBL_smpc_xyz_reference.png`
2. `copyBL_smpc_3d_trajectory.png`
3. `copyBL_smpc_noise_timeseries.png`
4. `copyBL_smpc_stage_cost_dual_gates.png`

The synchronized animation is `copyBL_dual_gate_summary_200frames.gif`. Its
status line shows the most recent feasibility/value decision, including value
rejections and committed updates.

## Reproduce

From the repository root with NumPy, SciPy, CVXPY/OSQP, h5py, Matplotlib,
Pillow, and pytest installed:

```bash
python experiments/copyBL_pelican_dual_gated_covariance_smpc/run_dual_gated_covariance_smpc.py
python experiments/copyBL_pelican_dual_gated_covariance_smpc/generate_figures.py
python experiments/copyBL_pelican_dual_gated_covariance_smpc/generate_summary_gif.py
python experiments/copyBL_pelican_dual_gated_covariance_smpc/verify_artifacts.py
python -m pytest tests/test_copybl_pelican_dual_gated_covariance_smpc.py -q
```

## Claim boundary

The value metric is a **Lin-inspired paired closed-loop realized-cost proxy**,
not Lin et al.'s strict prediction power and not an optimal-policy comparison.
The future disturbance segment is available because this is an offline,
single-seed model-in-the-loop feasibility experiment; it is not a deployable
online oracle. The D3RMPC-inspired backtracking gate does not transfer that
paper's recursive-feasibility, stability, or safety theorem; it is only a
feasibility/backtracking engineering proxy. `Sigma_eps` is only a composite
closed-loop latent innovation covariance. No probability
calibration, separated process/measurement/wind covariance, hardware, or flight
validation is claimed.
