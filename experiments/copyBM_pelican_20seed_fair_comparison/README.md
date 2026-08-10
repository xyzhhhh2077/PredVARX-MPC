# copyBM: fair Frozen / copyBJ / copyBL comparison

This is a same-disturbance paired benchmark. For each fixed pair
`innovation_seed=16+i, wind_seed=7016+i`, the driver creates one reference,
innovation, physical wind velocity, and wind latent increment, then passes
numerically identical copies to Frozen (`BH.simulate_controller`), copyBJ
(`BJ.simulate_online_controller`, shrinkage `online_weight=0.8`), and copyBL
(`BL.simulate_dual_gated_controller`). All use `initial_task=np.zeros(3)`.
The disturbance and equilibrium initial-state SHA-256 values are recorded and
the driver asserts the three disturbance hashes agree.

The 20 pairs are exploratory evidence, not Monte Carlo probability calibration.
Frozen/BJ have no gate, so their gate acceptance rate is `null`; BJ covariance
update attempts/counts are reported separately. BL gate counts come from its
event log, with value acceptance rate defined as
`value_accept_count/value_trigger_count`. BJ update rate and BL gate acceptance
rate are different quantities and must not be compared as if they were the
same.

Metrics are xyz RMSE, cumulative realized stage cost using
`BI.realized_stage_cost(S, reference[:completed], U, model["u_mean"])`, minimum
hard margin, QP count/failures/rate, and the controller-specific gate fields.
An incomplete or failed trajectory gets `null` for full-trajectory metrics; it
is never presented as a truncated success. Single-seed failures remain in
`results/per_seed.json` and are not silently deleted.

Run a short check with:

```text
.venv/bin/python run_fair_comparison.py --steps 120 --max-seeds 1
```

The full requested batch is `--max-seeds 20` and was not run during this
implementation turn. Results are atomically resumable after every successful
seed; resume rejects a changed configuration hash or seed-pair list.