# copyBF: long 3D boundary stress with external disturbances

This experiment reuses the frozen copyBE second-order Pelican model and compares
SMPC (`alpha_joint=0.05`) with deterministic MPC under exactly the same residual
innovation samples and prescribed external task-space gusts.

## Scope

The result is a 120 s presentation assembled from three independent 40 s
boundary scenarios. The controller/plant state is reset to a dynamically
consistent equilibrium at 40 s and 80 s; it is not one continuous 120 s flight.
This avoids treating an unreachable transfer between opposite boundary regimes
as controller failure.

- x: reachable-limit pressure. The identified model cannot hold the x hard
  boundary with the frozen input bounds, but the reference lies beyond the
  SMPC tightened boundary so the chance constraint is active.
- y: near-boundary pressure with zero-fallback calibration.
- z: near-boundary pressure calibrated below its no-fallback threshold. The
  chance constraint is active, but the tested realized trajectory does not show
  a larger minimum hard-bound margin than deterministic MPC.

The prescribed gust uses an outward pulse followed by an equal recovery pulse,
so each cycle has zero net impulse and does not create artificial long-term
drift. It is separate from the fitted residual innovation covariance.
Both controllers receive identical complete disturbance sequences. This is a
frozen-model-in-the-loop validation, not real-flight evidence and not a proof of
recursive feasibility or probabilistic calibration.

The x/y scenarios are deliberately strong pressure tests: realized input
saturation occupies roughly 71%-90% of samples across seeds. Therefore their
larger SMPC hard-bound margins demonstrate constraint activation under this
specific saturated regime, not general closed-loop superiority. In the z
scenario SMPC is active at every QP update, but its realized minimum hard-bound
margin is about 0.006 smaller than deterministic MPC for the tested seeds.

## Frozen controller contract

- latent dimension `d=10`, second-order companion state
- prediction horizon `N=18`
- `Q=80`, `Ru=0.18`
- normalized input increment bounds `[-6, 6]`
- common hard task bound from copyBE
- 100 Hz sampling

## Run and verify

```bash
python experiments/copyBF_pelican_long_3d_boundary_disturbance/run_long_3d_boundary_disturbance.py
python experiments/copyBF_pelican_long_3d_boundary_disturbance/verify_results.py
python -m pytest tests/test_copybf_pelican_long_boundary_disturbance.py -q
```

Outputs are written to `results/` as JSON metrics, NPZ arrays, and a PNG figure.
The PNG marks the two scenario resets explicitly.

`plot_3d_trajectories.py` additionally renders the three independent scenarios
in 3D task coordinates. This visualization must not be interpreted as one
continuous spatial flight because the state is reset between panels.

`run_continuous_3d_tightened_spiral.py` is the separate no-reset spatial-motion
experiment. Its 120 s reference moves smoothly along a rising elliptical helix
that approaches the axis-specific terminal tightened bounds. The external
wind channel is a per-step independent Gaussian process noise
``w_k ~ N(0, sigma_w^2 I_ell)`` injected in latent coordinates
(sigma_w = 0.045, fixed wind_seed = 7) — a copyT-style stochastic channel that
is independent of the identified residual innovation. The deterministic
zero-net-impulse gust helper has been removed; the figure title and result
metadata explicitly describe the channel as Gaussian process wind. Outputs
include JSON, NPZ, PNG, and a MAT file
`copyBF_continuous_3d_tightened_spiral.mat` (reference, SMPC/MPC trajectory
and control, residual noise, process wind, hard/tightened bounds, sigma_w,
seed, wind_seed, and core scalar metrics).
