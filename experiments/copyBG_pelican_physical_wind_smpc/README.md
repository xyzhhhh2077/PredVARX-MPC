# copyBG: Pelican Physical-Wind SMPC

This experiment is a **model-in-the-loop** comparison on a model identified
from the public Waterloo Pelican flight dataset. It is not real-flight or
hardware-in-the-loop validation.

## Scope

- Offline identification and noise estimation use Pelican training flights
  1-36 only.
- `Sigma_eps` is the covariance of the latent second-order dynamic innovation.
  It combines unmodeled dynamics, fitting error, and measurement-noise effects;
  it is not a separately identified physical process-noise covariance.
- `Sigma_obs_proxy` is the covariance of the output projection residual. It is
  saved for diagnostics only and never enters the controller QP.
- Closed-loop wind is an independent Gaussian velocity channel in meters per
  second. Each sample is multiplied by `dt=0.01 s`, converted to standardized
  xyz displacement, and mapped into the current latent innovation with `G_w`.
- No true sensor-noise covariance is used (`uses_true_Sigma_n=0`). The
  experiment does not identify separate `Q_w` and `R_v` objects.

## Control Contract

Both controllers use the same identified innovation array and the same physical
wind realization, verified by a shared SHA-256 digest. The controller state is
reconstructed from closed-loop output with `R.T @ (y-y_mean)` plus the previous
reconstructed latent block; the controller does not receive the plant's hidden
state or disturbance.

The fixed control contract is 100 Hz, one-sample control interval, `N=18`
(`0.18 s`), `Q=80`, `Ru=0.18`, and normalized input bounds `[-6,6]`. The SMPC
uses a joint chance budget of `0.05`; deterministic MPC keeps hard predicted
constraints but no chance tightening. A non-optimal QP stops that controller
immediately. There is no fallback and no synthesized post-failure trace.

## Reference paths

Default reference is a **rising circular helix** (120 s, 12000 samples at
100 Hz): 24 s bottom dwell at `z = -3.0`, three rising circular turns
`(rx cos θ, ry sin θ)` with `rx = ry = 2.0` and `z` rising linearly to `+3.0`,
then 24 s top dwell. This restores the copyBF circle-helix geometry (smooth
cylinder wall, no 90° corners) under the copyBG physical-wind channel contract.

Optional modes:

| `--reference` | Default pressures | Role |
|---|---|---|
| `circle-helix` (default) | `(2.0, 2.0, 3.0)` | Trackable cylindrical rising helix |
| `square-spiral` | `(3.2, 3.2, 3.2)` | L∞ edge corner-stress pressure test |
| `tour` | `(1.5, 1.5, 3.2)` | C1 center→+x→+y→+z→center tour |

Prior square-spiral calibrated artifacts are kept under
`results_square_spiral/`. Live default outputs go to `results/`.

## Run

```bash
python experiments/copyBG_pelican_physical_wind_smpc/run_physical_wind_smpc.py
python experiments/copyBG_pelican_physical_wind_smpc/verify_results.py
```

Useful flags: `--steps 12000`, `--sigma-wind-mps 0.02`,
`--reference {circle-helix,square-spiral,tour}`, `--pressures PX PY PZ`,
`--hard-bound 3.4`, `--output-dir PATH`.

The `results/` directory contains JSON, compressed NPZ, MATLAB MAT, and separate
SMPC and deterministic-MPC PNG diagnostics. Both figures show xyz position and
reference, hard bounds, the SMPC terminal mean-bound guide, the 3D trajectory,
tracking errors, all four inputs, and realized hard margin/activation.
