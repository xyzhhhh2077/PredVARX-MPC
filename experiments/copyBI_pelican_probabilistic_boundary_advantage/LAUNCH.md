# copyBI - Pelican near-boundary SMPC artifacts

This launch note covers the retained SMPC-only figures and synchronized GIF.
See `README.md` for the data/model contract, numerical results, and claim
boundaries. Historical deterministic-MPC comparison outputs remain available,
but they are not the current presentation result.

## Run

From the repository root:

```bash
python -u experiments/copyBI_pelican_probabilistic_boundary_advantage/generate_smpc_figures.py
```

The current SMPC-only outputs are written to `results_smpc/`:

- `copyBI_smpc_xyz_reference.png` — x/y/z task and reference time series
- `copyBI_smpc_3d_trajectory.png` — 3-D task and reference trajectory
- `copyBI_smpc_noise_timeseries.png` — composite latent innovation and physical wind
- `copyBI_smpc_stage_cost.png` — realized one-step tracking and input cost
- `copyBI_smpc_figure_data.npz` — plotted arrays
- `copyBI_smpc_figure_summary.json` — run metrics and plot definitions

To combine those data into a synchronized 200-frame GIF and audit the encoded
artifact:

```bash
python -u experiments/copyBI_pelican_probabilistic_boundary_advantage/generate_smpc_summary_gif.py
python experiments/copyBI_pelican_probabilistic_boundary_advantage/verify_smpc_summary_gif.py
```

This additionally writes:

- `copyBI_smpc_summary_200frames.gif`
- `copyBI_smpc_summary_keyframes.png`
- `copyBI_smpc_summary_gif_audit.json`

The top row places the task-coordinate and 3D panels side by side. Three
full-width rows below show the four standardized control inputs with their
`+/-7` box constraints, the two noise objects, and the realized stage cost. The
resulting animation is 1152x1035 pixels.

The noise panels do not claim separately identified process-noise and
measurement-noise covariances. They show the sampled composite latent
innovation and the independently injected physical-wind realization.

## Retained SMPC run

- steps: 18,000
- input cap: +/-7
- reference pressures: x/y/z = 3.46/3.46/3.46 standardized units
- copyBI experiment hard bound: x/y/z = +/-3.80 standardized units
- wind seed: 7016
- innovation seed: 16

## Interpretation boundary

The SMPC chance covariance remains the identified `Sigma_eps_aug`. The physical
wind is injected into both simulated plants but is not added to that covariance.
Consequently, the retained run demonstrates one complete SMPC trajectory with
zero hard-bound violations and zero fallback calls. It does not establish
wind-risk calibration, recursive feasibility, closed-loop stability, or a
robust hard guarantee across random seeds.
