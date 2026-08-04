# copyBI: Pelican near-boundary SMPC visual audit

copyBI is the current presentation endpoint of the copyBA--copyBI Pelican
experiment line. It loads the frozen second-order latent VARX model identified
from the Waterloo WAVELab AscTec Pelican dataset, closes the SMPC loop on a new
six-face reference, and injects a new physical-wind realization into the
simulated plant.

This is a closed-loop simulation with a frozen real-data model. It is not a
replay of a recorded flight and is not hardware or flight validation.

## Data and model contract

- Source: [Waterloo WAVELab Pelican dataset](https://github.com/wavelab/pelican_dataset)
- Sampling: 100 Hz, 54 complete indoor Vicon flights
- Split: flights 1--36 for identification; flights 37--54 held out
- Measured output: `y in R^10` = position (3, m), Euler angles (3, rad), and
  measured motor speeds (4)
- Manipulated input: `u in R^4` = commanded motor speeds
- Latent coordinate: `z in R^5` = three fixed standardized `xyz` task
  directions and two free statistical directions
- MPC state: `[z_k; z_{k-1}] in R^10` for the second-order latent VARX model

The two free latent directions are statistical subspace coordinates and are not
assigned one-to-one physical names.

## Retained SMPC run

| Quantity | Value |
|---|---:|
| Plant samples | 18,000 (180 s) |
| MPC update interval | 10 samples |
| QP decisions | 1,800 |
| Chance-active QP decisions | 564 (31.3%) |
| Standardized task hard bound | `+/-3.80` |
| Standardized reference pressure | `3.46` on x, y, and z |
| Standardized input box | `-7 <= u_i <= 7` |
| Hard-bound violations | 0 |
| Fallback calls | 0 |
| Minimum hard-bound margin | `0.00715` |
| Input-cap samples | 7,400 of 18,000 (41.1%) |
| Mean realized stage cost | `38.8386` |
| Maximum realized stage cost | `670.9595` |

The standardized `+/-7` command box is inherited as an experiment-design
constraint. It is not a manufacturer hardware rating and partly extends beyond
the command support observed in the identification data.

## Best synchronized GIF

![Synchronized copyBI Pelican SMPC run](results_smpc/copyBI_smpc_summary_200frames.gif)

The 1152x1035 animation contains 200 frames and plays for 18 seconds. It
compresses the 180-second experiment by a factor of ten. Its four rows show:

1. synchronized `xyz` tracking and the 3-D reference/trajectory;
2. four standardized control commands and their absolute `+/-7` bounds;
3. the composite latent innovation and injected three-axis physical wind;
4. the realized one-step stage cost.

The realized cost shown in the animation is
`80 ||s-r||_2^2 + 0.18 ||u-u_mean||_2^2`. It is a one-step diagnostic, not the
full finite-horizon QP objective.

## Noise interpretation

`Sigma_eps` is estimated from five-dimensional training innovations. It mixes
unmodeled dynamics, model-fit error, and measurement-noise effects; it is not a
separate estimate of process covariance `Q_w` and measurement covariance `R_v`.

The independently sampled physical wind has per-axis standard deviation
`0.04 m/s` and is injected into the simulated plant. Its covariance is not
included in the current chance-tightening covariance. Therefore this run does
not establish wind-risk calibration, recursive feasibility, closed-loop
stability, or a hard guarantee for arbitrary winds and random seeds.

## Reproduce and verify

From the repository root:

```bash
python -u experiments/copyBI_pelican_probabilistic_boundary_advantage/generate_smpc_figures.py
python -u experiments/copyBI_pelican_probabilistic_boundary_advantage/generate_smpc_summary_gif.py
python experiments/copyBI_pelican_probabilistic_boundary_advantage/verify_smpc_summary_gif.py
python -m pytest tests/test_copybi_probabilistic_boundary_advantage.py -q
```

Machine-readable metrics and artifact checks are retained in:

- `results_smpc/copyBI_smpc_figure_summary.json`
- `results_smpc/copyBI_smpc_summary_gif_audit.json`
- `results_smpc/copyBI_smpc_figure_data.npz`

Historical deterministic-MPC comparison outputs remain under `results/` and
named diagnostic directories. They are not used to claim a controller-wide
probability advantage because the exploratory seed scan included non-optimal
QP terminations and was not a valid Monte Carlo experiment.
