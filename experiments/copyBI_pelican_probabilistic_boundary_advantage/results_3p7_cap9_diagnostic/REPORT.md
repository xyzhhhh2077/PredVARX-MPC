# copyBI ±3.7 reference, ±9 input-cap diagnostic

## Fixed setup

- Reference pressures: x/y/z = ±3.7 standardized
- Input-command cap: ±9 standardized
- Requested length: 18,000 plant steps
- Hard state bound: ±3.912411728 standardized
- Wind sigma: 0.04 m/s
- Wind seed: 7016
- Innovation seed: 16
- Same controller implementation and disturbance realization as copyBI

## Result

| controller | completed steps | QP failure step | hard violations before failure | minimum hard margin | max input infinity norm | fallback |
|---|---:|---:|---:|---:|---:|---:|
| SMPC | 4,450 | 4,450 | 0 | +0.072487 | 9.0 | 0 |
| hard-constrained nominal MPC | 1,660 | 1,660 | 0 | +0.029033 | 9.0 | 0 |

Peak absolute task coordinates before failure:

- SMPC: [3.839925, 1.024242, 1.021112]
- hard MPC: [3.883379, 0.890832, 0.541917]

## Interpretation boundary

The ±9 input cap does not make the ±3.7 face-petal trajectory complete. Both controllers saturate at the expanded cap and later encounter a non-optimal QP before completing 18,000 steps. Because neither trajectory completes, this run is a feasibility diagnostic and not a valid sampled SMPC-versus-MPC advantage comparison. No success figure is generated.
