# copyBI ±3.6 reference, ±7 input-cap diagnostic

## Setup

- Reference pressures: x/y/z = ±3.6 standardized
- Input-command cap: ±7 standardized
- Requested length: 18,000 plant steps
- Hard state bound: ±3.912411728 standardized
- Wind sigma: 0.04 m/s
- Wind seed: 7016
- Innovation seed: 16

## Result

| controller | completed steps | QP failure step | violations before failure | minimum hard margin | max input | fallback |
|---|---:|---:|---:|---:|---:|---:|
| SMPC | 9,690 | 9,690 | 0 | +0.039121 | 7.0 | 0 |
| hard MPC | 1,680 | 1,680 | 0 | +0.174825 | 7.0 | 0 |

Peak absolute task coordinates before failure:

- SMPC: [3.843838, 3.873290, 1.943978]
- hard MPC: [3.737587, 1.332134, 0.502589]

## Interpretation

The uniform ±3.6 face pressure improves neither controller enough to finish the 18,000-step face-petal trajectory. The likely feasibility bottleneck is the compound face-plus-slide demand, not merely the axial pressure or input cap. This is a failed feasibility diagnostic, so no sampled-advantage figure is generated.
