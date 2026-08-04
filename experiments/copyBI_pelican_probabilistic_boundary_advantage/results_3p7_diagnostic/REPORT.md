# copyBI ±3.7 diagnostic run

## Requested reference

- x: ±3.7 standardized
- y: ±3.7 standardized
- z: ±3.7 standardized
- hard bound: ±3.912411728
- input cap: ±7
- steps requested: 18,000
- wind seed: 7016
- innovation seed: 16

## Result

| controller | completed steps | QP failure | status | hard violations before failure | minimum hard margin | saturated plant steps |
|---|---:|---:|---|---:|---:|---:|
| SMPC | 9,690 | 9,690 | infeasible | 0 | +0.024425 | 6,250 |
| hard MPC | 1,680 | 1,680 | infeasible | 0 | +0.160685 | 1,370 |

The fixed six-face trajectory at ±3.7 cannot be completed by either controller
with the current frozen model, face-petal slide geometry, input cap ±7, and
shared disturbance realization. No fallback was used. This run therefore does
not provide a valid full-length SMPC-vs-MPC hard-boundary comparison.

At the SMPC failure step, the reference was approximately
`[-0.6233, -3.7000, 1.8000]`; the realized task was approximately
`[-0.4701, -3.8880, 1.9374]`.

At the hard-MPC failure step, the reference was approximately
`[3.7000, -0.7203, 0.4500]`; the realized task was approximately
`[3.7517, -0.8342, 0.4957]`.
