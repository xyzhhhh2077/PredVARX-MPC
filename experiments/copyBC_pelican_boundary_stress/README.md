# copyBC Pelican boundary-stress test

This experiment tests whether the fixed-anchor SMPC chance constraints affect
closed-loop behavior near the position boundary. It does not replace the
normal held-out-flight experiment in `copyBB`.

## Fair comparison

Both controllers use the same frozen Pelican model, fixed position anchor,
initial state, critical reference, input bounds, `Q = 80 E_c E_c'`,
`Ru = 0.18 I`, prediction horizon, and disturbance sample at every plant step.

The only controller difference is:

- SMPC: `alpha_joint = 0.05`, with propagated covariance tightening;
- deterministic MPC: zero covariance tightening, with the same hard bound.

Two independent scenarios place either the positive-x or positive-y reference
0.005 standardized units inside the corresponding hard bound. This is about
4.54 mm on x and 4.59 mm on y, or less than one one-step identified-model
position-noise standard deviation. Each comparison starts at its critical task
position.

## Result

For seeds 7, 19, and 31, each run contains 1200 plant steps and 120 QP solves.

| Scenario | Controller | Violations by seed | Total violations | QP fallback |
|---|---|---:|---:|---:|
| x critical | Fixed-anchor SMPC | 0, 0, 0 | 0 / 3600 | 0 |
| x critical | Deterministic MPC | 423, 226, 461 | 1110 / 3600 | 0 |
| y critical | Fixed-anchor SMPC | 0, 0, 0 | 0 / 3600 | 0 |
| y critical | Deterministic MPC | 159, 174, 131 | 464 / 3600 | 0 |

The SMPC chance constraint was active in all 120 QP solves for every seed and
both critical axes. The deterministic MPC crossed the hard boundary by up to
0.0532 standardized units on x and 0.0381 on y, whereas the SMPC remained
inside it.

## Moving boundary trajectory

`copyBC_moving_boundary.png` adds a moving-reference test without replacing
the fixed-point evidence. Its reference keeps x at 0.005 standardized units
inside the positive hard bound while y follows a sinusoid of amplitude 1.0
and period 600 plant steps. The 1200-step run therefore contains two complete
cycles along the boundary.

With the same disturbances, the SMPC violation counts are 0, 0, and 0 for
seeds 7, 19, and 31. Deterministic MPC produces 424, 228, and 456 violating
steps. The figure includes x(t), y(t), and the two-dimensional x-y path.

## Evidence boundary

This is a deliberately constructed identified-model-in-the-loop stress test.
It demonstrates that covariance tightening changes the control action and
prevents boundary crossings under the tested model and disturbances. It is
not real-flight or hardware-in-the-loop evidence, and `0 / 3600` is not a
general probability, recursive-feasibility, or stability certificate.

Run:

```bash
python experiments/copyBC_pelican_boundary_stress/run_boundary_stress.py
```

Outputs are saved under `results/` as JSON, NPZ, and PNG.
