# copyBB: fixed-anchor CRTE on Pelican data with model-in-the-loop SMPC

This copy separates real-data identification from simulated control validation.

- **Training:** Waterloo Pelican flights 1-36, measured at 100 Hz.
- **Offline validation:** held-out flights 37-54.
- **Controlled directions:** prescribed standardized x/y position axes; data do
  not learn or rotate these axes.
- **Free latent complement:** copyAU/copyBA CRTE spectral construction with
  fixed `mu=0.10` and paper trace normalization.
- **Control validation:** the frozen identified stochastic model is the plant.
  Held-out flights 37-41 supply only initial states and reference trajectories.
  Their measured next states are never injected into the simulated plant.

## Relation to the article

This copy follows the article's key controlled-direction premise more closely
than copyBA: `E_c` is specified before identification and is preserved exactly.
Only the free complement is selected from data. The multirate noisy SMPC is an
engineering validation layer, not part of the original identification formulas.

The controller uses the copyAU objective and input weight:

```text
Q = 80 E_c E_c'
Ru = 0.18 I
```

The 100 Hz identified plant is lifted to a 10 Hz decision model with
`d=10`, `N=18`. The decision-step disturbance covariance is computed from the
physical-step model:

```text
Sigma_d = sum(i=0..d-1) A^i Sigma_eps (A^i)'
```

The SMPC and open-loop comparator receive the same noise sample at every
physical step.

## Identification result

```text
dual error                 1.46e-15
anchor preservation error  0
spectral radius             1.00005164
held-out x/y RMSE           [0.007789, 0.007338]
persistence RMSE            [0.007746, 0.007284]
```

The near-unit spectrum is consistent with slowly varying position, but the
one-step model is about 0.6-0.8% worse than persistence. The high R2 values
therefore do not establish strong one-step dynamics.

## Noisy closed-loop result

Five held-out reference flights were simulated for 1200 physical steps (12 s),
using seeds 7, 19, and 31. Across all 15 runs:

```text
x RMSE range                0.073-0.163 training standard deviations
y RMSE range                0.084-0.179 training standard deviations
mean improvement vs open    83.1% / 83.1% for x / y
minimum improvement         73.5%
QP fallback count           0
chance-bound violations     0
input saturation            0-17.5% of physical steps
```

For the formal seed-7 run, reference-normalized RMSE is 0.22-0.72. Thus the
controller tracks all five references substantially better than the same noisy
plant in open loop, but the smallest-amplitude x references remain relatively
hard.

## Honest boundary

This is evidence that the fixed-anchor controller works on a stochastic model
identified from real Pelican flights. It is not hardware-in-the-loop or a real
Pelican closed-loop experiment. The simulated controller also uses larger input
variation than the recorded expert on several channels (roughly 2-8x in the
formal run), and saturation reaches 17.5% for one seed. It should not be claimed
that the control action reproduces the pilot or that recursive feasibility and
stability are proved.

## Reproduce

```bash
cd experiments/copyBB_pelican_fixed_anchor_smpc
"E:/MATLABinhere/bin/matlab.exe" -batch "copyBB_pelican_fixed_anchor"
python run_fixed_anchor_closed_loop.py
python verify_results.py
python -m pytest ../../tests/test_copybb_pelican_fixed_anchor.py -q
```

Main artifacts:

- `results/copyBB_pelican_fixed_anchor_data.mat`
- `results/copyBB_pelican_fixed_anchor_metrics.txt`
- `results/copyBB_pelican_fixed_anchor_prediction.png`
- `results/copyBB_fixed_anchor_closed_loop.json`
- `results/copyBB_fixed_anchor_closed_loop.npz`
- `results/copyBB_fixed_anchor_closed_loop.png`
