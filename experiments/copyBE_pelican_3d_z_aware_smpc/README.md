# copyBE: Pelican 3D z-aware fixed-anchor SMPC

## Purpose

copyBE keeps the exact fixed task anchor

\[
E_c=[e_x,e_y,e_z]
\]

for standardized Pelican position outputs. It changes the latent dynamics from the
copyBD first-order model to a second-order VARX model:

\[
z_{k+1}=A_1z_k+A_2z_{k-1}+B(u_k-\bar u)+\varepsilon_{k+1}.
\]

The augmented MPC state is `[z_k; z_{k-1}]`. This integrates the missing
position-rate memory without changing the fixed xyz task geometry or adding a
new training set.

## What is fixed and what is diagnosed

- `P R' E_c = E_c` remains exact, so x, y, and z are all protected anchors.
- All four recorded motor commands remain independent inputs.
- The normalized collective motor direction `[1,1,1,1]'/2` is evaluated as a
  diagnostic. It is not an extra input and does not inject new information.
- `Q=80 E_c E_c'`, `R_u=0.18 I`, `d=10`, and `N=18` are unchanged from copyBD.
- copyBD is retained as the first-order comparison; copyBE is a new experiment.

## Identification evidence

Held-out one-step task RMSE in standardized coordinates is:

| axis | copyBE RMSE | persistence RMSE |
|---|---:|---:|
| x | 0.0001440 | 0.0077466 |
| y | 0.0001315 | 0.0072840 |
| z | 0.0000754 | 0.0055495 |

The augmented spectral radius is `0.9996568`. The finite-horizon all-input gain
for `[x,y,z]` is `[0.1307, 0.1056, 0.5173]`. The collective-mode diagnostic is
primarily vertical: its horizon gain is `[0.00483, 0.00767, 0.21730]`.

## Held-out recorded-input rollout audit

Each rollout uses a true held-out output only for initialization. Future outputs
are predicted recursively using recorded motor commands; no true future output
is fed back.

| horizon | copyBD z RMSE | copyBE z RMSE | result |
|---:|---:|---:|---|
| 10 | 0.05507 | 0.04142 | 24.8% lower |
| 50 | 0.26260 | 0.21778 | 17.1% lower |
| 100 | 0.47545 | 0.47614 | effectively tied (+0.15%) |

copyBE does not dominate copyBD on every axis: x/y recursive rollout errors are
higher, especially at 100 steps. The defensible claim is therefore a short- and
medium-horizon z improvement, not a universally better predictor.

## Frozen-model-in-the-loop closed loop

The reference, noise generation, controller weights, constraints, and three
seeds match copyBD. Results are standardized-coordinate RMSE:

| seed | x | y | z | violations | fallback |
|---:|---:|---:|---:|---:|---:|
| 7 | 0.0480 | 0.0452 | 0.0128 | 0 | 0 |
| 19 | 0.0532 | 0.0559 | 0.0112 | 0 | 0 |
| 31 | 0.0829 | 0.0558 | 0.0135 | 0 | 0 |

For seed 7, z peak-to-peak recovery changes from copyBD's `62.8%` to copyBE's
`102.3%`. No motor command reaches the standardized `[-6,6]` bounds. SMPC and
nominal MPC coincide because this centered trajectory does not activate the
chance-constraint tightening; boundary-pressure evidence remains in copyBC.

## Scope limitation

This is a frozen identified-model simulation, not a Pelican flight replay under
closed-loop intervention and not hardware validation. The strong closed-loop z
result shows that the second-order identified model supports vertical control
inside its own dynamics. The held-out rollout audit supports improved z dynamics
through 50 samples, but the 100-step tie and x/y degradation prevent a broader
claim that copyBE is globally superior to copyBD.

## Reproduction

```bash
"E:/MATLABinhere/bin/matlab.exe" -batch "cd('C:/Users/ROG/predvarx-repo/experiments/copyBE_pelican_3d_z_aware_smpc'); copyBE_pelican_3d_z_aware"
python experiments/copyBE_pelican_3d_z_aware_smpc/run_3d_z_aware_trajectory.py --seeds 7 19 31 --steps 1200
python experiments/copyBE_pelican_3d_z_aware_smpc/verify_results.py
python -m pytest tests/test_copybe_pelican_3d_z_aware.py tests/test_copybd_pelican_3d.py -q
```
