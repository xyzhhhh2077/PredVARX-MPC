# copyBD: three-dimensional fixed-anchor Pelican SMPC trial

This independent copy tests whether the real-data-trained Pelican model can use
three prescribed controlled directions rather than the two directions in
copyBB/copyBC.

## Design

- Training data: Waterloo Pelican flights 1-36 at 100 Hz.
- Held-out identification validation: flights 37-54.
- Fixed anchor: standardized `pos_x_mm`, `pos_y_mm`, and `pos_z_mm`.
- Latent dimension: `ell=5`; three fixed task axes plus two CRTE-selected free
  directions.
- Controller: `d=10`, `N=18`, `Q=80 E_c E_c'`, `Ru=0.18 I`.
- Joint risk: `alpha=0.05`, allocated over both faces, three task axes, and all
  18 prediction steps.
- Comparator: deterministic MPC with the covariance tightening removed. Both
  controllers receive the same innovation sample at every physical step.

The controller receives the future reference at each lifted prediction step;
the original trial incorrectly repeated the current sample over the full
horizon. The revised reference is a centered, one-cycle, closed 3D curve. It
varies on all axes and has rank three after centering:

```text
standardized x range -0.2500 to 0.2500
standardized y range -0.5500 to 0.5500
standardized z range -0.3200 to 0.3200
```

The corresponding physical reference ranges, using training-only scaling, are
approximately `x=0.240-0.694 m`, `y=-0.281-0.729 m`, and `z=0.988-1.556 m`.

## Identification evidence

```text
dual error                  1.76e-15
anchor preservation error  0
spectral radius             0.99993990
controllability rank        5/5
held-out x/y/z RMSE         [0.007788, 0.007338, 0.005528]
persistence RMSE            [0.007746, 0.007284, 0.005549]
```

The one-step model is slightly worse than persistence on x/y and slightly
better on z. The near-one R2 values reflect the slowly varying position
signals; they are not evidence by themselves of strong dynamics.

## Three-dimensional closed-loop result

Three seeds (`7, 19, 31`) were run for 1200 physical steps each. All six
controller runs completed with zero QP fallback and zero hard-bound violation.

```text
seed 7  SMPC RMSE [0.0536, 0.1270, 0.1954]
seed 19 SMPC RMSE [0.0588, 0.1236, 0.2225]
seed 31 SMPC RMSE [0.0686, 0.1398, 0.2199]
```

For representative seed 7, the simulated SMPC trajectory has nonzero ranges
`[0.5712, 1.0763, 0.4019]` on x/y/z. Thus this is an actual 3D trajectory, not
a 2D curve with a decorative z coordinate.

## Honest boundary

The chance constraints are not under critical pressure in this first 3D trial.
SMPC and deterministic MPC are numerically indistinguishable in this
noncritical trajectory. Therefore the result
shows that the 3D fixed-anchor extension can be identified and executed, but it
does not establish a 3D covariance-tightening advantage.

With unchanged `Q/Ru`, no seed-7 physical step reaches an input bound. The
future-reference fix and trackable reference design reduce seed-7 RMSE from
`[0.1318, 0.4055, 0.4802]` to `[0.0536, 0.1270, 0.1954]`. The z-axis output
range is still only about 63% of the requested range, so this is improved but
not exact 3D tracking. Boundary-pressure evidence remains in copyBC rather than
being mixed into this 3D trajectory-quality trial.

This remains a noisy frozen-model-in-the-loop simulation trained on real flight
data. It is not a real Pelican closed-loop flight, hardware-in-the-loop test,
recursive-feasibility proof, stability proof, or calibrated probability
certificate.

## Reproduce

```bash
cd experiments/copyBD_pelican_3d_fixed_anchor_smpc
"E:/MATLABinhere/bin/matlab.exe" -batch "copyBD_pelican_3d_fixed_anchor"
python run_3d_boundary_trajectory.py --seeds 7 19 31 --steps 1200
python verify_results.py
python -m pytest ../../tests/test_copybd_pelican_3d.py -q
```

Main artifacts:

- `results/copyBD_pelican_3d_fixed_anchor_data.mat`
- `results/copyBD_pelican_3d_fixed_anchor_prediction.png`
- `results/copyBD_3d_boundary_trajectory.json`
- `results/copyBD_3d_boundary_trajectory.npz`
- `results/copyBD_3d_boundary_trajectory.png`
