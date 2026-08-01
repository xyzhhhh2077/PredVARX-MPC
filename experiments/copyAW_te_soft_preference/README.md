# copyAW: copyAU soft preference on Tennessee Eastman data

This experiment applies the copyAU construction to native high-dimensional
Tennessee Eastman Process (TEP) data. It uses 41 measured variables, 11
manipulated variables, and a five-dimensional latent model. The 41-dimensional
output is not created by lag stacking.

## Data provenance

- Dataset: Rieth et al., *Additional Tennessee Eastman Process Simulation Data*
- DOI: <https://doi.org/10.7910/DVN/6C3JR1>
- Dataverse source file ID: `3031241`
- Source file: `TEP_FaultFree_Training.RData`
- Source MD5: `ec126484534331f85001d8c4ebce6d17`
- Repository snapshot: complete simulation runs 1 through 5 only

The committed MAT snapshot contains 2,500 samples: `41 x 2500` measurements
and `11 x 2500` inputs. Runs 1-3 (1,500 samples) train the model; runs 4-5
(1,000 samples) are held out. Run boundaries are excluded from all lagged
regressions and prediction metrics.

## Relation to copyAU

The experiment retains copyAU's normalized mixture of soft output preference
and finite-horizon input authority, followed by the anchored CRTE free block
with fixed `mu=0.10`, `Ntr` epsilon `1e-5`, `q=2`, and `ell=5`.

Unlike copyAU's synthetic column-decay weights, copyAW gives elevated soft
weights to TEP measurements 7-9, 11-13, 15-16, and 18. These are reactor,
separator, and stripper operating variables. Every other measured variable
retains a nonzero weight, so this remains a soft preference rather than a hard
axis lock.

## Claim boundary

This repository experiment performs offline one-step prediction on two unseen
normal simulation runs. It tests whether the copyAU identification construction
can be exercised on real public `41 x 11` benchmark trajectories. It does not
run COSTEP, replace the plant controller, or establish closed-loop tracking,
constraint satisfaction, stability, or recursive feasibility.

Run in MATLAB R2024a:

```matlab
run('experiments/copyAW_te_soft_preference/copyAW_te_soft_preference.m')
```

Results are written to `results/`.

## R2024a run result

- Dual-basis error: `2.96e-15`
- Identified spectral radius: `0.9080`
- Held-out transitions: `998`
- Learned-task RMSE: `[0.0366, 0.0827]`
- Persistence RMSE: `[0.0531, 0.1142]`
- Learned-task R2: `[-0.0050, 0.0608]`
- Full standardized-output RMSE: `0.8531`

The two learned-task RMSE values are lower than the persistence baseline, but
the R2 values are weak and the plotted predictions are visibly smoother than
the measurements. This is evidence that the copyAU construction runs on the
native high-dimensional benchmark, not evidence of an accurate TEP model or
an effective closed-loop controller.