# Vector double-projection pilot (qkag125 context)

Run from this directory:

```sh
octave --quiet --eval "addpath(pwd); run_pilot(20)"
```

`results/paired_metrics.csv` contains three methods for each of 20 paired seeds. The Octave `rand` and `randn` states are reset to the same seed at the start of each replicate; identification and test trajectories are reused across methods. Both learned candidates use the same `P,R,Pbar,Rbar` reconstruction and the same centered one-step VARX refit. All output metrics are held-out and do not select hyperparameters.

The baseline calls this repository's `experiments/copyR_moqin_oblique/predvarx_identify_moqin.m`. It is a Mo--Qin Algorithm-1-inspired implementation with a separate VARX extension, not a verified byte-for-byte port of the authors' `PredVAR.m`. The two candidate refinements are *vector analogies*, not the tensor-CP algorithm of Chang, Huang, Yao and Yu, arXiv:2606.08560. The second step residualizes the lagged target against contemporaneous other estimated factors, then applies coordinate hard-thresholding before QR and reconstructing the dual bases. The threshold is fixed a priori for this pilot, not tuned to test data.

This pilot evaluates loading subspace recovery and one-step observation prediction. It has no MPC, chance constraints or closed-loop control, and is not evidence of control advantage. The 20-seed table belongs to an exploratory experiment only. See `PredVAR_double_projection_pilot.tex` in the linked Overleaf project for results and limits.
