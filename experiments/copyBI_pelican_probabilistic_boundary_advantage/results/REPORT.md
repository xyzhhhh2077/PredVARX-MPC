# copyBI — sampled hard-boundary advantage

## What this run verifies

- Same innovation and physical-wind realization: `7671c9ab8a456370278ed9052267a126890699a23d8731eeb42e38fe263b3375`.
- Both controllers completed **18,000** plant steps; fallback = 0.
- SMPC: **0** hard-violation steps, minimum margin **0.019400**, active chance QPs **375**.
- Hard-constrained nominal MPC: **6** hard-violation steps, minimum margin **-0.002570**.

## Claim boundary

This seed pair was deliberately selected after an exploratory scan of 10 consecutive pairs. Six deterministic-MPC probes ended at a non-optimal QP, so that scan is not a valid violation-frequency estimate. The retained run is one full-length sampled crossing event, not a Monte Carlo estimate of the violation probability. The SMPC prediction covariance is the identified `Sigma_eps_aug`; physical wind is injected into both plants but is not added to that chance covariance. Therefore this run demonstrates the mechanism and a same-noise contrast, not probability calibration, recursive feasibility, or a robust hard guarantee under physical wind.
