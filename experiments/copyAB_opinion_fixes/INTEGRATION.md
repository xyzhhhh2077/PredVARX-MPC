# copyAB_opinion_fixes — integrated result

## Unified MATLAB test summary (parent re-run)

```text
SUMMARY fails=0 of 11
PASS test_copyAB_core_geometry
PASS test_opinion01_geometry
PASS test_opinion02_noise_opt
PASS test_opinion03_M_vs_tau
PASS test_opinion04_ivr_trace
PASS test_opinion05_sigma_denoms
PASS test_opinion06_cross_cov
PASS test_opinion07_sigma_obs
PASS test_opinion08_fallback_flag
PASS test_opinion09_terminal_default_off
PASS test_opinion10_input_residualize
```

## Per-opinion implementation

| # | Fix | Default | Status |
|---|---|---|---|
| 1 | Geometry assert left/right/column/dual | enforce ON | done |
| 2 | Declared Σ_n free-noise flags + isotropic collapse | always | done |
| 3 | Comments: M vs τ (latent energy) | n/a | done |
| 4 | Always export ivr_iter/trace/delta | always | done |
| 5 | Sigma_eps multi-denom; default return T-2 | T-2 default | done |
| 6 | cross_cov_diagnostics.m offline | offline print | done |
| 7 | sigma_obs_support_diag.m online last diag | diagnostic only | done |
| 8 | fallback cert `qp` / `uncertified_fallback` | always | done |
| 9 | optional terminal cost on z_N | OFF | done |
| 10 | optional input residualize IVR | OFF | done |

## Hard non-claims preserved

- Reconstruction geometry ≠ closed-loop tracking zero error
- Declared Σ_n free dual ≠ global unknown-noise optimum
- Drop-cross Σ_y still used in SMPC (diag only for opinion 6)
- Uncertified fallback ≠ chance guarantee
- Strict convexity ≠ recursive feasibility/stability
- Optional residualize-on-U ≠ proved input-conditional PredVARX
