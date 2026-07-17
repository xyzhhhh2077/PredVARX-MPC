# Opinion 2 patch report

status: done

## Opinion (pure-math review)

The free dual is optimal only for a **declared** sensor-noise covariance
`Sigma_n` with fixed free loading `V`. Reporting a free-block noise
objective without declaring this scope is misleading. Also, the isotropic
case `Sigma_n ∝ I` is a known degeneration: the noise-optimal free dual
must collapse to the orthogonal special case `R = P`.

## Design logic

| Object | Role |
|---|---|
| `trace(Rf' Sigma_n Rf)` s.t. `Rf' Pf = I` | free-block noise objective |
| `free_noise_baseline` | same objective with orthogonal free dual `Rf = Pf` |
| `free_noise_objective` | optimized free dual value |
| `free_noise_improvement` | baseline − objective (≥ 0 by construction) |
| `free_noise_is_declared_Sigma_n` | **true**: metric uses argument `Sigma_n`, not residual cov |
| `isotropic_collapsed_to_R_eq_P` | **true** iff `‖R−P‖_F < 1e-8` (isotropic degeneration) |

No change to the free dual formula itself (already correct in split_*).

## Files changed (only under `experiments/copyAB_opinion_fixes/`)

1. `split_control_free_ivr_varx.m`
   - Confirmed existing: `free_noise_baseline`, `free_noise_objective`,
     `free_noise_improvement`.
   - Added: `stats.free_noise_is_declared_Sigma_n = true`.
   - Added: `stats.isotropic_collapsed_to_R_eq_P = (norm(R-P,'fro') < 1e-8)`.
2. `tests/test_opinion02_noise_opt.m` (new)
   - Isotropic: `‖R−P‖_F < 1e-8` and `isotropic_collapsed_to_R_eq_P == true`.
   - Heteroscedastic: `free_noise_improvement > 0` and collapse flag false.
   - Both cases: `free_noise_is_declared_Sigma_n == true`.
3. `copyAB_opinion_fixes.m`
   - Metrics text now logs the two new flags.
4. `patches/opinion_2.md` (this file).

## Hard boundary (honest)

- Free dual optimality is **only** w.r.t. declared `Sigma_n` and fixed free
  loading `V` from IVR. It is not an estimated residual-noise optimum.
- Collapse flag is a numerical check (`‖R−P‖_F < 1e-8`), not a symbolic
  proof that `Sigma_n` is isotropic; it is the observable degeneration test.

## Verification

```matlab
cd('.../experiments/copyAB_opinion_fixes');
addpath(pwd); addpath(fullfile(pwd,'tests'));
test_opinion02_noise_opt;
test_copyAB_core_geometry;
```

Real MATLAB stdout (2026-07-18):

```text
PASS opinion02: iso ||R-P||=8.335e-16 collapsed=1; het improvement=1.563e-01 collapsed=0; declared_Sigma_n=1/1
PASS split extractor: iso ||R-P||=9.065e-16; het free-oblique=3.020e+00; noise objective 0.429293 -> 0.146887 (improvement 2.824e-01)
```

Gates:
- isotropic: `‖R−P‖_F < 1e-8` ✓ (8.3e-16)
- isotropic: `isotropic_collapsed_to_R_eq_P == true` ✓
- heteroscedastic: `free_noise_improvement > 0` ✓ (0.156)
- both: `free_noise_is_declared_Sigma_n == true` ✓
