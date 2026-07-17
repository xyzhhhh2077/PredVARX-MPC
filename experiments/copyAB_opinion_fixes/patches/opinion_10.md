# Opinion 10 patch report

status: done

## Goal

Add an **optional** input-residualization IVR mode for free-direction selection, **default OFF**, so the default path stays output-history IVR and all existing geometry assertions still pass.

## Design

Default (`input_residualize=false`):

- Free IVR runs on `Yp = Nperp'*yc` lags only (no U).
- Free span = predictable-in-output-history; U enters only in the coupled VARX step.

Optional candidate scheme 1 (`'input_residualize', true`):

1. Project to free coordinates: `Yp = Nperp'*yc`.
2. Ridge-regress `Yp` on centered `U`:
   `Bu = (U U' + λ I) \ (U Yp')`, `Yp_ivr = Yp - Bu' U`.
3. Run the existing iterative IVR on residual series `Yp_ivr`.
4. Free dual / tracked split / coupled VARX refit unchanged after free loading `V` is chosen.

### Explicit non-claim

`stats.ivr_input_conditional = true` means only that this **candidate** residualize-on-U path was used. It does **not** claim a proved optimal input-conditional PredVARX free subspace.

## Code change

File: `split_control_free_ivr_varx.m`

- Signature: `...(y,u,ell,tracked,Sigma_n,varargin)`
- Name-value: `'input_residualize'` (logical, default `false`)
- Stats field: `stats.ivr_input_conditional` (always present, false by default)

Does **not** modify `main/`, `copyAA`, or `copyX`.

## Tests

`tests/test_opinion10_input_residualize.m`:

1. Default / explicit-false: geometry (dual, left/right coverage, tracked columns) PASS; `ivr_input_conditional==false`.
2. Residualize mode: runs without crash; same geometry PASS; `ivr_input_conditional==true`; finite IVR trace.

Also re-run `test_copyAB_core_geometry` to confirm default path regression.

```matlab
cd('.../experiments/copyAB_opinion_fixes');
addpath(pwd); addpath(fullfile(pwd,'tests'));
test_copyAB_core_geometry;
test_opinion10_input_residualize;
```

## Boundary

Residualize-on-U is an engineering candidate for “input-aware free directions,” not a theorem. Geometry guarantees remain reconstruction-layer only.
