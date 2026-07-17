# Opinion 4 patch report

status: done

## Goal

Export IVR iteration diagnostics so that **stop ≠ proved convergence**:
callers can inspect the numerical stationarity trajectory instead of treating
`break` as a monotone/global optimum certificate.

## Math / semantics

Free-direction IVR iterates a predictable-energy quantity

$$
\tau_k = \frac{1}{N}\operatorname{tr}(X_{\mathrm{pred}} X_{\mathrm{pred}}^\top)
$$

and a projector change

$$
\delta_k = \|V_{\mathrm{new}}V_{\mathrm{new}}^\top - V V^\top\|_F.
$$

**Stop rule** (existing code):

$$
|\tau_k - \tau_{k-1}| < \mathrm{tol}\cdot\max(|\tau_k|,1)
\quad\text{and}\quad
\delta_k < \mathrm{tol}.
$$

This is **numerical stationarity**, not a proof that \(\tau_k\) is monotone or
that the free subspace is globally optimal. Export the full \(\tau\) history so
the claim boundary is auditable.

## Code change

File: `split_control_free_ivr_varx.m` only (copyAB).

1. Comment at the `break` site:
   `stop = numerical stationarity (trace + subspace delta), not a proof of monotone convergence`.
2. Always set (including `r = 0` / tracked-only):
   - `stats.ivr_iter` — iterations executed (`0` if no free directions)
   - `stats.ivr_trace` — per-iter \(\tau\) history (empty if `r=0`)
   - `stats.ivr_subspace_delta` — final projector Frobenius change (`0` if `r=0`)

No change to the stop rule itself, to IVR updates, or to the free dual /
VARX fit.

## Test

`tests/test_opinion04_ivr_trace.m`:

1. Free case (`ell > q`): all three fields exist; `ivr_iter >= 1`;
   `numel(ivr_trace) == ivr_iter`; finite nonneg. `ivr_subspace_delta`.
2. Tracked-only (`ell = q`): fields still exist with
   `ivr_iter=0`, empty `ivr_trace`, `ivr_subspace_delta=0`.

```matlab
cd('.../experiments/copyAB_opinion_fixes');
addpath(pwd); addpath(fullfile(pwd,'tests'));
test_opinion04_ivr_trace;
```

## Boundary

- Exporting the trace does **not** prove monotone convergence.
- Stationary stop does **not** certify global optimality of the free span.
- Output-history IVR remains not a proved input-conditional PredVARX subspace
  (see Opinion 10 / `input_residualize` candidate flag).
