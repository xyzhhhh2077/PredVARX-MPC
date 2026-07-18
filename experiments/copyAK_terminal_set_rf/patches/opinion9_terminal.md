# Opinion 9 terminal / RF probe (copyAK)

## Status
MATLAB real run completed (parent completion after subagent timeout).

## Probe results
- rho(Ahat)=0.942272, Schur=1
- dlyap residual F=4.743e-14, V-decrease OK=1.000
- RF chain: C1–C3 PROBED-PASS; C4–C6 UNPROVED
- alpha_term calibrated ≈ 448

## Closed-loop (T_cl=400)
Medium y_max=0.70 and tight 0.50: OFF vs COST vs SET show essentially identical cover/MAE/qp/soft.
Terminal cost/set act as soft regularizers here, not RF certificates.

## Non-claims
No chance-constraint recursive feasibility theorem; no closed-loop stability theorem;
spectral box ≠ exact ellipsoid Xf.
