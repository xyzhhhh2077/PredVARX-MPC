# copyBH — Pelican input-cap-only ablation

## Purpose

Test whether modestly relaxing the inherited standardized motor-command box from
`±6` to `±6.5` or `±7` changes x/y near-boundary tracking.

This is an **experiment-configuration ablation**, not a theory or controller
change. `copyBG` remains untouched. copyBH keeps the same frozen identified
model, `Sigma_eps`, physical-wind model (`sigma=0.04 m/s`), seeds, SMPC/MPC
objective, `Q=80`, `Ru=0.18`, `H=18`, `d=10`, chance tightening, and QP.

`±6` is an inherited experiment-design cap, not a Pelican hardware rating.
`±6.5` and `±7` are extrapolation tests beyond parts of the observed training
command support.

## Fair sweep

All three runs start at the center and use the same copyBH experimental
`xy-face-dwell` reference:

`center -> x=3.57 (hold 24 s) -> center -> y=3.57 (hold 24 s) -> center`

This avoids the cap-dependent initial equilibrium produced when `boundary-tour`
starts directly at `x=3.60`, and it tests sustained faces rather than a brief
near-hard crossing. The reference generator is experimental input only; it does
not alter the controller or model.

```bash
COPYBH_INPUT_CAP=6.0 python experiments/copyBH_pelican_u_cap_ablation/run_u_cap_ablation.py \
  --steps 12000 --reference xy-face-dwell --hard-bound 3.9124117280376636 \
  --pressures 3.57 3.57 0.0 \
  --output-dir experiments/copyBH_pelican_u_cap_ablation/results_xy_dwell_u6

COPYBH_INPUT_CAP=6.5 python experiments/copyBH_pelican_u_cap_ablation/run_u_cap_ablation.py \
  --steps 12000 --reference xy-face-dwell --hard-bound 3.9124117280376636 \
  --pressures 3.57 3.57 0.0 \
  --output-dir experiments/copyBH_pelican_u_cap_ablation/results_xy_dwell_u65

COPYBH_INPUT_CAP=7.0 python experiments/copyBH_pelican_u_cap_ablation/run_u_cap_ablation.py \
  --steps 12000 --reference xy-face-dwell --hard-bound 3.9124117280376636 \
  --pressures 3.57 3.57 0.0 \
  --output-dir experiments/copyBH_pelican_u_cap_ablation/results_xy_dwell_u7
```

Only `6.0`, `6.5`, and `7.0` are accepted by `COPYBH_INPUT_CAP`.

## Face-petal multi-direction stress (reference-only)

Six-face demo path: center → face → on-face diagonal slide → center for
`+x, -x, +y, -y, +z, -z`. Default pressures `(3.40, 3.40, 3.57)` sit above the
H=18 chance mean bounds and inside the source hard box so risk constraints can
activate. Controller / covariance / chance law unchanged.

```bash
COPYBH_INPUT_CAP=7.0 python experiments/copyBH_pelican_u_cap_ablation/run_u_cap_ablation.py \
  --reference face-petal-stress \
  --output-dir experiments/copyBH_pelican_u_cap_ablation/results_face_petal
```

Default length is 18000 samples (180 s). Compare:
`results_face_petal/copyBH_face_petal_stress_comparison.png` and `REPORT.md`.

## Interpretation boundary

A larger cap can remove an input-reachability limitation. It does not alter the
frozen horizon-18 chance mean bounds. Therefore, a larger cap alone is not
expected to make x/y sustain `|s| >= hard_bound - 0.4` when chance tightening is
the active constraint.
