# copyBH face-petal-stress REPORT

Reference-only multi-face stress test on frozen copyBG controller.

## Setup
- reference: face-petal-stress (center -> face -> on-face slide -> center x6)
- pressures: [3.4, 3.4, 3.57]
- hard: 3.912412
- terminal chance mean bound: [2.935428335885822, 3.148504892150341, 3.583293760582082]
- input cap: [-7.0, 7.0]
- steps: 18000 @ 100 Hz, wind sigma=0.04 m/s
- seeds: wind=7007, innovation=7

## Results (SMPC vs deterministic MPC, shared noise/wind)

| metric | SMPC | det. MPC |
|---|---:|---:|
| completed_steps | 18000 | 18000 |
| hard_violation_steps | 0 | 0 |
| hard_violation_rate | 0.0 | 0.0 |
| active_qp_steps | 303 | 0 |
| qp_failure_count | 0 | 0 |
| fallback_count | 0 | 0 |
| input_saturation_steps | 6630 | 6610 |
| min_hard_margin | 0.09553828968633837 | 0.04873073546126605 |
| peak |x|,|y|,|z| | [3.817, 3.765, 3.736] | [3.864, 3.765, 3.736] |
| RMSE xyz | [0.3445, 0.2526, 0.0811] | [0.3478, 0.252, 0.0809] |

## Readout
- Both controllers finish 18000 steps with **0 hard violations**.
- SMPC chance QP active on **303** decision steps; MPC active=0 (no chance layer).
- SMPC keeps a larger minimum hard margin (**0.0955** vs **0.0487**), mainly on x.
- Peak |x|: SMPC **3.817** < MPC **3.864** (SMPC pulls back under risk constraints).
- This is a **directional stress demo**, not a sustained near-hard cruise (that remains v4.0 z-face).

## Artifacts
- `copyBH_face_petal_stress_comparison.png`
- `copyBH_pelican_u_cap_ablation_smpc.png`
- `copyBH_pelican_u_cap_ablation_deterministic_mpc.png`
- `copyBH_pelican_u_cap_ablation.{json,npz,mat}`
- `preview_reference.png`
