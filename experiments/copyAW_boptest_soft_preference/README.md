# copyAW: copyAU soft preference on BOPTEST HVAC data

This experiment replaces the earlier Tennessee Eastman offline copyAW. It
applies the copyAU soft-preference CRTE construction to the official BOPTEST
`multizone_office_complex_air` testcase.

## Real data contract

- Source: `https://api.boptest.net`
- Testcase: `multizone_office_complex_air`
- Outputs: 15 zone temperatures plus 15 zone CO2 measurements (`30 x T`)
- Inputs: supply-air temperature and duct-pressure setpoints on three floors
  (`6 x T`)
- Communication step: 900 seconds
- Excitation: deterministic bounded two-level PRBS base held for four samples,
  with small samplewise dither to separate the six input channels
- Pilot segments 1-2 train the model; segment 3 is held out

Each segment is reinitialized at a different seasonal start time with zero
warmup because the public API cannot complete the seven-day warmup for this
large testcase within its request timeout. Segment boundaries are excluded
from all lagged regressions. The collector saves the exact API testcase, UTC
collection time, control policy, input bounds, output metadata, and transient
BOPTEST test IDs.

Each saved record pairs the control sent to `/advance` with the measurement
returned by that same call. Consequently the identified transition uses
`y_k`, `u_{k+1}` to predict `y_{k+1}`; the MATLAB contract test locks this API
alignment and excludes cross-segment transitions.

This first Git version contains 24 samples per segment (46 valid training
transitions and 23 validation transitions). A fourth instance was attempted
repeatedly, but the public service stopped allocating new instances while its
health endpoints remained available. The three completed instances are kept
as a clearly labeled pilot, not presented as a final statistical experiment.

## Relation to copyAU

The experiment retains copyAU's normalized soft output-preference and
finite-horizon input-authority objective. It learns two task directions, then
anchors them while learning three CRTE free directions. Parameters remain
`ell=5`, `mu=0.10`, reach horizon 18, preference strength 0.78, and paper trace
normalization epsilon `1e-5`.

All 30 physical outputs have nonzero preference. The 15 temperature channels
use weight 1.0 and the 15 CO2 channels use weight 0.85. These are soft weights,
not fixed task axes and not a claim that temperature is globally more important.

## Run

```bash
python experiments/copyAW_boptest_soft_preference/collect_boptest_data.py
```

Then in MATLAB R2024a:

```matlab
run('experiments/copyAW_boptest_soft_preference/copyAW_boptest_soft_preference.m')
```

## Claim boundary

The run uses active control overrides on a public building simulator, then
performs held-out offline one-step prediction. It does not yet run MPC in the
loop. Weather and occupancy are present in the simulator but are not included
as explicit disturbance regressors in this first version. Therefore this copy
does not establish closed-loop tracking, constraint satisfaction, stability,
recursive feasibility, or a complete HVAC control model.

## Pilot result

The real MATLAB run completed with dual error `4.64e-15` and identified-model
spectral radius `0.9994`. On the 23 held-out transitions, task RMSE was
`[0.01421, 0.03840]`, while persistence RMSE was lower at
`[0.005448, 0.01891]`; task R2 was `[-1.136, 0.3257]`. The pilot therefore
demonstrates that the BOPTEST-to-copyAU pipeline runs, but it does not show a
prediction advantage. The short training set and omitted weather/occupancy
regressors are unresolved limitations, not evidence that tuning has succeeded.