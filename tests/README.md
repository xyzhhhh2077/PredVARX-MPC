# PredVARX-MPC unit tests

Canonical root suite using MathWorks `matlab.unittest` (class-based).

## Layout

| File | Covers | Source under test |
|------|--------|-------------------|
| `tPredvarxIdentifyOblique.m` | dual-basis identities `R'P=I` etc. | `experiments/copyO_oblique/predvarx_identify_oblique.m` |
| `tCenteredSmpcStep.m` | centered prediction + Boole tightening | `experiments/copyP_centered_smpc/centered_smpc_step.m` |
| `tControlAwareSubspace.m` | tracked-axis coverage | `experiments/copyQ_control_aware/control_aware_subspace_varx.m` |
| `localpaths.m` | path fixture | repo roots |
| `run_all_tests.m` | suite entry | all of the above |

Old script-style root `test_*.m` files were replaced by these class-based tests.

## Run

From repo root in MATLAB:

```matlab
cd('tests')
run_all_tests
```

Or:

```matlab
results = runtests('tests', 'IncludeSubfolders', false);
assertSuccess(results);
```

## Conventions

- Class-based only (`matlab.unittest.TestCase`)
- File names use `t` prefix
- Arrange / Act / Assert
- No logic branches inside test methods
- Paths added in `TestClassSetup`, removed via teardown
- Experiment-local opinion tests stay under each `experiments/copy*/tests/`
