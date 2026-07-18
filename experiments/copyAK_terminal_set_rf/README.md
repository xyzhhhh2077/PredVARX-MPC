# copyAK_terminal_set_rf

## Purpose

Opinion 9 **numerical probes** for terminal cost / coarse terminal set / RF condition chain.

- Documents minimal RF chain C1–C6 and marks what is **unproved**
- Optional `use_terminal_cost` (default OFF): `Pterm ≈ dlyap(A',Qf)`
- Optional `use_terminal_set` (default OFF): spectral-box linear rows for
  `z_N' Pterm z_N ≤ α_term` (α calibrated or scanned)
- Closed-loop OFF vs COST vs SET at medium (`y_max=0.70`) and tight (`0.50`), multi-seed

**Non-claims:** no chance-constraint recursive feasibility theorem; no closed-loop stability theorem; box ≠ exact ellipsoid Xf.

Does not modify `main/` or other experiment copies.

## Run

```matlab
cd('.../experiments/copyAK_terminal_set_rf');
addpath(pwd); addpath('lib');
run('copyAK_terminal_set_rf.m');
```

Focused test:

```matlab
cd('.../experiments/copyAK_terminal_set_rf/tests');
runtests('test_opinion9_terminal');
```

## Outputs

- `results/copyAK_metrics.csv`
- `results/copyAK_terminal_set_rf_report.md`
- `results/copyAK_terminal_set_rf_data.mat`
- `patches/opinion9_terminal.md`
