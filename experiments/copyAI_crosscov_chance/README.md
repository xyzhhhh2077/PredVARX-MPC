# copyAI_crosscov_chance

Opinion 6 experiment: inject latent/residual **cross-covariance** into chance-constraint `Sigma_y`, with fair OFF vs ON comparison.

## Math (drop-cross gap)

With `y_c = P z + o`, `z = R' y_c`, dual `R' P = I`:

```text
Sigma_y_full = P Sigma_z P' + P Sigma_zo + Sigma_zo' P' + Sigma_o
Sigma_y_drop = P Sigma_z P' + Sigma_o          % default SMPC
```

`R' o = 0` samplewise does **not** imply `Sigma_zo = 0`.

## Layout

| Path | Role |
|---|---|
| `lib/cross_cov_diagnostics.m` | z,o,Σz,Σo,Σzo, drop/full recon, optional H variances |
| `lib/sigma_y_chance_blocks.m` | build full/drop Σy + hq'Σy hq ratios |
| `centered_smpc_step.m` | `opt.use_cross_cov` (default **false**) |
| `run_crosscov_config.m` | one closed-loop config |
| `copyAI_crosscov_chance.m` | fair multi-seed OFF vs ON runner |
| `tests/test_opinion06_cross_cov_chance.m` | focused unit test |
| `results/` | CSV, REPORT.md, MAT |
| `patches/opinion6.md` | PASS/FAIL + stdout |

## Non-claims

- Does **not** re-prove Boole risk allocation under full cross blocks.
- Offline `Sigma_zo` only; `Sigma_obs` still a proxy.

## Run

Use ASCII launcher + `char([20195 30721])` for `代码` path on Chinese Windows.
