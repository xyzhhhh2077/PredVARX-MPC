# PredVAR + MPC MATLAB experiments — organized workspace

> **2026-07-18 归档位置**：`代码/archive/organized_legacy_2026-07-10/`  
> **角色**：旧谱系实验归档（01–10）。**现行主仓库**请用 `代码/PredVARX-MPC/`。  
> 旧字母体系（Q/S/Z 等）与现行 Git 同名实验含义不同，见 `代码/PredVARX-MPC/LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md`。

## How to run an experiment

From MATLAB, add this archive workspace to the path once:

```matlab
root = 'E:\academic_files\phd-learning\代码\archive\organized_legacy_2026-07-10';
addpath(root);
addpath(fullfile(root, '02_identification'));
```

Then use the wrapper directly. It adds the canonical identification routine and runs a selected script in that script's own folder:

```matlab
root = 'E:\academic_files\phd-learning\代码\archive\organized_legacy_2026-07-10';
run_experiment(fullfile(root, '06_delta_mpc', 'copyD', 'copyD_deltampc.m'));
```

Outputs created by an experiment stay beside that experiment under its `results/` directory whenever the source script has been updated to use it. **Existing copied MAT/PNG files are historical artifacts; do not treat them as a fresh rerun.**

## Folder map

| Folder | Contents | Role |
|---|---|---|
| `01_baseline/` | Original `main/` source | Clean, read-only baseline snapshot |
| `02_identification/` | Canonical IVR/PredVARX identifier and quick interface test | Shared identification implementation |
| `03_mpc_variants/` | Q, QZ, and PredVAR-aligned MPC variants | Deterministic MPC experiment lineage |
| `04_soft_constraints/` | Z/Zfix and sigma-soft-constraint variants | Sliding-window and soft-constraint experiments |
| `05_smpc/` | SMPC, oracle, PCA, and retrospective tests | SMPC diagnosis and comparison experiments |
| `06_delta_mpc/` | Incremental-control Δ-MPC experiment | Δu formulation experiment |
| `07_reference_cases/` | Final C and Lorenz reference cases | Reference/sanity simulations |
| `08_presentation/` | Slide-generation asset | Presentation support |
| `09_metadata/` | Exported metadata | Auxiliary non-MATLAB data |

## Important conventions

- `01_baseline/main/` must not be changed for new designs.
- Each named copy is historical experimental evidence, not automatically a validated final method.
- A few old scripts write generic names such as `copyZ_data.mat` or `copyQ_data.mat`. They were separated into individual experiment folders to stop collisions. Run them through `run_experiment.m` so their relative `load`/`save` paths resolve locally.
- `predvarx_identify.m` in `02_identification/` is the canonical shared copy. It is byte-identical to the original root and baseline copy at the time of organization.
- The parent `matlab/README.md` is stale: it references `copyK_final.m`, which is not present in this directory snapshot.

## Known data lineage conflicts preserved as historical evidence

- `copyQ_dual_dinkla.m` and `copyQ_save.m` both use `copyQ_data.mat`.
- `copyQZ_soft.m` and `copyQZ_v2.m` both use `copyQZ_data.mat` but use different identifier implementations.
- `copyZfix.m`, `copyZfix_nosoft.m`, `copyZfitest.m`, and baseline `copyS_nosoft.m` use the `copyZfix_data.mat` name.

They have therefore been placed in separate leaf folders. This preserves the original artifacts without allowing one script to silently overwrite another experiment's data.
