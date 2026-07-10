# MATLAB 实验目录总览

> **仓库**：`matlab 2026,7,9`
> **原则**：每个版本目录只放本版本的主代码、辅助函数、局部测试、结果和说明；基础版本不被实验版本覆盖。

## 目录结构

```text
main/                              # 基础版本
├── copyS_nosoft.m
├── predvarx_identify.m
├── README.md
└── results/

experiments/
├── copyO_oblique/                 # 斜投影双基诊断副本
│   ├── copyO_oblique.m
│   ├── predvarx_identify_oblique.m
│   ├── tests/
│   ├── results/
│   └── README.md
├── copyP_centered_smpc/           # 全阶：控制逻辑正确性验证
│   ├── copyP_centered_smpc.m
│   ├── control_ready_subspace_varx.m
│   ├── centered_smpc_step.m
│   ├── tests/
│   ├── results/
│   └── README.md
└── copyQ_control_aware/           # 低阶：控制感知子空间版本
    ├── copyQ_control_aware.m
    ├── control_aware_subspace_varx.m
    ├── tests/
    ├── results/
    └── README.md
```

## 基础版本与副本差异

| 版本 | 相对基础版本的主要变化 | 目的 | 本次实测结论 |
|:---|:---|:---|:---|
| `main` | 原始 QR 正交化、差分 MPC、`alpha_cc=0.45` | 保存原始对照 | 只作历史基线；MAT 现在额外保存辨识快照，便于诊断 |
| `copyO_oblique` | SVD 双基，$R^TP=I$；静态残差改用 $\bar R^T$ | 验证斜投影数学身份 | 双基恒等式正确，但旧在线中心化/SMPC 逻辑仍导致反向跟踪 |
| `copyP_centered_smpc` | 全阶 $\ell=n$；中心化坐标；绝对预测 MPC；Boole 联合机会约束 | 先验证控制逻辑本身 | $y_1,y_2$ MAE 为 0.093/0.087，零越界，QP 1200/1200 可行 |
| `copyQ_control_aware` | 低阶 $\ell=4$；强制保留 $y_1,y_2$ 输出轴；其余方向按数据取子空间 | 在降维下保留控制方向 | MAE 0.114/0.093，零越界，QP 1200/1200 可行 |

## 如何运行

在仓库根目录：

```matlab
run('main/copyS_nosoft.m')
run('experiments/copyO_oblique/copyO_oblique.m')
run('experiments/copyP_centered_smpc/copyP_centered_smpc.m')
run('experiments/copyQ_control_aware/copyQ_control_aware.m')
```

> 现在主脚本会直接将新的 `.mat` 和 `.png` 写入各自 `results/`；这里保存的是可覆写的运行快照，避免生成物混入代码目录根部。

## 测试入口

```matlab
addpath('experiments/copyO_oblique'); addpath('experiments/copyO_oblique/tests');
test_predvarx_identify_oblique

addpath('experiments/copyP_centered_smpc'); addpath('experiments/copyP_centered_smpc/tests');
test_centered_smpc_step

addpath('experiments/copyQ_control_aware'); addpath('experiments/copyQ_control_aware/tests');
test_control_aware_subspace
```
