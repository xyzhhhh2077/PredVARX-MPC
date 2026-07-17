# Opinion 1 patch report

**状态：PASS**

## 意见内容

几何层（reconstruction-layer）应对 tracked 左右/列/对偶残差做强制断言，
确保 `tracked_left` / `tracked_right` / `tracked_column` / `dual` 数值为零，
而非仅事后打印诊断。

## 做了什么

仅修改 `experiments/copyAB_opinion_fixes/`，未触碰 `main/`、`copyAA/`、`copyX/`。

### 1. `split_control_free_ivr_varx.m`

在函数末尾（stats 几何残差字段已计算之后）增加 Opinion 1 强制断言：

| 残差字段 | 含义 | 容差 |
|---|---|---|
| `stats.tracked_left_error` | `‖E'PR' − E'‖_F` | `1e-8` |
| `stats.tracked_right_error` | `‖PR'E − E‖_F` | `1e-8` |
| `stats.tracked_column_error` | `‖R(:,1:q) − E‖_F` | `1e-8` |
| `stats.dual_error` | `‖R'P − I‖_F` | `1e-8` |

可选开关：

- `stats.enforce_geometry`（默认 `true`）
- 同时支持 name-value 入参：`'enforce_geometry', true|false`
- 默认开启；关闭时仍报告残差，但不 `assert`

### 2. `tests/test_opinion01_geometry.m`

- 合成异构相关传感器噪声数据，调用 `split_control_free_ivr_varx`
- 验证 `enforce_geometry` 默认为 `true`
- 显式断言 **left=0**、**column=0**（以及 right/dual）
- 校验 tracked 列 `R(:,1:q)=E`、`P(:,1:q)=E`
- 验证可关掉开关：`'enforce_geometry', false` 不抛错且残差仍数值为零

## MATLAB 运行

```text
命令:
/e/MATLABinhere/bin/matlab -batch \
  "cd('E:/academic_files/phd-learning/代码/PredVARX-MPC/experiments/copyAB_opinion_fixes'); \
   addpath(pwd); addpath(fullfile(pwd,'tests')); \
   test_opinion01_geometry"

stdout:
---RUN---
PASS opinion01 geometry: left=0.000e+00 right=0.000e+00 col=0.000e+00 dual=2.720e-15 (enforce=1, free_oblique=2.165e+00)

exit code: 0
```

## 结论

| 检查项 | 结果 |
|---|---|
| left residual ≈ 0 | **PASS** (`0.000e+00`) |
| column residual ≈ 0 | **PASS** (`0.000e+00`) |
| right residual ≈ 0 | **PASS** (`0.000e+00`) |
| dual residual ≤ 1e-8 | **PASS** (`2.720e-15`) |
| 内部 assert 默认开启 | **PASS** (`enforce=1`) |
| 可选开关可关闭 | **PASS** |

**总体：PASS**

## 边界说明

几何保证仅是 **重建层**（`P,R` 与 tracked 基 `E` 的代数关系），
**不** 意味着闭环跟踪 RMSE=0。
