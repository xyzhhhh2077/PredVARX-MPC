# PredVARX-MPC

面向 **PredVAR / PredVARX、斜投影潜变量辨识与随机模型预测控制（SMPC）** 的 MATLAB 实验仓库。

仓库采用“一个 `copy*` 目录对应一个独立实验版本”的组织方式。每个副本保留自己的算法、测试、指标和图片，便于区分：

- 原始构造；
- 理论或实现修复；
- 诊断性消融；
- 新方法试验；
- 尚未完成的研究草稿。

> 不要仅按字母判断算法先后，也不要把不同副本的指标直接排名。只有 plant、训练数据、参考、控制器、噪声与闭环时域一致时，数值才可公平比较。

## 当前研究主线

### 1. PredVARX 几何与 SMPC

研究低维潜变量模型

$$
z_{k+1}=A z_k+B(u_k-\bar u)+\varepsilon_{k+1},
\qquad \hat y_k=\bar y+Pz_k,
$$

以及斜投影双基关系

$$
R^\top P=I,
$$

如何与控制目标、预测误差协方差和机会约束保持一致。

主要目录：

- `copyO`–`copyZ`：正交/斜投影辨识、中心化 SMPC、控制感知子空间和公平消融；
- `copyAA`–`copyAM`：分块双基、意见修复、闭环诊断、协方差与保证层探针；
- `copyOPQSTR_unified`、`copyALL_unified`：统一工况对照。

### 2. CRTE 文稿实验算法

[`copyAR_crte_paper_spectral_validation_unknown_noise`](experiments/copyAR_crte_paper_spectral_validation_unknown_noise/) 是当前文稿 Sec. 3.3–3.4 对应的实验路径：

- 固定谱 surrogate；
- paper trace normalization；
- percentile 门控；
- validation-only 选型；
- 自由输出一步残差协方差 proxy；
- `uses_true_Sigma_n = 0`。

它不是 `copyAO` 的 profiled-teacher/min-teacher 路径，也不是已知真实 $\Sigma_n$ 的 Oracle。

单种子1200步结果：

| MAE $y_1$ | MAE $y_2$ | QP成功率 | fallback/不可行 |
|---:|---:|---:|---:|
| 0.06048 | 0.05706 | 100% | 0 |

### 3. AR 之后：学习任务方向与偏好

这些实验全部复用 copyAR 的1500个离线样本，不增加训练集。它们是 CRTE 之外的扩展或诊断，不应写成文稿原算法。

| 副本 | 学习/选择对象 | 1200步结果 | 证据边界 |
|---|---|---|---|
| [`copyAS`](experiments/copyAS_learned_task_anchor_smpc/) | 从完整输出任务参考学习混合任务锚 $E_T$ | MAE `[0.1962, 0.4163]`；QP 100%；fallback 0 | 学习锚与固定任务平面主角度 `[55.73°, 87.86°]`，未改善原任务跟踪 |
| [`copyAT`](experiments/copyAT_learned_output_directions/) | fixed / supervised / input-authority 三分支公平比较 | fixed=sup=`[0.0605, 0.0573]`；authority=`[0.2969, 0.3294]` | 任务输出与最易被输入推动的方向不是同一概念 |
| [`copyAU`](experiments/copyAU_soft_preference_output/) | 软偏好矩阵与输入权威共同学习两个输出组合 | MAE `[1.0600, 1.2582]`；QP 100%；fallback 0 | 数值求解成功，但没有实现给定任务跟踪；负面实验 |
| [`copyAV`](experiments/copyAV_hard_preference_output/) | 硬锁权重最高的原始输出轴，当前为 $y_1,y_2$ | MAE `[0.05973, 0.06185]`；QP 100%；fallback 0 | 单种子恢复良好跟踪，不证明普遍最优 |

copyAT统一使用：

- 相同旧 `(u,y)`；
- 相同随机种子和闭环扰动；
- 1200步、5段参考、每段240步；
- 相同控制器与约束；
- 原始物理输出 `y1,y2` 评价；
- `new_training_samples = 0`。

![copyAT三种输出方向比较](experiments/copyAT_learned_output_directions/results/copyAT_learned_output_directions_fig.png)

![copyAU软偏好结果](experiments/copyAU_soft_preference_output/results/copyAU_soft_preference_output_fig.png)

![copyAV硬偏好结果](experiments/copyAV_hard_preference_output/results/copyAV_hard_preference_output_fig.png)

## CRTE 路径不要混用

| 路径 | 目录 | 定位 |
|---|---|---|
| 固定谱 + validation | `copyAN_crte_fixed_surrogate` → `copyAR_crte_paper_spectral_validation_unknown_noise` | 文稿实验算法 |
| profiled teacher | `copyAO_crte_teacher_profiled_unknown_noise` | 结构研究线，不是文稿 validation 选型 |
| 多步 task | `copyAP_crte_multistep_task_20x4` | $t+1:t+H$ task 与 blocked proxy 扩展 |
| 高阶 VARX | `copyAQ_crte_varx_order_20x4` | $s>1$ 草稿；不是完成的20×4全网格 |
| 学习输出方向 | `copyAS`–`copyAV` | 非原文扩展与偏好实验 |

## 仓库结构

```text
PredVARX-MPC/
├── main/                         # 历史基线
├── experiments/                  # 现行独立实验副本
│   ├── copyO_oblique/ ... copyAM_tracked_cov_only/
│   ├── copyAN_crte_fixed_surrogate/
│   ├── copyAO_crte_teacher_profiled_unknown_noise/
│   ├── copyAP_crte_multistep_task_20x4/
│   ├── copyAQ_crte_varx_order_20x4/
│   ├── copyAR_crte_paper_spectral_validation_unknown_noise/
│   ├── copyAS_learned_task_anchor_smpc/
│   ├── copyAT_learned_output_directions/
│   ├── copyAU_soft_preference_output/
│   ├── copyAV_hard_preference_output/
│   ├── copyOPQSTR_unified/
│   └── copyALL_unified/
├── archive/                      # 非现行历史版本
├── tests/                        # 根目录 matlab.unittest
├── docs/version-summary/         # 版本图文说明
└── README.md
```

`copyAN_closed_loop_audit_grade` 只是空占位目录，不代表完成实验。

## 环境

- MATLAB R2024a 或更高版本；
- Optimization Toolbox（SMPC 的 `quadprog`）；
- Git LFS（结果 `.mat`、`.png`、`.fig` 等大文件）。

克隆：

```bash
git lfs install
git clone https://github.com/xyzhhhh2077/PredVARX-MPC.git
cd PredVARX-MPC
```

未安装 Git LFS 时，大文件只会下载指针。

## 运行 AR–AV

在 MATLAB 中进入仓库根目录：

```matlab
repo = pwd;

% copyAR：文稿实验算法
cd(fullfile(repo,'experiments','copyAR_crte_paper_spectral_validation_unknown_noise'))
run('tests/test_crte_paper_spectral_varx.m')
copyAR_crte_paper_spectral_validation_unknown_noise

% copyAS：学习式任务锚
cd(fullfile(repo,'experiments','copyAS_learned_task_anchor_smpc'))
run('tests/test_learned_task_anchor_varx.m')
copyAS_learned_task_anchor_smpc

% copyAT：三种输出方向公平比较
cd(fullfile(repo,'experiments','copyAT_learned_output_directions'))
run('tests/test_learn_output_directions.m')
copyAT_learned_output_directions

% copyAU / copyAV：软偏好与硬偏好
cd(fullfile(repo,'experiments','copyAU_soft_preference_output'))
runtests('tests/learnPreferredOutputDirectionsTest.m')
copyAU_soft_preference_output

cd(fullfile(repo,'experiments','copyAV_hard_preference_output'))
runtests('tests/selectHardPreferenceOutputsTest.m')
copyAV_hard_preference_output
```

各副本把结果写入自己的 `results/`：

- `*_metrics.txt`：可审计指标；
- `*_data.mat`：运行数据；
- `*_fig.png`：结果图。

## 根目录测试

```matlab
cd('tests')
run_all_tests
```

当前根测试覆盖：

- `tPredvarxIdentifyOblique`：斜投影双基恒等式；
- `tCenteredSmpcStep`：中心化预测与机会约束收紧；
- `tControlAwareSubspace`：降阶空间中的被控轴保持。

各 `copy*` 的专用合同测试位于对应 `tests/`。

## 结果解释规则

1. **数学证明、代码合同和单次仿真分开表述。**
2. `QP成功率=100%`、`fallback=0` 只说明该轨迹由主QP路径完成，不证明递归可行或稳定。
3. residual covariance 是未知噪声下的工程 proxy，不是真实传感器噪声。
4. Oracle 参数只能用于仿真对照，不能包装成未知噪声算法结论。
5. 输入权威高表示方向容易被输入推动，不表示它就是任务需要跟踪的输出。
6. 单种子结果不证明跨工况泛化、全局最优或普遍性能提升。

## 其他文档

| 文档 | 内容 |
|---|---|
| [`VERSION_EVOLUTION_AND_METRICS.md`](VERSION_EVOLUTION_AND_METRICS.md) | O–Z 演化与指标 |
| [`LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md`](LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md) | 旧字母复用与现行目录辨析 |
| [`docs/version-summary/PREDVARX_MPC_VERSION_SUMMARY.md`](docs/version-summary/PREDVARX_MPC_VERSION_SUMMARY.md) | 带图版本总览 |
| [`archive/README.md`](archive/README.md) | 历史归档说明 |

## 数据与版本管理

- `.mat`、`.png`、`.fig`、`.pdf` 由 Git LFS 跟踪；
- 新试验应建立新副本，不覆盖已有实验身份；
- 不同协议的指标不得直接合并为排名；
- 汇报时需注明 copy、随机种子、训练数据、噪声口径、闭环长度及是否完整运行。
