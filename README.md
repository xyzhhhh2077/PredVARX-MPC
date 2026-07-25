# PredVARX-MPC

面向 **PredVAR / PredVARX + SMPC** 的公开 MATLAB 实验副本库。

- **GitHub**：https://github.com/xyzhhhh2077/PredVARX-MPC
- **可见性**：public（公开）
- **默认分支**：`main`
- **大文件**：实验 `.mat` / `.png` / `.fig` 通过 **Git LFS** 跟踪

> **原则**：每个 `copy*` 目录自成一版（代码、辅助函数、局部测试、结果、说明）。  
> **不要**把不同 copy 混成一个算法。归档区里的 Q/S/Z 等字母可能与现行目录**同名不同义**。

---

## 克隆

```bash
git lfs install
git clone https://github.com/xyzhhhh2077/PredVARX-MPC.git
cd PredVARX-MPC
```

未安装 LFS 时，`.mat`/图片只会得到指针文件，不是真实数据。

---

## 仓库结构

```text
PredVARX-MPC/
├── main/                         # 历史基线（copyS_nosoft）
├── experiments/                  # 现行 copy 套件（优先入口）
│   ├── copyO_oblique/ ...
│   ├── copyAO_crte_teacher_profiled_unknown_noise/
│   ├── copyAP_crte_multistep_task_20x4/
│   └── ...
├── archive/                      # 非现行历史 MATLAB 版本
│   ├── early_snapshot_2026-07-09/
│   ├── organized_legacy_2026-07-10/
│   └── routeC_delta_inc/
├── tests/                        # 根目录 matlab.unittest 套件
├── docs/version-summary/         # 带图的版本说明
├── VERSION_EVOLUTION_AND_METRICS.md
├── LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md
└── README.md                     # 本文件
```

---

## 相关文档

| 文档 | 内容 |
|---|---|
| [`VERSION_EVOLUTION_AND_METRICS.md`](VERSION_EVOLUTION_AND_METRICS.md) | 现行 O–Z 演化、公平比较口径、指标总表 |
| [`LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md`](LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md) | 旧工作区字母复用 vs 现行 Git 目录名 |
| [`docs/version-summary/PREDVARX_MPC_VERSION_SUMMARY.md`](docs/version-summary/PREDVARX_MPC_VERSION_SUMMARY.md) | 逐版说明 + PNG |
| [`archive/README.md`](archive/README.md) | 归档内容与使用边界 |

---

## 现行实验（`experiments/`）

新工作优先放这里。下列目录均为已发布树中的现行副本。

### A. 辨识 / 几何 / 控制主链

| 目录 | 作用 |
|---|---|
| `main/` | 历史基线：`copyS_nosoft` + `predvarx_identify` |
| `copyO_oblique` | 斜投影双基恒等式诊断 |
| `copyP_centered_smpc` | 中心化绝对预测 SMPC（全阶控制逻辑验证） |
| `copyQ_control_aware` | 低阶、强制保留被控轴的控制感知正交子空间 |
| `copyR_moqin_oblique` | Mo–Qin Algorithm-1 风格白化 realization 基线 |
| `copyT_process_lv_smpc` | 高维 process-like plant + 在线协方差尺度 |
| `copyU_smooth_noise_smpc` | 平滑时变噪声压力测试（copyT 风格） |
| `copyV_iterative_ivr` | 被控补空间迭代 IVR（正交） |
| `copyW_fair_identifier_compare` | 固定 plant/controller 的多辨识器公平消融 |
| `copyX_control_aware_oblique` | 控制感知、轻度斜读取 $R_\alpha$（无白化） |
| `copyY_no_whitening_direct_update` | 显式标注 copyX 真实路径；**不是新算法** |
| `copyZ_strict_whitened_copyX_plant` | 在 copyX plant 上的白化消融（不是完整 Alg.1 主张） |

### B. 意见修复 / 闭环与保证层探针

| 目录 | 作用 |
|---|---|
| `copyAA_split_control_free_oblique` | 控制/自由分块斜投影构造 |
| `copyAB_opinion_fixes` | 意见硬修合集 + 单元测试 |
| `copyAC_open_problems_trial` | 开放问题试验，并与 AB 公平对照 |
| `copyAD_closure_ladder` | 闭环/证书阶梯诊断 |
| `copyAE_stress_calibration` | 压力与校准探针 |
| `copyAF_op9_and_p1pp` | 意见 9 / Prop-1 相关探针 |
| `copyAG_original_alpha_cert` | 原始 alpha 证书实验 |
| `copyAH_multistep_alpha_backup` | 主路径失败后的多步 alpha 备份 |
| `copyAI_crosscov_chance` | 交叉协方差机会约束探针 |
| `copyAJ_safety_filter_cert` | 安全滤波器证书探针 |
| `copyAK_terminal_set_rf` | 终端集 / 递归可行探针 |
| `copyAL_split_empirical_cov_oblique` | 经验协方差分块斜投影变体 |
| `copyAM_tracked_cov_only` | 仅被控轴协方差的机会约束路径 |
| `copyAN_crte_fixed_surrogate` | CRTE 固定谱代理（自由补空间） |
| `copyAN_closed_loop_audit_grade` | 占位目录（无完整运行；仅保留分支名连续性） |

### C. 统一对照包

| 目录 | 作用 |
|---|---|
| `copyOPQSTR_unified` | 同一 plant/噪声/SMPC 下还原 O/P/Q/R/S/T 算法差异 |
| `copyALL_unified` | 更广的多 copy 统一对照 + 分图 |

### D. CRTE profiled-teacher 线（当前研究前沿）

| 目录 | 作用 | 状态 |
|---|---|---|
| `copyAO_crte_teacher_profiled_unknown_noise` | 完整 profiled teacher、unknown-noise proxy、严格 FWL SVD 支撑、单次全长运行 | **CRTE 主结构结果** |
| `copyAP_crte_multistep_task_20x4` | 多步 task 堆叠 $t+1:t+H$、$\mu=1$、blocked forward noise proxy | 结构测试 + smoke + 一次 $H=3$ 全长确认 |
| `copyAQ_crte_varx_order_20x4` | VARX 阶数 / companion 草稿，用于 $s>1$ 消融 | 核心与测试草稿；**不是**完成的 20×4 全网格 |

**CRTE 阅读顺序**

1. 几何 / FWL / teacher 结构 → `copyAO`
2. 多步 task + blocked proxy 合同 → `copyAP`
3. 更高阶潜变量 AR 草稿 → `copyAQ`（未完成）

没有配对协议时，不要把 AO/AP/AQ 数字揉成“最优方法”结论。

---

## 归档（`archive/`）

非现行历史 MATLAB。新工作请只改 `experiments/` 与 `main/`。

| 路径 | 内容 |
|---|---|
| `archive/early_snapshot_2026-07-09/` | 早期独立快照：main + O/P/Q |
| `archive/organized_legacy_2026-07-10/` | 按 01–10 主题归档的旧实验（旧 Q/S/Z/SM/D/oracle 字母） |
| `archive/routeC_delta_inc/` | v10 Δ-tracking 加速线（C/E/F/G/H/I/P1–P4） |

旧字母复用见 `LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md`。  
例如：归档里的 `copyQ_*` **不等于** 现行 `copyQ_control_aware`。

---

## 如何运行（MATLAB R2024a+）

在仓库根目录：

```matlab
% 历史基线
run('main/copyS_nosoft.m')

% 现行示例
run('experiments/copyP_centered_smpc/copyP_centered_smpc.m')
run('experiments/copyX_control_aware_oblique/copyX_control_aware_oblique.m')

% CRTE 主结构
cd('experiments/copyAO_crte_teacher_profiled_unknown_noise')
test_crte_profiled_teacher_unknown_noise
copyAO_crte_teacher_profiled_unknown_noise

% CRTE 多步 task 测试
cd('../copyAP_crte_multistep_task_20x4')
runtests('tests/test_copyAP_multistep_task.m')
```

各现行 copy 会把快照写到各自 `results/`（`.mat` / `.png` / 指标文本）。

---

## 根目录测试

```matlab
cd('tests')
run_all_tests
```

| 测试类 | 覆盖 |
|---|---|
| `tPredvarxIdentifyOblique` | 斜投影双基恒等式 |
| `tCenteredSmpcStep` | 中心化预测 + Boole 收紧 |
| `tControlAwareSubspace` | 降维下被控轴覆盖 |

各 copy 自有的意见/诊断测试在 `experiments/copy*/tests/`。

---

## 数据与图片策略

- **Git LFS 跟踪**：实验/归档结果下的 `*.mat`、`*.png`、`*.fig`
- **不跟踪**：文献 PDF 缓存、散落的 `experiments/*.pdf`、本地编辑器小工具
- 结果图是**运行证据**，不等于自动可投稿的“最终主图”
- 除非 plant、噪声、时域、控制器协议一致，否则**不能**跨 copy 直接比 MAE 排名

---

## 先看哪里

| 目标 | 打开 |
|---|---|
| 了解 O–Z 演化 | `VERSION_EVOLUTION_AND_METRICS.md` |
| 公平 OPQSTR 对照 | `experiments/copyOPQSTR_unified/` |
| 当前 CRTE teacher 结构 | `experiments/copyAO_crte_teacher_profiled_unknown_noise/` |
| 多步 task 扩展 | `experiments/copyAP_crte_multistep_task_20x4/` |
| 旧字母考古 | `LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md` + `archive/` |
| 带图总览 | `docs/version-summary/PREDVARX_MPC_VERSION_SUMMARY.md` |

---

## 备注

- 多数 copy 需要把对应目录 `addpath` 进 MATLAB（有时还要加 `tests/`）。
- SMPC 主路径通常需要 Optimization Toolbox 的 `quadprog`。
- 并行池可选；部分后续脚本更倾向多进程 launcher，而不是超大 parpool。
- 汇报结果时务必写明**确切 copy 目录**，以及图/mat 属于 smoke、单种子还是完整战役。
