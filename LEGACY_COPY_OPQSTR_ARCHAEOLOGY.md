# 旧版 copyO / P / Q / S / T / R 版本考古与现行 Git 对照

> 目的：区分两套曾经重用字母命名的实验。第一套是早期 MATLAB 工作区中的 Q/S/SM/QZ 等控制副本；第二套是 2026-07-09 Git 仓库中重新定义的 copyO–copyZ 实验。相同字母不一定表示同一算法。

## 1. 为什么此前总表看起来少了 OPQSTR

现行 Git 已经包含：

- `copyO_oblique`
- `copyP_centered_smpc`
- `copyQ_control_aware`
- `main/copyS_nosoft`
- `copyT_process_lv_smpc`
- `copyR_moqin_oblique`

但旧工作区 `E:/academic_files/phd-learning/PredVAR+MPC/matlab/` 还保存了另一批更早的 Q/S/SM/QZ 副本。这些副本和现行 Git 中同字母版本不是同一含义。因此不能简单把“copyQ”写成唯一版本。

## 2. OPQSTR 主版本对照

| 字母 | 现行 Git 版本 | 核心目的 | 与旧工作区的关系 |
|---|---|---|---|
| O | `experiments/copyO_oblique` | 修复斜投影双基，验证 $R^TP=I$ | 有完整源码、README、测试和结果；没有发现另一份不同算法的独立旧 O 主源码。 |
| P | `experiments/copyP_centered_smpc` | 中心化全阶子空间、绝对预测、Boole 机会约束 | 有完整源码与理论笔记。历史 session 提到过 `copyP_fast_reid.m`，但当前磁盘没有找到该源码，不能按真实版本补造。 |
| Q | `experiments/copyQ_control_aware` | 低阶控制感知子空间，强制保留 $e_1,e_2$ | 与旧工作区的 `copyQ_save`、`copyQ_dual_dinkla` 不是同一算法；旧 Q 是在线重辨识/噪声估计路线。 |
| S | `main/copyS_nosoft.m` | 原始基线：滑动窗口重辨识、Pbar 同步和旧差分 MPC/SMPC | 旧工作区 `01_baseline/main/copyS_nosoft.m` 是同一路线；另外还有 `copyS_delta` 软约束分支。 |
| T | `experiments/copyT_process_lv_smpc` | p=30 工业过程、控制感知正交子空间、在线协方差更新 | 是后期新定义的 T；未发现另一份独立旧 T 主源码。 |
| R | `experiments/copyR_moqin_oblique` | Mo–Qin 白化—IVR—反白化诊断基线 | 工作区 `10_strict_moqin_predvar/copyR_moqin_oblique` 是同步副本。 |

## 3. 旧 Q 系列

### 3.1 `copyQ_save`

路径：

```text
matlab/03_mpc_variants/copyQ_save/copyQ_save.m
```

主要设定：

- $p=15,\ell=2,T_{cl}=10000$；
- 每 50 步重新辨识；
- 噪声每 1000 步分段变化；
- 预测时域按阶段取 5、12、30；
- 保存 `copyQ_data.mat`；
- 只做输入边界 QP，没有现行 copyQ 的 control-aware tracked-axis 构造。

该版本名为 Q，但它不是 `copyQ_control_aware`。

### 3.2 `copyQ_dual_dinkla`

路径：

```text
matlab/03_mpc_variants/copyQ_dual_dinkla/copyQ_dual_dinkla.m
```

在旧 Q 上增加：

- 50 步在线重辨识；
- 潜空间创新协方差 $\Sigma_\varepsilon$ 的窗口估计；
- 静态残差协方差 $\Sigma_{\bar e}$ 的窗口估计；
- 输出协方差
  $$
  \Sigma_y=P\Sigma_zP^T+\bar P\Sigma_{\bar e}\bar P^T;
  $$
- Li/Dinkla 风格 chance constraints。

它是“在线双协方差估计 Q”，不是现行“控制方向覆盖 Q”。

### 3.3 `copyQZ_v2` 与 `copyQZ_soft`

路径：

```text
matlab/03_mpc_variants/copyQZ_v2/copyQZ_v2.m
matlab/03_mpc_variants/copyQZ_soft/copyQZ_soft.m
```

这两项比较旧 Q 与 Z 型软约束：

- Q：纯 MPC；
- QZ：根据当前输出、噪声尺度或非对称规则修改第一步目标；
- 目标是降低 $y_2>y_{max}$ 的越界。

它们属于启发式软约束试验，不是现行 copyZ 白化消融。

## 4. 旧 S 系列

### 4.1 `copyS_nosoft`

现行 Git 路径：

```text
main/copyS_nosoft.m
```

工作区同步路径：

```text
matlab/01_baseline/main/copyS_nosoft.m
```

核心特征：

- $p=15,\ell=2,T_{cl}=2000$；
- 比较累计重辨识 C 与最近窗口 SMPC-only/Zfix；
- 每 30 步重辨识；
- Z 分支使用离线数据加最近在线数据；
- 重辨识后同步更新 $P, R, A, B$ 和预测矩阵 $M_j,N_j$；
- 同步更新 $P_{bar}$；
- 仍保留旧差分 MPC 与 `alpha_cc=0.45` 负分位数问题。

所以 S 是历史在线重辨识基线，不是现代 copyT/X 的固定结构加在线协方差更新。

### 4.2 `copyS_delta`

路径：

```text
matlab/04_soft_constraints/copyS_delta_sigma/copyS_delta.m
```

在 S 路线上加入：

$$
\Delta\sigma_k=\sigma_k-\sigma_{k-1},
\qquad
soft_k=\max(0,\Delta\sigma_k).
$$

噪声估计上升时，通过启发式目标偏移进行收紧。它不是严格概率机会约束。

## 5. 旧 SM / SMPC 分支

### `copySM_smpc`

路径：

```text
matlab/05_smpc/copySM/copySM_smpc.m
```

核心变化：

- 使用滑动窗口重新辨识；
- 同步更新 $P_{bar}$；
- 传播潜空间协方差；
- 加入双侧机会约束；
- `alpha_cc=0.84`，分位数为正；
- 取消启发式 soft，令 `soft=0`。

它是从旧 S/QZ 路线向正式 SMPC 过渡的版本，但仍采用旧差分目标和旧坐标体系，不能与后来的 copyP/T/X 等同。

## 6. O/P/Q/R/T 的理论笔记证据

现有理论笔记：

```text
论文笔记/.../01_v10理论主文档/
├── copyO_oblique：斜投影双基实验.md
├── copyP：居中全阶子空间SMPC验证.md
├── copyQ：低阶控制感知子空间SMPC.md
├── copyR：严格Mo-Qin斜投影PredVARX-SMPC基线.md
└── copyT 工业过程 PredVARX-SMPC 理论推导.md
```

这些笔记对应现行 Git 的 O/P/Q/R/T，而不是所有旧工作区同名字母实验。

## 7. 字母重名的最终解释

| 名称 | 旧工作区含义 | 现行 Git 含义 |
|---|---|---|
| Q | 50 步重辨识、时变噪声、双 Dinkla 或 QZ 比较 | 低阶控制感知正交子空间 |
| S | 滑动窗口重辨识/Zfix 基线；另有 $\Delta\sigma$ soft | `main` 历史基线 |
| Z | 多个软约束/Zfix/QZ 启发式版本 | copyZ：copyX plant 上的白化消融 |
| O | 斜投影双基诊断 | 同一条主线，现已进入 Git |
| P | 历史 session 曾提及 fast re-id，但源码未找到 | 中心化全阶 SMPC |
| R | Mo–Qin 白化对照 | 同一条主线，现已进入 Git |
| T | 未发现独立旧 T | p=30 工业过程验证 |

因此，今后引用版本时应同时写完整目录名，例如：

- `legacy copyQ_dual_dinkla`；
- `current copyQ_control_aware`；
- `legacy copyZfix`；
- `current copyZ_strict_whitened_copyX_plant`。

只写“copyQ”或“copyZ”会产生歧义。

## 8. 证据边界

- O/P/Q/R/T 现行版本：源码、README、Git 和理论笔记均存在。
- S：源码存在于 Git `main/` 和旧工作区同步目录。
- 旧 Q/S/SM/QZ：源码和部分 MAT 产物存在，但未逐一重新运行，本报告只记录可从源码确认的算法结构。
- `copyP_fast_reid.m`：只在历史 session 中出现过文件名，当前磁盘未找到；因此不把它伪装成可恢复源码或有效指标版本。
- 任何旧版性能数字在重新运行并持久化 metrics 前，不进入现行版本性能排名。
