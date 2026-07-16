# PredVARX–MPC Git 版本演化、差异与指标总表

> 审计日期：2026-07-16
>
> 仓库：`git_repo_2026-07-09`
>
> 口径：优先采用 Git 中持久化的 README、metrics TXT/CSV 和实际重跑输出；`.mat/.png` 只作运行证据，不纳入 Git。不同 plant、维数、噪声、时域或控制器的行不可直接按 MAE 排名。

> 旧工作区还存在另一套重用 Q/S/Z 等字母的实验。其版本考古和与现行 Git 的对照见 [`LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md`](LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md)。本文以下的 O–Z 默认指 2026-07-09 Git 仓库中的现行目录名。

## 1. 先给结论

版本演化不是一条“每个字母都更好”的单线，而是四条问题线：

1. **历史基线与斜投影身份**：`main → copyO → copyR`；
2. **修正 MPC 和控制覆盖**：`copyP → copyQ → copyT → copyU`；
3. **同一高维 plant 上的辨识器演化**：`copyU → copyV → copyX → copyY / copyZ`；
4. **公平辨识器消融**：`copyW`。

最重要的事实：

- `copyX` 只中心化，不白化；它直接使用 $Y^c$ 的 tracked 正交补做 IVR。
- `copyY` 与 `copyX` 数值等价，只把“无白化 + 直接方向替换”显式标记出来，不是新算法。
- `copyZ` 才执行白化与直接反白化，但它是满秩、一阶消融，不是完整 Algorithm 1；在同一 copyX plant 上控制显著退化。
- `copyV`、`copyW` 原先存在于工作区和历史记忆但未进 Git，本次已重跑、持久化并提交。

## 2. Git 演化节点

| Git | 版本/作用 | 核心变化 |
|---|---|---|
| `f0cefcf` | `main` 原始基线 | 原始 `copyS_nosoft` 与 `predvarx_identify`；历史对照，不作为当前正确控制基准。 |
| `e84ffb4` | `copyO` | 恢复斜投影双基 $R^TP=I$，但保留旧控制逻辑。 |
| `ce79ce3` | `copyP` | 中心化、绝对输出预测、Boole 双侧机会约束；正交全阶控制验证。 |
| `c5b077a` | `copyQ` | 低阶且强制保留 $e_1,e_2$ tracked axes。 |
| `a3efee3` | `copyR` | 论文归一化/反归一化式的白化 PredVAR 对照；无控制覆盖。 |
| `4997a36`–`50c8f27` | `copyT` | 高维 process-like plant、在线噪声尺度更新和约束激活诊断。 |
| `a486770`–`716bd96` | `copyU` | 平滑时变高斯噪声；修正中心化 SMPC 一致性。 |
| `2e815b1` | `copyV` | 在 tracked 补空间用迭代 IVR 替代单次 SVD；本次补入 Git。 |
| `e069a6e` | `copyW` | 同一 plant/controller 下四种 identifier 公平消融；本次补入 Git。 |
| `82e638e` | `copyX` | 基于 copyV 的 $P$，加入正则化斜读取器 $R_\alpha$，$\alpha=0.02$。 |
| `80d4fe3` | `copyY` | 显式标注 copyX 的真实路径：无白化、原坐标直接更新；数值与 copyX 一致。 |
| `f652dca`、`09b0322` | `copyZ` | copyX plant 上的满秩 $s=1$ 白化—IVR—直接反白化消融，并校正“完整 Algorithm 1”过强标签。 |
| `f2384d7` | copyX 边界修复 | 修复 $\ell=q$ 时 `iter` 未定义；当前 $\ell=5,q=2$ 指标不变。 |

## 3. 各版本算法差异

| 版本 | 数据预处理 | 子空间/读取器 | 控制覆盖 | MPC/噪声特点 | 正确定位 |
|---|---|---|---|---|---|
| `main` | 历史实现 | IVR 后 QR，近似令 $R=P$ | 不保证 | 旧差分 MPC、旧风险口径 | 只作历史基线 |
| `copyO` | 旧流程中的归一化 IVR | SVD 对齐双基，满足 $R^TP=I$ | 不保证 | 旧控制逻辑仍存在 | 验证斜投影代数，不验证控制优越性 |
| `copyP` | $Y^c,U^c$，无白化 | SVD/PCA 正交子空间，$R=P$，$\ell=6$ | 全阶 plant 状态映射下可控 | 绝对预测 + Boole 联合机会约束 | 控制逻辑正确性基线 |
| `copyQ` | $Y^c,U^c$，无白化 | $P=[e_1,e_2,Q_\perp]$，$R=P$ | $PP^Te_i=e_i$ | 与 copyP 相同 | 低阶控制感知正交版本 |
| `copyR` | EVD 白化 | 直接反白化 $P,R,\bar P,\bar R$ | 不保证 | copyP 型 SMPC | p=15 的论文白化对照，不与 copyX 直接量化比较 |
| `copyT` | 中心化，无白化 | tracked axes + 补空间 SVD，$R=P$ | 精确 | 固定噪声 plant；在线更新 $\Sigma_\varepsilon,\Sigma_{obs}$ | 高维工业过程基线 |
| `copyU` | 中心化，无白化 | 同 copyT | 精确 | 平滑时变噪声 | copyT 的时变噪声压力测试 |
| `copyV` | 中心化，无白化 | tracked axes + 补空间迭代 IVR，$R=P$ | 精确 | 与 copyU plant/controller 相同 | 正交控制感知 IVR 基线 |
| `copyW` | 与 copyV 相同 | 四种 identifier，仅替换辨识器 | 依 identifier 而异 | controller/plant 固定 | 公平消融，不是单一新控制器 |
| `copyX` | **中心化，无白化** | $P=[E_c,N_\perp V]$；$R_\alpha=P+0.02(R_{full}-P)$ | $PR^TE_c=E_c$ | 在线只更新协方差尺度 | 控制感知、轻微斜投影工程原型 |
| `copyY` | **中心化，无白化** | 与 copyX 同一数值更新，只显式命名 | 精确 | 与 copyX 完全相同 | copyX 的“事实标签版本”，不是新算法 |
| `copyZ` | **EVD 白化** | 白化空间 IVR，直接反白化；无 tracked 约束、无 $R_\alpha$ | 不保证 | 同 copyX plant/controller | 满秩 $s=1$ 白化消融，不是完整 Algorithm 1 |

## 4. 可直接读取的闭环指标

### 4.1 早期 p=15 实验（只在各自设置内解释）

| 版本 | MAE $(y_1,y_2)$ | RMSE $(y_1,y_2)$ | QP/越界 | 其他关键量 | 结论 |
|---|---:|---:|---|---|---|
| `copyO` | README 只记录 $y_2$: C=3.2488，Oblique=6.3173 | — | Oblique 越界 26 次 | 双基测试通过 | 斜投影代数正确，但旧控制逻辑失败 |
| `copyP` | (0.0931, 0.0867) | (0.1156, 0.1082) | QP 1200/1200；零越界 | max chance residual $-2.296\times10^{-1}$ | 修正后的全阶控制基线有效 |
| `copyQ` | (0.1142, 0.0930) | (0.1424, 0.1176) | QP 1200/1200；零越界 | coverage=0；重构残差 0.3023 | 低阶但保留控制方向 |
| `copyR` | 未形成有效闭环 MAE | — | 第 1 步 QP 不可行 | coverage error 约 (2.34, 10.82) | 预测子空间不等于控制覆盖 |

### 4.2 p=30 process-like 演化链

| 版本 | MAE $(y_1,y_2)$ | RMSE $(y_1,y_2)$ | 平均 $J$ | QP 成功率 / fallback | 重构残差 | coverage/dual | 说明 |
|---|---:|---:|---:|---:|---:|---|---|
| `copyT` | (0.0812, 0.0609) | (0.1047, 0.0763) | 678.60 | 100% / 0 | 0.1103 | tracked=0 | 固定噪声设置，代价口径/参考条件与 U–Z 不宜直接排名 |
| `copyU` | (0.0959, 0.0959) | (0.1425, 0.1385) | 110.36 | 100% / 0 | 0.1103 | tracked=0 | 平滑时变噪声，单次 SVD 补空间 |
| `copyV` | (0.0963, 0.0953) | (0.1429, 0.1378) | 107.91 | 100% / 0 | 0.1190 | tracked=0 | 迭代 IVR，较 U 只小幅改变结果 |
| `copyX` | (0.0973, 0.0950) | (0.1444, 0.1375) | 109.16 | 100% / 0 | 0.1191 | tracked $1.03\times10^{-14}$；dual $1.04\times10^{-14}$ | 无白化，$\alpha=0.02$ 轻微斜读取器 |
| `copyY` | (0.0973, 0.0950) | (0.1444, 0.1375) | 109.16 | 100% / 0 | 0.1191 | 与 X 相同 | 与 copyX 逐位一致，只有标签更明确 |
| `copyZ` | (0.5172, 0.7246) | (0.7024, 0.9747) | 786.36 | 39.67% / 724 | 0.2622 | coverage 1.6327；dual $8.01\times10^{-15}$ | 白化几何正确，但没有控制覆盖，闭环显著退化 |

共同安全指标：T–Z 的持久化结果中，上界和绝对值越界率均为 0；但 copyZ 大量 fallback，不能把零越界解释为成功的概率校准。

## 5. copyW 公平 identifier 消融

| Identifier | MAE $(y_1,y_2)$ | prediction RMSE | reconstruction | tracked coverage error | 平均 $J$ | QP / fallback |
|---|---:|---:|---:|---:|---:|---:|
| PCA/SVD orthogonal | (0.0770, 0.0788) | 0.0948 | 0.0610 | 1.2995 | 74.75 | 100% / 0 |
| main QR-IVR | (0.0770, 0.0843) | 0.0963 | 0.0822 | 1.2989 | 78.31 | 100% / 0 |
| copyO free oblique IVR | (0.9685, 1.1231) | 0.2842 | 0.3524 | 1.9038 | 2103.07 | 13% / 1044 |
| copyV control-aware IVR | (0.0917, 0.0947) | 0.0995 | 0.1190 | 0 | 101.55 | 100% / 0 |

这张表显示：在该 plant 上，最低输出预测 RMSE 并不自动等于最强控制覆盖；control-aware 版本牺牲少量预测/重构指标换取 tracked axes 的严格表示。copyO 的失败也不能仅归因于“斜投影”三个字，而是自由预测子空间与控制目标严重错位。

## 6. 最终总表（建议以后优先引用）

| 版本 | 是否白化 | 是否控制感知 | $R=P$？ | 与谁公平可比 | 核心结果 | 最终判断 |
|---|---|---|---|---|---|---|
| main | 历史实现标签不宜继续外推 | 否 | QR 后近似是 | 仅历史 | 旧 MPC/风险口径 | 存档，不作为结论基线 |
| O | 是/论文式归一化路径，但有后验对齐 | 否 | 否 | main/O 自身 | 双基正确，控制差 | 代数诊断版 |
| P | 否 | 否；但全阶 | 是 | P/Q/R 的 p=15 plant | MAE 0.093/0.087，QP 100% | 正确控制基线 |
| Q | 否 | 是 | 是 | P/Q/R | MAE 0.114/0.093，coverage 0 | 低阶控制感知基线 |
| R | 是 | 否 | 否 | P/Q/R（同 plant），不与 X/Z | 首步 QP 不可行 | 论文预测基线不自动可控 |
| T | 否 | 是 | 是 | 主要自身 | MAE 0.081/0.061 | 固定噪声 process 基线 |
| U | 否 | 是 | 是 | V/X/Y/Z | MAE 0.096/0.096 | 时变噪声 SVD 基线 |
| V | 否 | 是 | 是 | U/X/Y/Z | MAE 0.096/0.095 | 正交 IVR 基线 |
| W | 依 identifier | 依 identifier | 依 identifier | W 内部四项 | coverage 对 QP 影响显著 | 公平消融证据 |
| X | **否** | **是** | 否，$\alpha=0.02$ | V/Y/Z | MAE 0.097/0.095，QP 100% | 可运行的轻微斜投影原型；未优于 V |
| Y | **否** | **是** | 同 X | X | 与 X 逐位一致 | 事实标签版本，不是新算法 |
| Z | **是** | **否** | 否 | X/Z 最公平 | MAE 0.517/0.725，QP 39.67% | 白化几何正确但控制覆盖失败 |

## 7. Git 补充记录

本次把原先只存在于工作区/记忆、未纳入 Git 的版本补齐：

- `2e815b1 Add copyV iterative IVR experiment`
- `e069a6e Add copyW fair identifier comparison`

未提交 `.mat/.png`。用户原有未提交修改（`copyT_process_lv_smpc.m`、其 metrics、`main/copyS_nosoft.m`）未纳入这些提交。

## 8. 一句话研究结论

> 目前最可靠的研究链不是“白化一定更好”或“斜投影一定更好”，而是：**预测子空间的双基几何、tracked control coverage、坐标一致的 VARX 和可行的 SMPC 必须同时成立。** copyX 证明轻微斜读取器可以在控制覆盖不丢失时运行；copyZ 和 copyW 则证明，仅有白化/双基/预测性并不足以保证闭环控制能力。
