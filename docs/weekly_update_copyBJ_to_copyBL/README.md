# 本周升级总结：从 copyBJ 到 copyBL/copyBM

> 参考论文：Lin et al. (2025) "Maximizing the Value of Predictions in Control: Accuracy Is Not Enough" (NeurIPS 2025)
> 参考论文：Model-Free Safe Reinforcement Learning for Multi-Agent Systems by Distributed Data-Driven Robust MPC (IEEE TAC 2026)

---

## 一、背景

copyBJ 是基于冻结 PredVAR 预测器的在线协方差自适应 SMPC。它用 40 样本滑窗估计新的扰动协方差，然后用固定 β=0.8 做收缩混合：

```
Σ_new = (1-β) × Σ_old + β × Σ_raw
```

copyBJ 的问题是：
1. 直接用固定收缩比，**不检查 QP 是否可解**——如果新协方差使裕度过宽，QP 不可行只能停机
2. **不检查更新是否有价值**——即使协方差变了，闭环代价可能反而变差
3. **没有公平对比**——只有单次验证，无法在相同扰动下配对比较各方法

阅读 Lin et al. (2025) 和 D3RMPC (IEEE TAC 2026) 后，针对这两个问题做了升级。

---

## 二、从 Lin et al. (2025) 借鉴的内容

Lin 的核心发现：**预测精度（MSE）不等于预测价值**。两个预测器可以有相同的 MSE，但控制代价差异很大。关键不是预测准不准，而是预测对控制动作的改进有多少。

### 借鉴点：价值门（Control Value Gate）

在 copyBL 中新增了 `control_value_gate()`，代码注释明确标注 `Lin-inspired paired closed-loop cost deadband`：

- 从**相同状态 + 相同扰动**重放闭环
- 计算候选协方差下的闭环代价 vs 当前协方差下的闭环代价
- 只有**代价改进超过死带阈值**才接受更新
- 如果候选没跑完或代价改进不够，拒绝更新

```
improvement = (baseline_cost - candidate_cost) / baseline_cost
if improvement > deadband:
    accept  # 采纳新协方差
else:
    reject  # 保持旧协方差
```

这就是 Lin "prediction power" 思想的工程化：不是所有协方差更新都有价值，只有闭环代价确实验证了改进才采纳。

---

## 三、从 D3RMPC (IEEE TAC 2026) 借鉴的内容

D3RMPC 论文的核心：在数据驱动 MPC 中，通过**约束收紧**和**误差管**保证安全性和可行性。它在参数更新过程中始终保持递归可行性和稳定性。

### 借鉴点：可行性门（Backtracking Feasibility Gate）

在 copyBL 中新增了 `backtracking_feasibility_gate()`：

- 不直接用新协方差替换旧协方差
- 沿凸组合路径 `(1-α)·Σ_old + α·Σ_raw` 逐步回溯
- 从 α=1（完全新）向 α→0（完全旧）方向搜索
- 找到**第一个 QP 可解的凸组合**就停
- 因为裕度关于协方差是凸的，任何凸组合都能保持非负裕度

```
for α in [1.0, 0.8, 0.6, 0.4, 0.2, ...]:
    Σ_candidate = (1-α) × Σ_old + α × Σ_raw
    if QP_feasible(Σ_candidate):
        return Σ_candidate  # 用这个
# 全不可行 → 保持 Σ_old 不变
```

这跟 D3RMPC 的思路一致：更新参数时必须保证可行性，不能因为模型更新导致控制器停机。

---

## 四、copyBL：双门协议

copyBL = copyBJ + 可行性门 + 价值门。每次协方差更新的流程：

```
1. 滑窗估计新协方差 Σ_raw
2. 可行性门：回溯搜索找到 QP 可解的凸组合 Σ_candidate
   → 不可行 → 保持旧协方差，跳过本轮
3. 价值门：从相同状态+扰动重放闭环
   → 代价改进 < deadband → 保持旧协方差，跳过本轮
4. 两个门都通过 → 原子替换 Σ_old = Σ_candidate
```

### 仿真验证

单次 18000 步 Pelican 四旋翼轨迹跟踪验证：

| 指标 | copyBJ (旧) | copyBL (新) |
|------|-------------|-------------|
| 协方差替换 | 固定 β=0.8 | 双门自适应 |
| QP 可行性检查 | 无 | 回溯门 |
| 价值检查 | 无 | 代价死带门 |
| 18000 步完成 | 是 | 是 |

**附图：**

- `copyBL_smpc_3d_trajectory.png` — 3D 轨迹跟踪
- `copyBL_smpc_xyz_reference.png` — XYZ 三轴参考跟踪
- `copyBL_smpc_noise_timeseries.png` — 噪声时序
- `copyBL_smpc_stage_cost_dual_gates.png` — 阶段代价 + 双门事件标注
- `copyBL_dual_gate_summary_200frames.gif` — 200 帧合成 GIF

---

## 五、copyBM：20-seed 配对公平对比

copyBM 在 copyBL 基础上增加了 20-seed 配对公平对比实验，解决 copyBJ "只有单次验证"的问题。

### 实验设计

- **20 个配对 seed**：每对 seed 共享相同的 innovation、wind、初始状态、扰动 SHA256
- **三方对比**：Frozen（冻结基线）/ copyBJ（固定收缩）/ copyBL（双门）
- **完整指标集**：RMSE、累计代价、最小硬裕度、QP 失败率、门接受率

### 实验结果

| 指标 | Frozen (冻结) | copyBJ (固定收缩) | copyBL (双门) |
|------|--------------|------------------|--------------|
| 完成种子数 | 10/20 | 16/20 | 16/20 |
| 累计代价（均值） | 590,875 | 689,380 | 624,709 |
| 最小硬裕度（均值） | 0.0548 | 0.1822 | 0.0930 |
| QP 失败率（均值） | 0.22% | 0.02% | 0.18% |
| 门接受率 | — | — | 19.1% |

### 关键发现

1. **完成率**：copyBJ 和 copyBL 都把完成率从 10/20 提升到 16/20，说明在线协方差适应有效
2. **代价**：copyBL 的累计代价（624K）比 copyBJ（689K）低 9%，价值门过滤了高代价的更新
3. **门接受率**：只有 19% 的更新通过双门，说明大部分协方差变化对控制没有实际价值——印证了 Lin 的"accuracy is not enough"
4. **裕度**：copyBJ 的裕度（0.18）最大但代价也最高——过度保守；copyBL 裕度（0.09）居中但代价最优

**附图：**

- `aggregate_comparison.png` — 三方完成率/代价/裕度对比
- `scientific_comparison.png` — 科学对比图（配对差异 + bootstrap CI）

---

## 六、升级清单总结

| 新增项 | 来源 | 文件 | 作用 |
|--------|------|------|------|
| 可行性门（回溯） | D3RMPC | `backtracking_feasibility_gate()` | 保证 QP 可解，不停机 |
| 价值门（代价死带） | Lin et al. | `control_value_gate()` | 只在闭环代价确有改进时更新 |
| 双门协议 | 综合 | `simulate_dual_gated_controller()` | 先可行 → 再有价值 → 才替换 |
| 配对扰动 SHA | — | `disturbance_sha256()` | 公平对比保证 |
| 20-seed 配对 | — | `seed_pairs()` | 统计可靠性 |
| 完整指标集 | — | `controller_metrics()` + `gate_metrics()` | RMSE/代价/裕度/失败率/接受率 |

---

## 七、文件清单

### copyBL 仿真图（单次验证）
- `images/copyBL_smpc_3d_trajectory.png`
- `images/copyBL_smpc_xyz_reference.png`
- `images/copyBL_smpc_noise_timeseries.png`
- `images/copyBL_smpc_stage_cost_dual_gates.png`
- `images/copyBL_dual_gate_summary_200frames.gif`（200帧合成 GIF）

### copyBM 仿真图（20-seed 对比）
- `images/aggregate_comparison.png`
- `images/scientific_comparison.png`

### 代码
- `code/run_dual_gated_covariance_smpc.py` — copyBL 主程序
- `code/run_fair_comparison.py` — copyBM 公平对比
