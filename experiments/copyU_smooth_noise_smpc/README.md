# copyU — 平滑时变噪声 PredVARX-SMPC

## 目的

从已跑通的 `copyT_process_lv_smpc` 独立复制，保持 `main/` 与 copyT 不变，测试在线协方差估计和机会约束能否跟随平滑增长/下降的噪声强度。

## 噪声定义

噪声样本仍是高斯白噪声，只有瞬时标准差使用 raised-cosine 包络：

$$
w_k=\sigma_w(k)\xi_k,\qquad v_k=\sigma_e(k)\eta_k,
$$

$$
\xi_k\sim\mathcal N(0,I_n),\qquad \eta_k\sim\mathcal N(0,I_p).
$$

当前设置：

- 周期：400 步，1200 步闭环经历 3 个周期；
- $\sigma_w(k)\in[0.020,0.090]$；
- $\sigma_e(k)\in[0.025,0.100]$；
- 测量噪声相位偏移：$\pi/3$；
- 离线辨识数据仍使用 copyT 的固定噪声尺度，专门测试训练—运行噪声分布变化。

这属于对 Mo–Qin 平稳高斯噪声模型的**平滑非平稳扩展**，不是原文 Algorithm 1 的噪声要求。

## 最初从 copyT 派生时的实验差异

1. 闭环过程噪声和测量噪声改用平滑标准差包络；
2. MAT 文件保存 `sigma_w_profile` 和 `sigma_e_profile`；
3. 噪声诊断图同时显示：
   - 两条指定包络；
   - 两条随机实现 RMS；
   - 两条在线滑窗估计；
4. 控制感知子空间、VARX 与 Boole 风险分配沿用 copyT 的总体结构。

> 后续审查发现 copyT/copyU 共用的旧 QP 构造存在输入边界长度、中心化和代价尺度问题。下面的修复仅落在独立 copyU 中，因此修复后的 copyU QP 已不再与当前 copyT 实现完全相同。

## 2026-07-11 控制器一致性修复

审查发现旧版 `lb/ub` 只有 $N=18$ 项，而决策向量有 $Nn_u=54$ 项，导致只约束前6个预测时刻。本版已修复并重新运行：

- 标量/逐通道输入界均扩展为完整 $Nn_u$ 维；
- 目标与机会约束统一采用 $U-U_0$；
- 输入代价采用 $(U-U_0)^T\bar R_u(U-U_0)$；
- `quadprog` 使用 $H_{qp}=2\widetilde H, f_{qp}=2\widetilde f$；
- `costJ` 记录严格的中心化原始 MPC 代价；
- 新增 `tests/test_copyU_controller_consistency.m` 回归测试。

修复后的主要指标见 `results/copyU_smooth_noise_smpc_metrics.txt`。旧结果视为修复前历史结果。

## 运行

```matlab
run('experiments/copyU_smooth_noise_smpc/tests/test_copyU_smooth_noise_profile.m')
run('experiments/copyU_smooth_noise_smpc/tests/test_copyU_controller_consistency.m')
run('experiments/copyU_smooth_noise_smpc/copyU_smooth_noise_smpc.m')
```

## 验收标准

- 包络测试通过；
- QP success rate = 1；
- fallback count = 0；
- `max(AU-b) <= 10^{-9}`（数值容差内）；
- 两个 tracked output 的经验越界率为 0；
- 图中包络平滑、随机 RMS 围绕包络抖动、40 步估计带延迟跟随。

## 文件

```text
copyU_smooth_noise_smpc/
├── copyU_smooth_noise_smpc.m
├── smooth_noise_profile.m
├── control_aware_subspace_varx.m
├── centered_smpc_step.m
├── tests/test_copyU_smooth_noise_profile.m
├── tests/test_copyU_controller_consistency.m
└── results/
```
