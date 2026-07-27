# copyAS — 学习式任务锚 PredVARX-SMPC 试验

## 身份

这是基于 CRTE 双基几何的新扩展，**不是原文算法**。原文要求人工指定标准基 `E_c`；本副本从完整输出任务参考中学习任务子空间 `E_T`，不接收 tracked-output 索引。

## 方法

1. 完整输出一步回归得到预测内容 `S_yu^(y)` 与残差噪声 proxy。
2. 全输出任务参考得到 `M_task`。
3. 全输出 VARX 的 Markov 响应得到 `W_y`。
4. 组成 `M_E = task + predict + reach - noise`（Frobenius 归一）。
5. 解 `M_E E_T = C_y E_T Lambda_E`。
6. 构造配对任务锚 `P_T=E_T`, `R_T=C_y E_T`。
7. 在 `null(R_T')` 中学习自由方向，建立完整双基、重拟合 VARX。
8. 用同类 centered SMPC 闭环验证；物理约束仍作用于原始输出。

## 运行

```matlab
cd('experiments/copyAS_learned_task_anchor_smpc')
run('tests/test_learned_task_anchor_varx.m')
copyAS_learned_task_anchor_smpc
```

## 验证边界

- 单元测试验证配对双基、任务子空间保持、PSD 协方差和无 tracked 索引合同。
- smoke 验证一次闭环轨迹、QP 成功率、fallback、原输出跟踪和约束。
- 不证明闭环稳定性、递归可行性、全局最优或优于固定 `E_c`。
- 真实 `Sigma_n` 只用于 plant 生成，从不进入学习算法。
