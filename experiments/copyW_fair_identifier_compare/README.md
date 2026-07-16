# copyW：同一 plant/controller 下的辨识器公平比较

`copyW_fair_identifier_compare` 固定 copyV 的 plant、离线数据、参考、噪声、控制器、约束和随机种子，只替换辨识器，比较：

1. `copyP_svd_ols`：输出 PCA/SVD 正交子空间；
2. `main_qr_ivr`：main 风格 IVR 后 QR 正交化；
3. `copyO_oblique_ivr`：自由斜投影 IVR 双基；
4. `copyV_control_aware_ivr`：固定 tracked axes 后在补空间迭代 IVR。

这不是按字母顺序继承 copyV 的新控制算法，而是 identifier-only ablation。关键判据除了跟踪误差，还包括输出预测 RMSE、重构残差、tracked coverage、QP 成功率和 fallback。

运行：

```matlab
run('experiments/copyW_fair_identifier_compare/copyW_fair_identifier_compare.m')
```

持久化表：`results/copyW_fair_identifier_compare_metrics.csv`。

当前固定种子结果表明，自由预测性子空间的预测 RMSE 并不自动转化为控制能力；`copyO_oblique_ivr` 的 tracked coverage 误差最大且 QP 成功率只有 0.13，而 control-aware IVR 覆盖误差为 0、QP 成功率为 1。不要把不同 identifier 的 MAE 差异单独解释成算法普遍优劣。生成的 `.mat/.png` 不提交 Git。
