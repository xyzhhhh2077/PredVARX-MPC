# copyO_oblique：斜投影双基 PredVARX-MPC 实验

> **目的**：在不改动 `main/` 基准的前提下，修正原实现把 `R_hat=P_hat` 的 QR 正交化简化，恢复 PredVAR 所需的双基关系。

## 文件

- `predvarx_identify_oblique.m`：OLS 输入剥离 + IVR + SVD 双基对齐 + VARX 回归。
- `copyO_oblique.m`：在 `main/copyS_nosoft.m` 基准上建立的闭环实验副本。
- `tests/test_predvarx_identify_oblique.m`：回归测试。
- `results/`：MAT 数据、PNG 图和运行日志快照。

## 核心构造

IVR 反白化后得到原始对：

$$
P_{\rm raw}=U D^{1/2}C^*,\qquad R_{\rm raw}=U D^{-1/2}C^*.
$$

令

$$
S=R_{\rm raw}^T P_{\rm raw}=U_s\Sigma_sV_s^T.
$$

使用：

$$
P=P_{\rm raw}V_s\Sigma_s^{-1},\qquad R=R_{\rm raw}U_s,
$$

从而保证：

$$
R^TP=I_\ell.
$$

再以 $[P,\bar P]$ 的完整方阵基的逆转置构造整个对偶基 $[R,\bar R]$，因此同时满足：

$$
R^TP=I,\quad R^T\bar P=0,\quad \bar R^TP=0,\quad \bar R^T\bar P=I.
$$

静态残差必须用 $\bar R^T$ 提取：

$$
\bar e_k=\bar R^T(y_k-Pz_k),\qquad z_k=R^Ty_k.
$$

不能沿用正交版本中的 `Pbar'*(y-P*x)`；当 $P,\bar P$ 不正交时那不是正确的对偶坐标。

## 已验证结果

MATLAB R2024a 运行：

```text
PASS: oblique PredVARX dual-basis invariants hold.
```

完整 2000 步闭环运行完成并生成：

- `copyO_oblique_data.mat`
- `copyO_oblique_fig.png`

本次固定随机种子下：C 的 $y_2$ MAE 为 `3.248796`，Oblique-SMPC 为 `6.317314`；C 越界 `0` 次，Oblique-SMPC 越界 `26` 次。

> 这说明双基理论约束已通过，但**不代表控制表现改善**。当前 SMPC 仍有独立问题：`alpha_cc=0.45` 使正态分位数为负，机会约束实际放松；协方差滚动携带和 C/Z 代价记录口径也仍需单独修正。
