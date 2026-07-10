# copyP：居中全阶子空间 SMPC

> **目标**：修正 copyO 中的坐标不一致、错误风险分位数与不一致的差分代价，构造一个可验证的控制实验。

## 控制模型

离线数据先分别中心化：

$$
y_k^c=y_k-\bar y,\qquad u_k^c=u_k-\bar u.
$$

由于本实验取 $\ell=n=6$，直接在中心化输出数据上取正交全阶子空间：

$$
z_k=P^T(y_k-\bar y),\qquad P^TP=I,
$$

并回归：

$$
z_{k+1}=Az_k+B(u_k-\bar u)+\varepsilon_k.
$$

这里 $R=P$ 是**明确的正交全阶特例**，不是把它伪称为一般斜投影 PredVAR。这样控制坐标与辨识坐标一致，且输出前两个通道在完整 $6$ 维状态空间中可控。

## 正确的绝对跟踪 MPC

使用：

$$
\mu_y(j)=\bar y+PA^jz_k+G_j(\mathbf U-\mathbf 1\otimes\bar u),
$$

$$
J=\sum_{j=1}^{N}\|\mu_y(j)-r_{k+j}\|_Q^2+
\sum_{j=0}^{N-1}\|u_{k+j}\|_{R_u}^2.
$$

没有把绝对输出约束塞进差分预测，也没有使用无法与 chance constraint 对齐的旧 `J_const` 记录方式。

## 联合机会约束

目标是双输出、双侧、整个预测域的联合安全概率：

$$
\Pr\{H y_{k+j}\le h,\ \forall j\}\ge1-\alpha_{\rm joint}.
$$

用 Boole 分配：

$$
\epsilon=\frac{\alpha_{\rm joint}}{2n_cN},
\qquad z_\epsilon=\Phi^{-1}(1-\epsilon)>0.
$$

本实验 $\alpha_{\rm joint}=0.20$、$n_c=2$、$N=15$，故：

$$
z_\epsilon=2.7131.
$$

每个上/下界约束用：

$$
h_i^T\mu_y(j)+z_\epsilon\sqrt{h_i^T\Sigma_y(j)h_i}\le h_i.
$$

## 实测结果

MATLAB R2024a 完整运行 1200 步：

| 指标 | $y_1$ | $y_2$ |
|:---|---:|---:|
| MAE | 0.0931 | 0.0867 |
| RMSE | 0.1156 | 0.1082 |
| 与参考相关系数 | 0.9867 | 0.9881 |
| 上界越界 | 0 | 0 |
| 双侧越界 | 0 | 0 |

所有 $1200/1200$ 个 QP 均成功求解；最大 chance-constraint 残差为 $-2.296\times10^{-1}<0$。

## 文件

- `copyP_centered_smpc.m`：主实验
- `control_ready_subspace_varx.m`：中心化全阶子空间 VARX
- `centered_smpc_step.m`：绝对预测、Boole 风险分配 QP
- `tests/test_centered_smpc_step.m`：中心化与 Boole 机会约束单元测试。
- `results/`：MAT 数据、PNG 图和运行日志快照。
