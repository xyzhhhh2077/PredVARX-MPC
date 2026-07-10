# copyQ：低阶控制感知子空间 SMPC

> **目的**：在 copyP 的中心化、绝对预测与正确 chance tightening 基础上，恢复低阶模型，同时保证受跟踪输出不会被降维子空间丢掉。

## 为什么不能只按可预测性选低阶 $P$

copyO 使用 $\ell=2$ 的 IVR 子空间。它可以是“可预测”的，却未必包含 $y_1,y_2$ 的受控方向，因此 MPC 可能无法改变它真正被罚的输出。

copyQ 取 $\ell=4$，并强制：

$$
e_1,e_2\in\operatorname{span}(P).
$$

具体构造：

$$
P=[e_1,e_2,Q_\perp],
$$

其中 $Q_\perp$ 是剩余输出空间上中心化离线数据的前两个主动态方向。故：

$$
P^TP=I,
\qquad PP^Te_i=e_i,\ i=1,2.
$$

这使受控输出轴的投影误差严格为零。

## 其他控制逻辑

与 copyP 相同：

- 辨识与在线均使用中心化 $z=P^T(y-\bar y)$；
- 使用绝对输出预测，而不是差分目标拼接；
- 用双侧 Boole 分配的联合机会约束；
- $\alpha_{\rm joint}=0.20$、$N=15$、正的正态分位数 $2.7131$；
- 输入范围 $[-4,4]$，输出安全边界 $y_{\max}=2.5$。

## 实测

MATLAB R2024a，1200 步：

| 指标 | $y_1$ | $y_2$ |
|:---|---:|---:|
| MAE | 0.1142 | 0.0930 |
| RMSE | 0.1424 | 0.1176 |
| 与参考相关性 | 0.9840 | 0.9872 |
| 上界越界 | 0 | 0 |
| 双侧越界 | 0 | 0 |

其他验证：

```text
QP 成功：1200/1200
最大 chance residual：-2.531e-13 <= 0
输入饱和次数：0
受控输出投影误差：0
降维重构残差：0.3023
```

> $0.3023$ 是四维模型丢弃未受控输出方向后的正常信息损失；它不影响 $y_1,y_2$ 的表示，因为二者被强制完整保留。

## 文件

- `control_aware_subspace_varx.m`
- `copyQ_control_aware.m`
- `tests/test_control_aware_subspace.m`
- `copyQ_control_aware_data.mat`
- `copyQ_control_aware_fig.png`
