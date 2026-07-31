# copyAU：软偏好学习最终输出

使用 copyAR 的同一离线 `(u,y)` 数据，不新增训练集。输出权重按
`w_i=exp(-0.18(i-1))` 逐渐减小，学习目标把归一化偏好矩阵与有限时域输入权威矩阵加权，取最大两个方向。

在学习到的软锚方向下，自由块使用固定内点 `mu=0.10` 的 CRTE 度量，并使用修正后的论文归一化
`Ntr(A)=A/max(abs(trace(A))/d,10e-6)`。MAT 文件保存了 `A_T`、三个实际归一化矩阵及辨识统计，测试会直接重算 `Ntr(A_T)`。

最终控制与评价对象是 `s=Etask'*y`，不是旧的 `y1,y2`。参考、Q和机会约束都变换到该二维方向。

本副本只提供1200步单种子实验，不证明方向全局最优、稳定性或递归可行性。

## 实跑结果

- 离线训练样本：1500；新增训练样本：0
- 闭环长度：1200步
- MAE：`[1.0305, 1.2664]`
- RMSE：`[1.1957, 1.3417]`
- QP成功率：`100%`；fallback：`0`
- 最大 QP 不等式余量 `max(AU-b)`：`-1.7672`

方向贡献主要集中在前两个输出，但学习输出几乎未跟踪分段参考。修正 `mu` 和 `Ntr` 后，当前“软偏好+输入权威”构造仍只证明 QP 路径可行，没有实现所需任务跟踪；该结果作为负面实验保留，不作为性能改进。

![copyAU软偏好1200步结果](results/copyAU_soft_preference_output_fig.png)

数据见 [`results/copyAU_soft_preference_output_data.mat`](results/copyAU_soft_preference_output_data.mat)，完整指标见 [`results/copyAU_soft_preference_output_metrics.txt`](results/copyAU_soft_preference_output_metrics.txt)。
