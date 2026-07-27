# copyAV：硬偏好选择最终输出

沿用与 copyAU 相同的递减权重，但直接选择权重最大的两个原始输出轴并硬锁定。本设置下最终输出严格为 `y1,y2`，其他输出不能旋转进任务子空间。

复用 copyAR 同一离线训练集和已验证模型，不新增训练数据；运行1200步同种子闭环，保存独立PNG和MAT。

该结果验证当前权重排序下的硬选择，不证明这两个物理输出对所有任务都最优。

## 实跑结果

- 离线训练样本：1500；新增训练样本：0
- 闭环长度：1200步
- 选择输出：`[y1, y2]`
- MAE：`[0.05973, 0.06185]`
- RMSE：`[0.08141, 0.08669]`
- QP成功率：`100%`；fallback：`0`

硬偏好阻止其他输出旋转进任务子空间，因此最终控制对象仍是原始任务输出。当前单种子结果恢复了与copyAR相近的跟踪水平，但不证明该选择对其他任务或权重排序普遍最优。

![copyAV硬偏好1200步结果](results/copyAV_hard_preference_output_fig.png)

数据见 [`results/copyAV_hard_preference_output_data.mat`](results/copyAV_hard_preference_output_data.mat)，完整指标见 [`results/copyAV_hard_preference_output_metrics.txt`](results/copyAV_hard_preference_output_metrics.txt)。
