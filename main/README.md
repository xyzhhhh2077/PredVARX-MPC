# main：原始基线

## 文件

- `copyS_nosoft.m`：原始 C 与 Zfix/SMPC-only 对比脚本。
- `predvarx_identify.m`：原始 OLS 输入剥离 + IVR + QR 正交化 + VARX 辨识器。
- `results/`：本次运行保存的 MAT、图和日志快照。

## 基线的作用

`main/` 是对照版本，不应用 copyO/copyP/copyQ 的算法重写。唯一增加的是更完整的 MAT 导出字段：真实系统、离线数据、在线模型快照和版本标记，以支持后续诊断。

## 已知与后续副本的差异

- 采用 QR 后 `R=P` 的正交化实现；
- 在线控制和辨识中心化坐标并不严格一致；
- `alpha_cc=0.45` 导致正态分位数为负，机会约束实际放松；
- MPC 差分目标和绝对输出机会约束并不在统一预测口径。

这些问题不在本目录修补，分别在 `experiments/copyO_oblique/`、`copyP_centered_smpc/` 和 `copyQ_control_aware/` 中逐步诊断和解决。
