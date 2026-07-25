# copyR：严格 Mo--Qin Algorithm-1 斜投影 PredVARX-SMPC

## 目的

`copyR_moqin_oblique` 是原始 PredVAR 理论的基线对齐实验。它不改动 `main/`、`copyO`、`copyP` 或 `copyQ`。

它严格采用 Mo & Qin (2025), *Automatica* 180, 112476 的归一化 IVR 与反归一化 realization：

\[
y_k^*=D^{-1/2}U^\top(y_k-\bar y),\qquad P^\top_*P^*=I,
\]

\[
P=UD^{1/2}P^*,\qquad R=UD^{-1/2}P^*.
\]

补空间使用原文 Eq. (34) 的同一反归一化结构：

\[
\bar P=UD^{1/2}\bar P^*,\qquad \bar R=UD^{-1/2}\bar P^*.
\]

因此自动满足：

\[
R^\top P=I,\quad R^\top\bar P=0,\quad
\bar R^\top P=0,\quad\bar R^\top\bar P=I.
\]

## 与 copyO 的区别

- **copyO** 在反归一化后额外做 SVD 对齐并以通用 dual completion 补齐基底。
- **copyR** 不做 QR、也不做后验 SVD 对齐；直接保留原文 Algorithm 1 的 $P=UD^{1/2}P^*$、$R=UD^{-1/2}P^*$ 与 Eq. (34) 补空间构造。
- 因而 copyR 是“原文 PredVAR realization 是否足以支持已修正 SMPC”的干净对照，而不是新的控制感知方法。

## PredVARX 与 SMPC 部分

原文 PredVAR 是 VAR，不包含控制输入。此实验明确把内层模型扩展为中心化 VARX：

\[
z_k=R^\top(y_k-\bar y),\qquad
z_{k+1}=Az_k+B(u_k-\bar u)+\varepsilon_k.
\]

SMPC 复用 copyP/copyQ 已验证的绝对输出预测和 Boole 联合机会约束。这个 VARX/SMPC 部分是本项目扩展，不归因于 Mo & Qin 原文。

## 解释边界

copyR **没有**强制控制覆盖：

\[
PR^\top e_i=e_i.
\]

所以若它的控制效果不如 copyQ，不能把问题归咎于“斜投影错误”；它展示的是原始 PredVAR 的预测性子空间与控制目标覆盖之间的实际差异。

## 文件

- `predvarx_identify_moqin.m`：Algorithm-1 对齐辨识器。
- `copyR_moqin_oblique.m`：中心化 VARX + Boole-SMPC 闭环实验。
- `tests/test_predvarx_identify_moqin.m`：反归一化、完整对偶基及 PSD 协方差回归测试。
- `results/`：可覆盖的 MAT、PNG 与运行日志，不提交到 Git。
