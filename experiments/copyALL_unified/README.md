# copyALL_unified：H/K/O/P/Q/R/S/T/U/V/X/Y/Z 全版本统一参数

## 目的

这是独立的新全版本公平套件，不覆盖已经跑通的 `copyOPQSTR_unified` 或各单独 copy。所有 13 个标签共享同一 plant、离线数据、参考、噪声样本、SMPC、指标和随机种子，只保留辨识几何与在线更新策略差异。

统一参数：

```text
seed=20260710
n=6, m=3, p=30, ell=5, tracked=[1 2]
T_off=1500, T_cl=1200, N=18
u in [-3,3], y_max=2.00, alpha_joint=0.10
smooth sigma_w=[0.020,0.090], sigma_e=[0.025,0.100]
```

## 版本定义

| 版本 | 全版本统一实现 |
|---|---|
| H | 最近 300 步闭环数据滑动窗口结构重辨识。 |
| K | 累积结构重辨识 + 50 步创新协方差 + EMA 偏差修正。 |
| O | 自由斜投影 IVR；不强制 tracked coverage。 |
| P | 全局 PCA/SVD 正交子空间。 |
| Q | tracked-axis 控制感知迭代 IVR，$R=P$。 |
| R | rank-full、$s=1$ 白化与直接反白化消融。 |
| S | 每 30 步用离线数据 + 最近 100 步重辨识。 |
| T | 固定 Q 结构，仅在线更新协方差。 |
| U | tracked-axis + 一次 SVD 补空间，不做迭代 IVR。 |
| V | tracked-axis + 补空间迭代 IVR；在本套件中与 Q 同定义。 |
| X | V 的 loading 子空间 + $\alpha=0.02$ 轻度斜双基 extractor。 |
| Y | 无白化、原始中心化坐标直接 eigenspace replacement；审计后与当前 X 数值同构。 |
| Z | 严格白化 rank-full、$s=1$ 消融；在本套件中与 R 同定义。 |

## 重复标签不是错误

- Q/V 数值相同：两者使用同一个 control-aware iterative-IVR identifier，只属于不同历史演化支线。
- X/Y 数值相同：当前 X 的源码审计显示其 IVR 本来就在 raw centered complement 中直接更新，并未实现论文白化；Y 只是把“无白化直接更新”显式命名。
- R/Z 数值相同：两者都调用同一 rank-full、$s=1$ strict-whitened ablation。Z 是后续 copyX plant 上的版本标签。

保留这些标签是为了完整呈现版本谱系，而不是制造虚假的算法差异。

## 运行

```matlab
run('experiments/copyALL_unified/tests/test_copyALL_unified.m')
run('experiments/copyALL_unified/copyALL_unified.m')
```

输出：

```text
results/copyALL_unified_metrics.csv
results/copyALL_unified_fig.png
results/individual/copyH_unified.png ... copyZ_unified.png
```

除一张总对比图外，入口脚本还会为 H/K/O/P/Q/R/S/T/U/V/X/Y/Z 各生成一张独立运行图。W 是比较器而非单一闭环算法，其独立图继续使用 `copyW_fair_identifier_compare` 的四辨识器公平比较图。

本轮真实 MATLAB 运行已通过结构测试并输出 13 行方法指标。PNG 大小为 467301 bytes；同步到 Obsidian 的 `15-copyALL-unified.png` 后，两份文件 SHA-256 均为 `f13f4c4219d934facb59657ea681065f3c8fd673eb7fa6bd6383482af5cefb09`。

## 本轮统一结果

| 版本 | MAE y1/y2 | prediction RMSE | coverage | QP成功率 | fallback | 重辨识 |
|---|---:|---:|---:|---:|---:|---:|
| H | 0.3577 / 0.1946 | 0.2227 | 0 | 82.25% | 213 | 4 |
| K | 0.0905 / 0.0964 | 0.0987 | 0 | 100% | 0 | 4 |
| O | 0.9685 / 1.1231 | 0.2842 | 1.9038 | 13.00% | 1044 | 0 |
| P | 0.0770 / 0.0788 | 0.0948 | 1.2995 | 100% | 0 | 0 |
| Q | 0.0917 / 0.0947 | 0.0995 | 0 | 100% | 0 | 0 |
| R | 0.4451 / 0.7320 | 0.1473 | 1.6327 | 43.92% | 673 | 0 |
| S | 0.0807 / 0.0801 | 0.0990 | 0 | 100% | 0 | 40 |
| T | 0.0917 / 0.0947 | 0.0995 | 0 | 100% | 0 | 0 |
| U | 0.0912 / 0.0949 | 0.1043 | 0 | 100% | 0 | 0 |
| V | 0.0917 / 0.0947 | 0.0995 | 0 | 100% | 0 | 0 |
| X | 0.0926 / 0.0947 | 0.0996 | 1.03e-14 | 100% | 0 | 0 |
| Y | 0.0926 / 0.0947 | 0.0996 | 1.03e-14 | 100% | 0 | 0 |
| Z | 0.4451 / 0.7320 | 0.1473 | 1.6327 | 43.92% | 673 | 0 |

## 源码边界

- 共享 SMPC 每拍由当前 `model.A/B/P/R` 现场重建预测矩阵。
- H/K/S 重辨识后新结构立即进入名义预测。
- R/Z 不是完整 Mo--Qin Algorithm 1，只是明确标注的 rank-full、$s=1$ 白化消融。
- W 本身是“比较器/实验容器”，不是第五种独立 identifier。因此 W 在全版本文档中作为公平比较层保留；其四个内部方法已由 P、O、V 以及 main-QR 基线覆盖。全版本主曲线不再把 W 当成一条伪算法曲线。