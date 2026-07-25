# 代码归档（非现行）

本目录存放**非现行** MATLAB 快照。新工作请只改仓库根下的 `experiments/` 与 `main/`。

| 子目录 | 内容 | 说明 |
|:---|:---|:---|
| `early_snapshot_2026-07-09/` | 早期独立快照：main + copyO/P/Q | 已被现行 `experiments/` 覆盖/超越，仅作历史对照 |
| `organized_legacy_2026-07-10/` | 旧实验按 01–10 主题归档 | 字母体系与现行不同；见根目录 `LEGACY_COPY_OPQSTR_ARCHAEOLOGY.md` |
| `routeC_delta_inc/` | v10 Δ-tracking 加速线（C/E/F/G/H/I/P1–P4） | 现行 `experiments/` 中无同名脚本 |

## 使用规则

1. **不要**把 archive 里的 Q/S/Z 与现行 `copyQ_control_aware` / `copyZ_strict_whitened_copyX_plant` 当成同一算法。
2. 比较指标前先确认 plant、噪声、时域、控制器是否一致。
3. 需要可复现的现行结果时，优先用 `experiments/copy*/results/`。

## 与 GitHub 发布的关系

该 archive 已随仓库一并公开发布（`.mat`/`.png` 走 Git LFS），目的是保留完整 MATLAB 版本谱系，而不是推荐继续在 archive 上开发。
