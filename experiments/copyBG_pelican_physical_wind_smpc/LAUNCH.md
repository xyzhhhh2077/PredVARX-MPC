# copyBG launch — majority near task hard (v4.0)

## Narrative (v4.0 default)

目标是让**实际闭环轨迹至少 50% 的时间贴近任务硬约束**，而不是只在四个轴尖短暂经过。

统一统计口径：

- task hard：`±3.912411728`
- near-hard 带：仍在 hard 内，且 `max_i |y_i| >= hard - 0.40 = 3.512411728`
- 越界样本不计入 near-hard，并单独要求 hard violation = 0

默认参考为 `z-face-large-arc-cruise`：在 `z=+3.57` 与 `z=-3.57` 两个任务硬约束面上持续飞水平大圆弧。

1. **Approach top** (15%, 1.5 圈)：`(0.90,0.90,1.20)` 平滑进入上表面
2. **Top cruise** (30%, 3 圈)：`z=+3.57`，水平大圆弧持续飞行
3. **Face transfer** (10%, 1 圈)：从 `+3.57` 平滑转移到 `-3.57`
4. **Bottom cruise** (30%, 3 圈)：`z=-3.57`，水平大圆弧持续飞行
5. **Return** (15%, 1.5 圈)：平滑返回初始高度/半径

水平半径为 `(1.20,1.20)`；仍用 `GROWING_ARC_CENTER_OFFSET=12` 的 inward large-arc diamond，因此不碰 `(±1.2,±1.2)` 多轴尖角。

### 为什么默认不让 x/y 面长期贴边

DC 诊断必须最小化峰值电机指令 `||u||_inf`，不能把最小二范数伪逆解当成该指标。修正后，在 near-hard 内缘 `3.5124`，x/y/z 单轴稳态输入分别为 `6.152 / 5.151 / 0.785`；在真正 hard `3.9124` 处分别为 `6.853 / 5.738 / 0.875`。因此 x 面确实超过 `±6`，但 y 面在输入意义下可以通过减慢运动、降低动态余量需求来接近 hard。

当前 SMPC 仍不能让预测均值长期停在 y hard：第 18 个预测步的 x/y/z 收紧量分别为 `0.977 / 0.764 / 0.329`，对应均值边界 `2.935 / 3.148 / 3.583`。该协方差收紧不随参考速度下降而消失，所以 y 单面即使输入可达，也不能仅靠减速达到“长期 near-hard”的当前 SMPC 合同。z 面的 `3.57 < 3.583`，且完整 `z=±3.57` + 水平大圆弧最坏稳态输入为 `u_inf=2.280 < 6`，因此默认选择 z 面。

### Physical wind（沿用 v3.8/v3.9）

- 长期 `sigma = 0.04 m/s`
- 慢变强度包络 depth=0.75 / ~40 s；白噪声驱动（无 AR 漂移）
- chance QP 仍只使用综合创新协方差 `Sigma_eps`；不访问真实 `Sigma_n`

## Run

```bash
python experiments/copyBG_pelican_physical_wind_smpc/run_physical_wind_smpc.py \
  --steps 48000 \
  --output-dir experiments/copyBG_pelican_physical_wind_smpc/results_growing_v40
```

## Official v4.0 run (`results_growing_v40/`, 48k full)

| 指标 | SMPC | deterministic MPC |
|---|---:|---:|
| completed | **48000** | **48000** |
| actual near-hard fraction | **52.3167%** (25112 步) | 52.3125% (25110 步) |
| hard violation | **0** | **0** |
| max `max_i |y_i|` | **3.82744** | 3.82744 |
| min hard margin | **0.08498** | 0.08498 |
| active chance QP | **18** | 0 |
| QP failure / fallback | **0 / 0** | **0 / 0** |
| max QP residual | `8.63e-11` | `0` |
| input saturation steps | 17030 (35.48%) | 17030 (35.48%) |
| RMSE x/y/z | 0.2437 / 0.2575 / 0.0707 | 0.2437 / 0.2574 / 0.0707 |

参考 near-hard 占比为 **63.4229%**，实际闭环为 **52.3167%**；两者都在 hard 内，满足“至少 50% 时间贴近任务硬约束且零越界”。

阈值敏感性（同一 SMPC 实际轨迹）：离 hard 不超过 `0.10 / 0.20 / 0.30 / 0.40` 时，占比分别为 `0.029% / 1.992% / 11.688% / 52.317%`。因此 v4.0 **只在预先声明的 0.40 带宽口径下达到 50%**；不能写成“50% 时间都在 3.8 以上”。0.40 带内样本的 `max_i|y_i|` 中位数为 `3.5741`。

输入也没有超过 `±6`，但有 **35.48%** 的样本恰好达到饱和边；因此这里可称“零输入越界”，不能称“执行器远离约束”。这里的 `u` 是标准化电机指令

```text
u_std = (Motors_CMD - u_offset) / u_scale
```

`±6` 是从 copyBA/copyBE 沿用的实验设计边界，约对应训练均值的 `±6` 个标准差，**不是 Pelican 硬件铭牌极限**。按本次训练统计量映回原始无量纲 `Motors_CMD`，四通道约为 `[52.9,197.7] / [47.8,199.3] / [45.8,188.7] / [42.0,181.4]`。数据说明给出的命令域是 `[0,218]`，但本副本实际 54 段飞行的逐通道观测范围更窄，约为 `M1 [45.8,198.9] / M2 [45.0,183.4] / M3 [46.2,183.6] / M4 [44.8,184.4]`；因此扩大控制边界必须标成模型外推试验，不能直接称为数据内或硬件已验证。

> 口径边界：v4.0 的主要目标是**持续贴 hard**。chance QP 确实激活 18 步，但 SMPC 与 MPC 轨迹差异很小；它不是比 v3.9/boundary-tour 更强的 SMPC-vs-MPC 分离演示。若要突出 chance tightening 的差异，仍用 `--reference boundary-tour` 作为专门对照。

## Artifacts

- `copyBG_pelican_physical_wind_smpc.json`：完整指标与 near-hard contract
- `copyBG_pelican_physical_wind_smpc.npz` / `.mat`：同一 48k 实跑数据
- `copyBG_pelican_physical_wind_smpc_smpc.png`
- `copyBG_pelican_physical_wind_smpc_deterministic_mpc.png`
- `copyBG_pelican_physical_wind_smpc_noise_actual_vs_estimated.png`
- `copyBG_pelican_physical_wind_smpc_smpc_playback.gif`：200 帧 × 120 ms 固定机位

## Scope

- 未改 CETR/CRTE 论文理论
- 未改 `copyAU/`
- 未改 `copyBF_pelican_long_3d_boundary_disturbance/run_continuous_3d_boundary_circuit.py`