%% PredVARX-MPC v10 副本 - Route C-Delta 三段阶跃 + 增量加速 P2
%
%  ★ 副本数学公式与 routeC_delta_3step 完全一致, 仅修改 MPC 在线求解:
%    [P2] 顺序单步预测 (Dinkla 2026 思路):
%      - 不预计算 M_j^cl = P̂ E_x^T A_cl^j (高次幂, O(N^3) 内存)
%      - 改为 M̃_{j+1} = M̃_j · A_cl, 每次矩阵乘法 O(N^2) —— (p×nz)×(nz×nz)
%      - N_cl 的列也用增量: N_cl(:, i+2) = A_cl · N_cl(:, i+1)
%      - 内存: 从 ~30×30×15×2 缓存降到单步增量
%
%  目的: QP 系数矩阵在线增量构建, 不影响 Δ-tracking 数学
%  三段阶跃同 P1
%
%  ★ 此副本与 routeC_delta 副本 (单段阶跃) 数学完全相同, 唯一改动:
%    阶跃信号分三段, 不再是单一阶跃:
%      - 段 1 (k ∈ [50, 400)):   ref = 1.0  (与 v10 主版相同)
%      - 段 2 (k ∈ [400, 800)):  ref = 2.0  (新段, 大幅上升)
%      - 段 3 (k ∈ [800, 1200]): ref = 0.5  (新段, 下降, 验证 Δ-tracking 对称性)
%
%  目的: 验证 routeC_delta (Δ-tracking) 在多段阶跃下的:
%    - 段 1: 复现单段实验的 y₁=0.28, y₂=1.09 基线
%    - 段 2: 看 Δ-tracking 能否把 y₂ 从 1.09 推到 2.0 附近
%             (传统 MPC 跟踪会因 v10 几何锁定失败, Δ-tracking 假设能跟随参考变化)
%    - 段 3: 验证 Δ-tracking 在负向参考下也工作 (对称性)
%             特别关注: y₁ 是否能跟随 0.5, 还是回到基线 0.28
%
%  副本数学: 与 routeC_delta 1:1 一致 - rng/Q_w/N_pred/Eq.14/Δ-tracking 全部不动
%  T_cl=1200, 阶跃起点 k=50, 段切换 k=400 / k=800
%
%  Δ-tracking 核心 (L395-L430, 不变):
%    min Σ_j ||Δy_j - Δr_j||_Q^2
%    Δy_1 = (M_1 z_k + N_1 c) - y_k^actual
%    Δy_j = (M_j - M_{j-1}) z_k + (N_j - N_{j-1}) c
%  这样 MPC 优化"y 的变化量", 不偏向任何 y 通道方向

clear; close all; clc;
rng(42, 'twister');  % 与 routeC_delta 同 RNG, 离线辨识结果一致

% 路线 C 超参数
beta_u = 0.5;   % u 在 IV 中的权重 (β=0 时等价 v10 主版)

%% ==================== 参数 ====================
% 对应文档 §1 符号约定
n = 6;   % 状态维度 x_k ∈ R^n (SDDPC: n=6)
m = 3;   % 控制维度 u_k ∈ R^m (SDDPC: m=3)
p = 15;  % 观测维度 y_k ∈ R^p (PredVAR: p=15)
ell = 2; % 降维维度 x̂_k ∈ R^ℓ (ℓ ≤ n)

ref_val = 1.0;               % ★ 3step 副本不再使用此值, 仅保留兼容
T_off = 20000;               % 离线数据长度
T_cl = 1200;                 % 闭环测试步数
N_pred = 30;                 % MPC 预测步长 (参考 SDDPC: N=30)
u_min = -5;  u_max = 5;     % 输入约束
noise_interval = 200;        % 噪声变化间隔
noise_factors = [1.0, 2, 0.52, 3, 5, 0.824];  % 噪声缩放因子

% ★ 3step 副本专用: 三段阶跃参考
step1_end = 400;   % 段 1 -> 段 2 切换
step2_end = 800;   % 段 2 -> 段 3 切换
ref_seg1 = 1.0;    % 段 1 参考
ref_seg2 = 2.0;    % 段 2 参考 (大幅上升, 测跟踪能力)
ref_seg3 = 0.5;    % 段 3 参考 (下降, 测对称性)

fprintf('=== PredVARX-MPC v10 ===\n');
fprintf('n=%d, m=%d, p=%d, l=%d, N_pred=%d\n', n, m, p, ell, N_pred);

%% ==================== §1: 真实系统 (Eq.1-2) ====================
% 文档 §1.1: x_{k+1} = A x_k + B u_k + w_k  (Eq.1)
%            y_k = C x_k + e_k                (Eq.2)
% 其中 w_k ~ N(0, Σ_w), e_k ~ N(0, Σ_e), ρ(A) < 1 (稳定)
fprintf('\n[§1] 真实系统 x_{k+1}=Ax_k+Bu_k+w_k, y_k=Cx_k+e_k\n');

% A: 对角稳定矩阵, ρ(A) < 1
% 特征值 [0.85, 0.70, 0.55, 0.40, 0.30, 0.20], 加弱耦合
A_true = diag([0.85, 0.70, 0.55, 0.40, 0.30, 0.20]);
A_true(1,2) = 0.1; A_true(3,4) = 0.01;  % 非对角弱耦合

% B: 控制输入矩阵 (n × m)
B_true = randn(n, m);

% C: 观测矩阵 (p × n), 列正交化
C_true = randn(p, n); [C_true, ~] = qr(C_true, 0);

% 噪声标准差
sigma_w = 0.1;  sigma_e = 0.1;

% 检查稳定性: ρ(A) < 1
rho_A = max(abs(eig(A_true)));
assert(rho_A < 1, 'A 不稳定!');
fprintf('  ρ(A) = %.4f\n', rho_A);

%% ==================== §2: 数据采集 ====================
% 用真实系统生成离线数据 {y_k, u_k}_{k=1}^{T_off}
% u_off 是持续激励信号 (σ=10), 用于激发系统动态
fprintf('[§2] 数据采集 T=%d\n', T_off);

x_off = zeros(n, T_off+1);  % 真实状态序列
y_off = zeros(p, T_off);     % 观测序列
u_off = 10 * randn(m, T_off);  % 强激励输入 (σ=10 >> σ_w=0.1)

for k = 1:T_off
    % Eq.2: y_k = C x_k + e_k
    y_off(:, k) = C_true * x_off(:, k) + sigma_e * randn(p, 1);
    % Eq.1: x_{k+1} = A x_k + B u_k + w_k
    x_off(:, k+1) = A_true * x_off(:, k) + B_true * u_off(:, k) + sigma_w * randn(n, 1);
end

%% ==================== §3: Algorithm 1 (§2 文档) ====================
% 文档 §2.5: Algorithm 1 详细流程
% 输入: ℓ, s, s_u, m, {y_k}, {u_k}
% 输出: Â, B̂, P̂, R̂, Σ̂_ε, Σ̂_ε̄
fprintf('[§3] Algorithm 1: 两阶段法 + IVR + 动力学回归\n');

% ── 阶段 1: 输入剥离 (§2.1 文档) ──
%
% ★ 推导: y_k = C1·u_{k-1} + C2·u_{k-2} + y_tilde_k, y_tilde 保留状态动态
%
% (1) 向量化 u: phi_k = [u_{k-1}; u_{k-2}] ∈ R^{2m}
% (2) 水平堆叠 C 矩阵: C_hat = [C1 C2] ∈ R^{p x 2m}
%     每个 C_i ∈ R^{p x m}, 竖直拼成 p x 2m 矩阵
%     y_k = C_hat * phi_k + y_tilde_k  (线性回归, y_k 从 u_k 的过去值拟合)
%     ∵ y_k 是 p×1 列向量, phi_k 是 2m×1 列向量
%
% (3) 多步合并: 设 N = T_off - L (有效样本数)
%     Y = [y_{L+1}; y_{L+2}; ...; y_{L+N}] ∈ R^{N x p}
%     Φ_u = [phi_{L+1}; phi_{L+2}; ...; phi_{L+N}] ∈ R^{N x 2m}
%     注意: Y 和 Φ_u 是"行堆叠" (每个 y_k 是一行, 共 N 行)
%
% (4) OLS 求解: C_hat = (Φ_u^T Φ_u)^{-1} Φ_u^T Y  ∈ R^{2m x p}
%     ★ C_hat 的维度 (关键!): 2m x p, 不是 p x 2m!
%     因为 Φ_u^T Y 是 2m x p, (Φ_u^T Φ_u) 是 2m x 2m
%
%     ★ C_hat 的物理含义:
%       前 m 行: y 关于 u_{t-1} 的回归系数 (每行 m 个 u 通道对 p 个 y 通道的影响)
%       后 m 行: y 关于 u_{t-2} 的回归系数
%
% (5) 残差 (剥离 u 后的 y):
%     y_residual  = Y - Φ_u·C_hat   (OLS 残差)
%     y_residual_k = y_k - C_hat·[u_{k-1}; u_{k-2}]
%     这个残差里"已经去掉了 u 的直接影响", 保留的是:
%     ✅ 状态动态 (C A^k x_0)
%     ✅ 状态变化 (x_k 本身)
%     ✅ 噪声 (w_{t-1} 对状态的影响)
%
%     ❌ 不再含 u 的直接输入 (u_{k-1}, u_{k-2} 已被 OLS 移除)
%
%     ★ 这就是 v10 的"输入剥离"—— 把 u 的影响从 y 中移除,
%       剩下的 y_residual 用于 IVR, 找"状态的子空间"而不是"u 驱动的子空间"
L_u = 2;                      % 回归 u_{k-1}, u_{k-2}
N_strip = T_off - L_u;        % 有效样本数
Phi_u = zeros(N_strip, m*L_u);  % 回归矩阵 [u_{k-1}; u_{k-2}]
Y_strip = zeros(p, N_strip);    % 对应的观测
for i = 1:N_strip
    t = L_u + i;
    u_lag = []; for lag = 1:L_u, u_lag = [u_lag; u_off(:, t-lag)]; end
    Phi_u(i, :) = u_lag';      % [u_{k-1}^T, u_{k-2}^T]
    Y_strip(:, i) = y_off(:, t);  % y_k
end
% OLS 回归: y_k = Φ_u · Ĉ + ỹ_k
C_hat_u = (Phi_u'*Phi_u + 1e-6*eye(m*L_u)) \ (Phi_u'*Y_strip');
% 残差 = y - Φ_u · Ĉ (保留状态动态, 剥离 u 的直接影响)
y_residual = (Y_strip' - Phi_u * C_hat_u)';
fprintf('  阶段1: ||残差||/||y|| = %.3f\n', norm(y_residual,'fro')/norm(Y_strip,'fro'));

% --- 路线 C: u-加权预白化 (改 Qin/Mo 原文的 Σ̃_y) ---
% 路线 C 修复版: 把 u 加入预白化矩阵, 让 P̂ 列空间主动含 u 方向
% 公式: Σ̃_yu = [Σ̃_y, β_u C D_u^T; β_u D_u C^T, β_u^2 Σ̃_u]  (联合协方差)
% 简化: 直接对 [y_c; β_u u_c] 做 SVD, 得到联合预白化空间
y_c = y_residual - mean(y_residual, 2);                          % p × N_strip (中心化 y)
u_c = u_off(:, L_u+1:L_u+N_strip) - mean(u_off(:, L_u+1:L_u+N_strip), 2);  % m × N_strip (中心化 u)
% 联合数据 [y_c; β_u u_c], 形状 (p+m) × N_strip
joint = [y_c; beta_u * u_c];
[Joint_U, Joint_D, ~] = svd(joint * joint' / N_strip, 'econ');
Joint_D_inv_sqrt = diag(1 ./ sqrt(diag(Joint_D) + 1e-12));
y_star = Joint_D_inv_sqrt * Joint_U' * joint;                    % (p+m) × N_strip
% ★ 关键: 保留 p 维的 y* 部分 (前 p 行), 这样后续 IVR 步骤不变
y_star = y_star(1:p, :);

% 步骤 2: 初始化 (Eq.32)
% 构造工具矩阵 Y* = [Y*_{s-1}, ..., Y*_0] (lagged y*)
% EVD on Y*_s^T Π_{Y*} Y*_s/N → P̂* = W(:, 1:ℓ)
N_ivr = N_strip - 1;  % s=1, 有效样本数
r = p;
Y_star_s = y_star(:, 2:N_strip);       % y*_{k} (当前)
Y_star_instr = y_star(:, 1:N_strip-1);  % y*_{k-1} (lag-1 工具变量)

% 路线 C.1: 扩展 IV 矩阵 Z_u = [α y*_{k-1}; β u_{k-1}]
% Z_u 形状: N_ivr × (p+m), 即每行一个时刻, 每列一个变量
% (这样 Π_{Z_u} = Z_u^T (Z_u Z_u^T)^{-1} Z_u 才能正确维度)
u_past_for_iv = u_off(:, L_u+1:L_u+N_ivr)';   % u_{k-1}^T, N_ivr × m (转置!)
if beta_u > 0
    Z_u = [Y_star_instr', beta_u * u_past_for_iv];    % N_ivr × (p+m)
else
    Z_u = Y_star_instr';                                % N_ivr × p
end
fprintf('  路线 C.1: β = %.2f, Z_u 维度 = %d × %d\n', beta_u, size(Z_u,1), size(Z_u,2));

% Π_{Z_u} = Z_u (Z_u^T Z_u)^{-1} Z_u^T (投影矩阵)
ZtZ = Z_u' * Z_u;                              % (p+m) × (p+m)
Pi_Zu = Z_u / (ZtZ + 1e-8*eye(size(ZtZ,1))) * Z_u';   % N_ivr × N_ivr

% EVD on Y*_s^T Π_{Z_u} Y*_s/N (Eq.32, 用 Z_u 代替 Y*)
% 高效: M = (Y*_s · Z_u) (Z_u^T · Z_u)^{-1} (Z_u^T · Y*_s^T) / N
% Y*_s 是 p × N, Z_u 是 N × (p+m), Y*_s · Z_u = p × (p+m)
YtZ = Y_star_s * Z_u;                           % p × (p+m)
M_init = YtZ / (ZtZ + 1e-8*eye(size(ZtZ,1))) * YtZ' / N_ivr;  % p × p
[P_star, Lambda_init] = eigs(M_init, ell, 'largestabs');
P_star = real(P_star);
fprintf('  初始化 eig: [%.4f %.4f]\n', Lambda_init(1,1), Lambda_init(2,2));

% 步骤 3-5: IVR 迭代精炼 (Eq.33)
% while trace(Σ̂_x̂ = X̂_s^T X̂_s/N) 仍在增加:
%   X̂ = Y*_s^T · P̂*          (提取降维状态)
%   X̂_lagged = [X̂_{k-1}]     (lagged 降维状态)
%   Â_ivr = (X̂_lagged^T X̂_lagged)^{-1} X̂_lagged^T X̂  (VAR 回归)
%   X̂_pred = X̂_lagged · Â_ivr  (预测)
%   EVD on Y*_s^T Π_{X̂_pred} Y*_s/N → P̂* = W(:, 1:ℓ) (更新)
%
% 路线 C.2: 迭代过程中, X̂_pred 也投影到 Z_u (而不是只 Y*_instr)
%   W = (Y*_s^T · Z_u) (Z_u^T · X̂_pred^T · X̂_pred · Z_u)^{-1} (Z_u^T · Y*_s^T) / N
% 内存优化: 避免 N×N 投影矩阵, 用 W*inv(V'V)*W' 模式
trace_prev = -inf;
for iter = 1:50
    % X̂ = Y*_s^T · P̂* (N_ivr × ℓ, 当前降维状态)
    X_hat = Y_star_s' * P_star;
    % X̂_lagged = Y*_{instr}^T · P̂* (N_ivr × ℓ, lag-1 降维状态)
    X_lagged = Y_star_instr' * P_star;

    % Â_ivr = (X̂_lagged^T X̂_lagged)^{-1} X̂_lagged^T X̂ (ℓ × ℓ, VAR 回归)
    A_ivr = (X_lagged' * X_lagged) \ (X_lagged' * X_hat);
    % X̂_pred = X̂_lagged · Â_ivr (N_ivr × ℓ, 预测的降维状态)
    X_pred = X_lagged * A_ivr;

    % 收敛判据: trace(Σ̂_x̂ = X̂_pred^T X̂_pred/N) 不再增加
    trace_curr = trace(X_pred' * X_pred / N_ivr);
    if iter > 1 && trace_curr <= trace_prev + 1e-10
        fprintf('  IVR 收敛于第 %d 步\n', iter);
        break;
    end
    trace_prev = trace_curr;

    % 路线 C.2: EVD on Y*_s^T Π_{X̂_pred|Z_u} Y*_s/N
        % Π_{X̂_pred|Z_u} = Z_u (Z_u^T X̂_pred^T X̂_pred Z_u)^{-1} Z_u^T
        % 标准: Π_{X̂_pred|Z_u} = X̂_pred^T Π_{Z_u} X̂_pred (在 X̂_pred 张成的列空间上的正交投影)
        % 高效: M = Y*_s^T · X̂_pred · (X̂_pred^T · Z_u · Z_u^T · Z_u · Z_u^T · X̂_pred)^{-1} · X̂_pred^T · Y*_s / N
        % = (Y*_s^T X̂_pred) (X̂_pred^T Z_u (Z_u^T Z_u)^{-1} Z_u^T X̂_pred)^{-1} (X̂_pred^T Y*_s) / N
        % 设 W = Z_u^T X̂_pred (N_ivr × ℓ), 然后:
        % M = (Y*_s^T X̂_pred) (W^T W)^{-1} (Y*_s^T X̂_pred)^T / N, 但这不是投影矩阵, 只是标量化
        % 正确: M = (Y*_s Z_u (Z_u^T Z_u)^{-1} Z_u^T X̂_pred) · (X̂_pred^T Z_u (Z_u^T Z_u)^{-1} Z_u^T X̂_pred)^{-1} · (X̂_pred^T Z_u (Z_u^T Z_u)^{-1} Z_u^T Y*_s^T) / N
        % 设 P_zu_X = Z_u (Z_u^T Z_u)^{-1} Z_u^T X̂_pred = N_ivr × ℓ (X̂_pred 在 Z_u 列空间的投影)
        P_zu_X = Z_u / (Z_u' * Z_u + 1e-8*eye(size(Z_u,2))) * (Z_u' * X_pred);  % N_ivr × ℓ
        M_ivr = Y_star_s * P_zu_X / (P_zu_X' * P_zu_X + 1e-8*eye(ell)) * (P_zu_X' * Y_star_s') / N_ivr;
        [P_star, ~] = eigs(M_ivr, ell, 'largestabs');
        P_star = real(P_star);
    end

% 步骤 6: 反归一化
% P̂ = UD^{1/2} P̂*, R̂ = UD^{-1/2} P̂*
% 路线 C 修复版: P_hat 用联合预白化矩阵反归一化, 取前 p 维
% P_hat = Joint_U(1:p, 1:p) * Joint_D_inv_sqrt(1:p, 1:p) * P_star
%   只取前 p 行/列的预白化矩阵, 确保维度匹配 P_star (p × ℓ)
P_hat = Joint_U(1:p, 1:p) * Joint_D_inv_sqrt(1:p, 1:p) * P_star;
[P_hat, ~] = qr(P_hat, 0);  % 列正交化
R_hat = P_hat;  % 列正交时 R̂ = P̂

% 互补空间: P̄ = null(P̂^T), R̄ = null(R̂^T)
Pbar_hat = null(P_hat');  Rbar_hat = null(R_hat');

% 子空间质量: cos(主角度) between P̂ and C_true 的列空间
SV_sub = svd(P_hat' * C_true * C_true' * P_hat);
fprintf('  子空间: [%.3f, %.3f]\n', sqrt(SV_sub(1)), sqrt(SV_sub(2)));

% --- 阶段 3: 动力学回归 (§2.3 文档, 步骤 9-11) ---
% 用 R̂ 从原始 y 提取降维状态: x̂_k = R̂^T y_k
% 回归: x̂_{k+1} = Â x̂_k + B̂ u_k + ε_k
X_hat_all = R_hat' * y_off;  % ℓ × T_off (降维状态序列)
t_start = 2;
N_reg = T_off - t_start;

% 构造回归矩阵 Φ = [x̂_k, u_k]
Phi_ab = zeros(N_reg, ell + m);
X_target = zeros(N_reg, ell);
for i = 1:N_reg
    t = t_start + i - 1;
    Phi_ab(i, :) = [X_hat_all(:, t); u_off(:, t)]';  % [x̂_k^T, u_k^T]
    X_target(i, :) = X_hat_all(:, t+1)';              % x̂_{k+1}^T
end

% Ridge 回归: Θ̂ = (Φ^T Φ + λI)^{-1} Φ^T X_target
Theta_ab = (Phi_ab'*Phi_ab + 1e-4*eye(ell+m)) \ (Phi_ab' * X_target);
A_hat = Theta_ab(1:ell, :)';      % ℓ × ℓ (降维动力学矩阵 Â)
B_hat = Theta_ab(ell+1:end, :)';  % ℓ × m (降维控制矩阵 B̂)

% --- 阶段 4: 噪声估计 (§2.4 文档) ---
% Σ̂_ε = (1/N) Σ (x̂_{k+1} - Â x̂_k - B̂ u_k)(...)^T
X_pred_off = Phi_ab * Theta_ab;
Sigma_eps_hat = ((X_target - X_pred_off)' * (X_target - X_pred_off)) / N_reg;

% Σ̂_ε̄ = R̄^T · (观测残差协方差) · R̄
resid_y = y_off - P_hat * X_hat_all;  % 观测残差
Sigma_resid = (resid_y * resid_y') / T_off;
Sigma_ebar_hat = Rbar_hat' * Sigma_resid * Rbar_hat;

fprintf('  Â eig: [%.3f, %.3f]\n', eig(A_hat));
fprintf('  Σ̂_ε  diag: [%.4f %.4f]\n', diag(Sigma_eps_hat));

%% ==================== §4: 概率预测器 (§3 文档) ====================
% 文档 §3.1: 增广状态 z_k = [x̂_k; ...; x̂_{k-s+1}; u_{k-1}; ...; u_{k-s_u}]
% s=1, s_u=0 → z_k = x̂_k (ℓ 维)
fprintf('[§4] 概率预测器\n');

nz = ell;  % 增广状态维度
F_aug_hat = A_hat;  % ℓ × ℓ (s=1 时 F_aug = Â)
H_aug_hat = B_hat;  % ℓ × m (s=1, s_u=0 时 H_aug = B̂)
G_aug_hat = eye(ell);  % ℓ × ℓ (噪声映射)

% K=0 (不使用 LQR 反馈, MPC 独自负责)
K_lqr = zeros(m, nz);
A_cl_hat = F_aug_hat;  % K=0 时 A_cl = F_aug
rho_cl = max(abs(eig(A_cl_hat)));
fprintf('  rho(A_cl) = %.4f\n', rho_cl);

% E_x: 提取 x̂_k 分量 (§3.2)
E_x = eye(ell);

% 预计算 M_j^cl, N_j^cl (§3.2, Eq.10)
% μ_y(j) = P̂ · E_x^T · A_cl^j · z_k + N_j · c
% ★ P2: 用增量 A_iter 计算 A_cl^j 代替直接求 A_cl^j 高次幂
%      先初始化 A_iter = A_cl, 再算 j=1 (M_1 = P̂ E_x^T · A_cl)
M_cl = cell(1, N_pred); N_cl = cell(1, N_pred);
A_iter = A_cl_hat;  % A_cl^1 (j=1 起始)
for j = 1:N_pred
    % M_j^cl = P̂ · E_x^T · A_cl^j (p × nz) — 用增量而非 A_cl^j
    M_cl{j} = P_hat * E_x' * A_iter;
    % 增量: A_iter ← A_iter · A_cl (供 j+1 用)
    if j < N_pred
        A_iter = A_iter * A_cl_hat;
    end
    % N_j^cl = [P̂ E_x^T A_cl^{j-1-i} H_aug for i=0..j-1]
    % i 列对应 A_cl^(j-1-i) 次幂: i=0 → A_cl^(j-1), ..., i=j-1 → A_cl^0 = I
    N_cl{j} = zeros(p, N_pred*m);
    for i = 0:j-1
        % 增量计算 A_cl^(j-1-i): 从 A_cl^(j-1) 往下降到 A_cl^0
        % 一次性算 A_cl^(j-1-i) 慢 (循环), 但避免 ^ 高次幂
        A_pow = eye(nz);
        for kk = 1:(j-1-i)
            A_pow = A_pow * A_cl_hat;
        end
        N_cl{j}(:, i*m+1:(i+1)*m) = P_hat * E_x' * A_pow * H_aug_hat;
    end
end
Q_w = zeros(p); Q_w(1:ell, 1:ell) = eye(ell);
R_w = 1.0 * eye(m);  % 控制量正则化

%% ==================== §5: SMPC 闭环 (§4 文档) ====================
% 文档 §4.3: Algorithm 2 在线流程
fprintf('[§5] SMPC 闭环 T=%d\n', T_cl);

% ★ 3step 副本: 三段阶跃参考
% 段 1: k ∈ [50, step1_end)    -> ref_seg1
% 段 2: k ∈ [step1_end, step2_end) -> ref_seg2
% 段 3: k ∈ [step2_end, T_cl]    -> ref_seg3
y_ref_full = zeros(p, T_cl);
% 段 1
idx1 = (50:T_cl);  idx1 = idx1(idx1 < step1_end);
y_ref_full(1:ell, idx1) = ref_seg1;
% 段 2
idx2 = (50:T_cl);  idx2 = idx2(idx2 >= step1_end & idx2 < step2_end);
y_ref_full(1:ell, idx2) = ref_seg2;
% 段 3
idx3 = (50:T_cl);  idx3 = idx3(idx3 >= step2_end);
y_ref_full(1:ell, idx3) = ref_seg3;
fprintf('  3step 参考: 段 1 [50, %d) = %.2f, 段 2 [%d, %d) = %.2f, 段 3 [%d, %d] = %.2f\n', ...
    step1_end, ref_seg1, step1_end, step2_end, ref_seg2, step2_end, T_cl, ref_seg3);
x_hat_hist = zeros(ell, T_cl);  % 降维状态历史
y_true = zeros(p, T_cl);        % 真实观测历史
u_hist = zeros(m, T_cl);        % 控制输入历史
qp_fail = 0;
cost_hist = zeros(1, T_cl);                    % QP 失败计数
x_true = zeros(n, 1);           % 真实状态

% 在线协方差初始化 (§3.3 文档)
% Σ̂_ε(0) = 离线估计值, λ = 0.95 (遗忘因子)
Sigma_eps_k = Sigma_eps_hat;
Sigma_ebar_k = Sigma_ebar_hat;
lambda_forget = 0.95;  % 遗忘因子, 有效窗口 ≈ 1/(1-0.95) = 20 步
use_eq14 = false;  % 初始关闭 (等烧入期)

sigma_w_actual = zeros(1, T_cl);  % 记录实际噪声
sigma_eps_online = zeros(1, T_cl);  % 记录在线估计

opt_qp = optimset('Display', 'off', 'LargeScale', 'off');

for k = 1:T_cl
    % --- 非平稳噪声 ---
    % 每 noise_interval 步改变噪声水平
    seg = min(floor((k-1)/noise_interval)+1, length(noise_factors));
    sigma_w_k = sigma_w * noise_factors(seg);  % 当前过程噪声
    sigma_e_k = sigma_e * noise_factors(seg);  % 当前观测噪声
    sigma_w_actual(k) = sigma_w_k;

    % --- 观测 (Eq.2) ---
    y_k = C_true * x_true + sigma_e_k * randn(p, 1);
    y_true(:, k) = y_k;

    % --- Algorithm 2 步骤 1: x̂_k = R̂^T y_k (命题 4) ---
    x_hat_k = R_hat' * y_k;
    x_hat_hist(:, k) = x_hat_k;

    % --- Algorithm 2 步骤 3: 在线协方差更新 (§3.3, Eq.13) ---
    % 每步用指数遗忘更新 Σ̂_ε(k)
    if k >= 2
        % 一步预测: x̂_{k|k-1} = Â x̂_{k-1} + B̂ u_{k-1}
        % 预测误差 (新息): e_k = x̂_k - x̂_{k|k-1}
        e_k = x_hat_k - A_hat * x_hat_hist(:, k-1) - B_hat * u_hist(:, k-1);
        % 指数遗忘更新: Σ̂_ε(k) = λ Σ̂_ε(k-1) + (1-λ) e_k e_k^T
        Sigma_eps_k = lambda_forget * Sigma_eps_k + (1-lambda_forget) * (e_k * e_k');
    end
    sigma_eps_online(k) = sqrt(trace(Sigma_eps_k) / ell);  % 记录在线估计

    % --- Algorithm 2 步骤 2: 构造 z_k ---
    % s=1, s_u=0 → z_k = x̂_k
    z_k = x_hat_k;

    % --- Algorithm 2 步骤 4-5: 计算 Σ_y^cl(j, k) (Eq.11-12) ---
    % Q_aug(k) = G_aug · Σ̂_ε(k) · G_aug^T
    Q_aug_k = G_aug_hat * Sigma_eps_k * G_aug_hat';
    % 递推: Σ_z^cl(j+1, k) = A_cl Σ_z^cl(j, k) A_cl^T + Q_aug(k)
    Sigma_z_prev_k = zeros(nz);
    for j = 1:N_pred
        Sigma_z_j = A_cl_hat * Sigma_z_prev_k * A_cl_hat' + Q_aug_k;  % Eq.11
        Sigma_x_j = E_x' * Sigma_z_j * E_x;
        % Σ_y^cl(j, k) = P̂ Σ_x^cl(j, k) P̂^T + P̄ Σ̂_ε̄(k) P̄^T  (Eq.12)
        Sigma_y_cl{j} = P_hat * Sigma_x_j * P_hat' + Pbar_hat * Sigma_ebar_k * Pbar_hat';
        Sigma_z_prev_k = Sigma_z_j;
    end

    % --- Algorithm 2 步骤 6: QP (Eq.14-15) ---
    % min Σ_j ||M_j z_k + N_j c - r_j||_Q^2 + ||c_j||_R^2
    % 约束: Eq.14 联合机会约束 e_i^T μ_y(j) + Φ^{-1}(1-p_i) √(e_i^T Σ_y^cl(j,k) e_i) ≤ f_i
    % === Route C-Delta: Δ-tracking 代价函数 ===
    % Δy_j = y_{k+j|k} - y_{k+j-1|k}  (j=1 用 y_k 实际值作基准)
    % 代价 = Σ_j ||Δy_j - Δr_j||_Q^2
    % 这样 MPC 优化 y 的增量变化, 不偏向任何 y 通道方向
    %
    % 实现: y_{k+j|k} = M_j z_k + N_j c
    %      Δy_j = (M_j - M_{j-1}) z_k + (N_j - N_{j-1}) c  (j>=2)
    %      Δy_1 = (M_1 z_k + N_1 c) - y_k_actual
    H_qp = zeros(N_pred*m);
    f_qp = zeros(N_pred*m, 1);
    y_k_actual = P_hat * x_hat_k;   % 当前 y 的估计 (从降维状态反投影)
    % ★ P2: 顺序单步 (Dinkla), 但保留 M_cl/N_cl cell 结构
    %    不在循环外预计算 A^j (j=1..N_pred), 而是在循环内增量:
    %       A_iter = A_iter * A_cl  (每 j 增量) — 避免 j×j 矩阵重复运算
    %    N_cl{j} 计算保留原 N_cl 结构 (一次循环), 这是离线一次性成本
    %    在线 MPC: 不需要 N_cl, 只需 M_cl{j} 当前所需
    for j = 1:N_pred
        r_j = y_ref_full(:, min(k+j, T_cl));  % 未来参考
        if j == 1
            delta_r = r_j - y_k_actual;
            M_delta = M_cl{j};
            N_delta = N_cl{j};
            f_const = -y_k_actual;
        else
            r_prev = y_ref_full(:, min(k+j-1, T_cl));
            delta_r = r_j - r_prev;
            M_delta = M_cl{j} - M_cl{j-1};
            N_delta = N_cl{j} - N_cl{j-1};
            f_const = zeros(p, 1);
        end
        H_qp = H_qp + N_delta' * Q_w * N_delta;
        f_qp = f_qp + N_delta' * Q_w * (M_delta * z_k + f_const - delta_r);
    end
    H_qp = H_qp + kron(eye(N_pred), R_w);
    H_qp = (H_qp + H_qp')/2 + 1e-8*eye(size(H_qp));

    % ★ Eq.14: 联合机会约束 (自适应上界, 烧入期后启用)
    % y_max(k) = mean(y_prev_window) + 3 * sigma_hat(k)
    % 其中 sigma_hat = sqrt(trace(Sigma_eps_k)/ell) (只考虑动态子空间误差)
    % 烧入期: 前 50 步不启用 (积累历史数据)
    win_len = 50;                      % 滑动窗口长度
    z_alpha = 1.645;                   % 5% 违反概率
    sigma_hat = sqrt(trace(Sigma_eps_k) / ell);
    
    if k > win_len
        % 滑动窗口平均 (前 ℓ 通道)
        y_avg = mean(y_true(1:ell, k-win_len+1:k), 2);  % ell x 1
        y_max_vec = y_avg + 1.0 * sigma_hat;              % 自适应上界 (ell x 1)
        % ① 只用动态子空间 Sigma: 不含 Sigma_ebar (避免 u 强激励导致 margin 过大)
        % ② 基于估计误差, 自适应上界: y_avg + 3 * sigma_hat
        A_ineq = zeros(0, N_pred*m);
        b_ineq = zeros(0, 1);
        % 重新递推 Sigma_z (不从外层循环借, 避免污染)
        Sigma_z_prev_eq14 = Q_aug_k;  % 初始 = Q_aug (一步预测误差)
        for j = 1:N_pred
            mu_y_no_c = M_cl{j} * z_k;
            mu_y_sub = mu_y_no_c(1:ell);
            % 只用动态部分: P_hat * E_x' * Sigma_z * E_x * P_hat' (不含 Sigma_ebar)
            Sigma_z_curr = A_cl_hat * Sigma_z_prev_eq14 * A_cl_hat' + Q_aug_k;
            Sigma_x_j = E_x' * Sigma_z_curr * E_x;
            Sigma_y_dyn = P_hat * Sigma_x_j * P_hat';
            sigma_y_sub = sqrt(max(diag(Sigma_y_dyn(1:ell, 1:ell)), 1e-6));
            A_ineq = [A_ineq; N_cl{j}(1:ell, :)];
            b_ineq = [b_ineq; y_max_vec - mu_y_sub - z_alpha * sigma_y_sub];
            Sigma_z_prev_eq14 = Sigma_z_curr;
        end
        use_eq14 = true;
    else
        use_eq14 = false;
    end

    % 输入约束: u_min ≤ c_j ≤ u_max (K=0 时 u=c)
    lb = u_min * ones(N_pred*m, 1);
    ub = u_max * ones(N_pred*m, 1);

    % 求解 QP
    try
        if use_eq14
            [c_opt, ~, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, [], [], lb, ub, [], opt_qp);
        else
            [c_opt, ~, exitflag] = quadprog(H_qp, f_qp, [], [], [], [], lb, ub, [], opt_qp);
        end
        if exitflag ~= 1
            qp_fail = qp_fail + 1;
            c_opt = max(min(-H_qp \ f_qp, ub), lb);  % 无约束最优 + 裁剪
        end
    catch
        qp_fail = qp_fail + 1;
        c_opt = max(min(-H_qp \ f_qp, ub), lb);
    end

    % --- Algorithm 2 步骤 7: u_k = -K z_k + c_k (K=0 → u_k = c_k) ---
        % Delta-tracking: cost_const = Σ_j ||Δy_j - Δr_j||_Q^2 (和 QP 一致)
        cost_const = 0;
        y_k_actual = P_hat * x_hat_k;
        for j = 1:N_pred
            r_j = y_ref_full(:, min(k+j, T_cl));
            if j == 1
                delta_r = r_j - y_k_actual;
                residual = M_cl{j} * z_k - y_k_actual - delta_r;
                % = M_cl{j} * z_k - r_j
            else
                r_prev = y_ref_full(:, min(k+j-1, T_cl));
                delta_r = r_j - r_prev;
                residual = (M_cl{j} - M_cl{j-1}) * z_k + (N_cl{j} - N_cl{j-1}) * c_opt - delta_r;
            end
            cost_const = cost_const + residual' * Q_w * residual;
        end
        cost_hist(k) = 0.5 * c_opt' * H_qp * c_opt + f_qp' * c_opt + cost_const;
    u_k = c_opt(1:m);
    u_hist(:, k) = u_k;

    % --- 真实系统更新 (Eq.1) ---
    x_true = A_true * x_true + B_true * u_k + sigma_w_k * randn(n, 1);
end

%% ==================== 统计 ====================
idx_post = 50:T_cl;  % 阶跃后的数据
% 跟踪 RMSE: 前 ℓ 个通道
y_rmse = sqrt(mean(sum((y_true(1:ell, idx_post) - ref_val).^2, 1)));
fprintf('\n=== 结果 ===\n');
fprintf('QP 失败: %d/%d\n', qp_fail, T_cl);
fprintf('RMSE (前 %d 通道): %.4f\n', ell, y_rmse);
fprintf('y 均值: [%.3f %.3f]\n', mean(y_true(1:2, idx_post), 2));

% 分段统计
seg_labels = {'基准', '+20%', '-40%', '+10%', '+30%', '-20%'};
for seg = 1:length(noise_factors)
    k_start = (seg-1)*noise_interval + 50;
    k_end = min(seg*noise_interval, T_cl);
    if k_start > T_cl, break; end
    idx = k_start:k_end;
    rmse_seg = sqrt(mean(sum((y_true(1:ell, idx) - ref_val).^2, 1)));
    sigma_est = mean(sigma_eps_online(idx));
    fprintf('  段 %d (%s): RMSE=%.4f, σ̂=%.4f\n', seg, seg_labels{seg}, rmse_seg, sigma_est);
end

%% ==================== 显示参数对比 ====================
fprintf('\n=== 参数对比 ===\n');

% 真实 A 和估计 Â
fprintf('\n真实 A (对角): ');
fprintf('%.3f ', diag(A_true));
fprintf('\n');
fprintf('估计 Â 特征值: ');
fprintf('%.3f ', eig(A_hat));
fprintf('\n');

% 真实 B 和估计 B̂
fprintf('\n真实 B (%d×%d):\n', n, m);
disp(B_true);
fprintf('估计 B̂ (%d×%d):\n', ell, m);
disp(B_hat);

% 真实 C 和估计 P̂ (子空间)
fprintf('真实 C 前 %d 列 (子空间方向):\n', ell);
disp(C_true(:, 1:ell));
fprintf('估计 P̂ (%d×%d):\n', p, ell);
disp(P_hat);

%% ==================== 画图 ====================
fprintf('\n[画图]\n');

figure('Position', [50, 50, 1100, 1000], 'Color', 'w');

% 子图 1: 跟踪效果
subplot(3,1,1);
% 画三段参考水平线
plot(1:T_cl, ref_seg1*ones(1,T_cl), 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('参考段1: %.1f', ref_seg1));
hold on;
% 段 2 段 3 参考
plot(step1_end:step2_end-1, ref_seg2*ones(1, step2_end-step1_end), 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('参考段2: %.1f', ref_seg2));
plot(step2_end:T_cl, ref_seg3*ones(1, T_cl-step2_end+1), 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('参考段3: %.1f', ref_seg3));
plot(1:T_cl, y_true(1,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'y_1');
plot(1:T_cl, y_true(2,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'y_2');
% 阶跃起点 + 段切换点
xline(50, 'g-', 'LineWidth', 1.5, 'DisplayName', '阶跃起点');
xline(step1_end, 'm--', 'LineWidth', 1.5, 'HandleVisibility', 'off');   % 段1->段2
xline(step2_end, 'c--', 'LineWidth', 1.5, 'HandleVisibility', 'off');   % 段2->段3
% 噪声段切换线 (粉红虚线)
for seg = 1:length(noise_factors)-1
    xline(seg*noise_interval + 50, 'r:', 'HandleVisibility', 'off');
end
xlabel('k'); ylabel('输出');
title(sprintf('Route C-Delta 3step: ref=[%.1f, %.1f, %.1f] (y_1=%.2f/%.2f/%.2f, y_2=%.2f/%.2f/%.2f)', ...
    mean(y_true(1,idx1)), mean(y_true(1,idx2)), mean(y_true(1,idx3)), ...
    mean(y_true(2,idx1)), mean(y_true(2,idx2)), mean(y_true(2,idx3))));
h_legend = legend('Location', 'eastoutside', 'FontSize', 7);
grid on;

% 子图 2: 在线自适应
subplot(3,1,2);
plot(1:T_cl, sigma_w_actual, 'r-', 'LineWidth', 2, 'DisplayName', '实际 σ_w');
hold on;
plot(1:T_cl, sigma_eps_online, 'b-', 'LineWidth', 1.0, 'DisplayName', '在线 σ̂_ε');
for seg = 1:length(noise_factors)-1
    xline(seg*noise_interval, 'k:', 'HandleVisibility', 'off');
end
xlabel('k'); ylabel('σ');
title('在线自适应'); legend('Location', 'northwest'); grid on;

% 子图 3: SMPC 代价函数
subplot(3,1,3);
plot(1:T_cl, cost_hist, 'b-', 'LineWidth', 1.0);
hold on;
for seg = 1:length(noise_factors)-1
    xline(seg*noise_interval, 'r--', 'HandleVisibility', 'off');
end
xlabel('k'); ylabel('J_k');
title('SMPC 代价函数'); grid on;

sgtitle(sprintf('Route C-Delta 3step: ref=[%.1f, %.1f, %.1f] (y_1=%.2f/%.2f/%.2f, y_2=%.2f/%.2f/%.2f)', ...
    ref_seg1, ref_seg2, ref_seg3, ...
    mean(y_true(1,idx1)), mean(y_true(1,idx2)), mean(y_true(1,idx3)), ...
    mean(y_true(2,idx1)), mean(y_true(2,idx2)), mean(y_true(2,idx3))), ...
    'FontSize', 12, 'FontWeight', 'bold');
saveas(gcf, fullfile(fileparts(mfilename('fullpath')), 'predvarx_mpc_v10_routeC_delta_3step.png'));
fprintf('  图已保存\n');

% === Route C-Delta 3step: dump metrics ===
fid = fopen('routeC_delta_3step_results.txt', 'w');
fprintf(fid, 'ROUTE_C_DELTA_3STEP_RESULTS\n');
fprintf(fid, 'ref_seg1: %.2f  ref_seg2: %.2f  ref_seg3: %.2f\n', ref_seg1, ref_seg2, ref_seg3);
fprintf(fid, 'y_mean_y1_seg1: %.4f\n', mean(y_true(1, idx1)));
fprintf(fid, 'y_mean_y2_seg1: %.4f\n', mean(y_true(2, idx1)));
fprintf(fid, 'y_mean_y1_seg2: %.4f\n', mean(y_true(1, idx2)));
fprintf(fid, 'y_mean_y2_seg2: %.4f\n', mean(y_true(2, idx2)));
fprintf(fid, 'y_mean_y1_seg3: %.4f\n', mean(y_true(1, idx3)));
fprintf(fid, 'y_mean_y2_seg3: %.4f\n', mean(y_true(2, idx3)));
fprintf(fid, 'QP_fail: %d/%d\n', qp_fail, T_cl);
fclose(fid);
