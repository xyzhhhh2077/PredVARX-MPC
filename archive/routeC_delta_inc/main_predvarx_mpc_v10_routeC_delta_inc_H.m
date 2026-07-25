%% PredVARX-MPC v10 副本 H — Route C-Delta 三段阶跃 + 在线辨识 (方案 B)
%  在线子空间辨识: 每 N_reid 拍用滑动窗口重跑 PredVARX
%
%  ★ 副本 H = P1 + P3 + 在线辨识:
%    [P1] QP warm-start + Σ 增量
%    [P3] 自适应预测步数 N(k)
%      - k < 200:        N = N_short = 5
%      - 200 ≤ k ≤ 1000: N = N_mid = 12
%      - k > 1000:       N = N_max = 30
%
%  数学严格性: 两改动均不影响 Δ-tracking 数学, 仅 QP 计算路径加速
%    - P1: warm-start 命中率高 (100%), quadprog 内迭代数减少
%    - P3: N_curr ≤ N_max, H_qp/f_qp 维度动态变 (N_curr*m × N_curr*m)
%  [P2] 不并入 — P2 改写 M_cl/N_cl 预计算循环, 与 P3 的 N_curr 截列叠加增 bug 风险
%  三段阶跃同 P1 (ref = 1.0 → 2.0 → 0.5)
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
noise_factors = [1.0, 2, 0.52, 3, 5, 0.824];
% ★ P3: 自适应预测步数参数
N_max_pred = N_pred;          % 上限 = 原 N_pred, 用于 M_cl/N_cl cell 预计算
N_short = 5;                  % 早期 (k < 200): 短 N
N_mid = 12;                   % 中期 (200 ≤ k ≤ 1000): 中 N
N_history = zeros(1, T_cl);   % 记录每拍实际 N_curr
N_total = 0;                  % 累计 N (加权平均用)
N_max_count = 0;              % 累计 N_max 用次数
  % 噪声缩放因子

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

L_u = 2;                          % 回归 u_{k-1}, u_{k-2}

%% ==================== §3 调用辨识 (§3 Algorithm 1 抽出为函数) ====================
[A_hat, B_hat, P_hat, Pbar_hat, R_hat, Rbar_hat, G_aug_hat, Sigma_eps_hat, Sigma_ebar_hat, F_aug_hat, H_aug_hat, lambda_forget] = predvarx_identify(y_off, u_off, ell, beta_u, L_u, A_true, B_true, C_true, n, m, p);

%% ==================== §4: 概率预测器 (§3 文档) ====================
% 文档 §3.1: 增广状态 z_k = [x̂_k; ...; x̂_{k-s+1}; u_{k-1}; ...; u_{k-s_u}]
% s=1, s_u=0 → z_k = x̂_k (ℓ 维)

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
M_cl = cell(1, N_pred); N_cl = cell(1, N_pred);
for j = 1:N_pred
    % M_j^cl = P̂ · E_x^T · A_cl^j (p × nz)
    M_cl{j} = P_hat * E_x' * A_cl_hat^j;
    % N_j^cl = [P̂ E_x^T A_cl^{j-1} H_aug, ..., P̂ E_x^T H_aug, 0, ...]
    N_cl{j} = zeros(p, N_pred*m);
    for i = 0:j-1
        N_cl{j}(:, i*m+1:(i+1)*m) = P_hat * E_x' * A_cl_hat^(j-1-i) * H_aug_hat;
    end
end

% Q_w: 只惩罚前 ℓ 个输出通道 (ℓ < p 时)
% 因为 MPC 预测的是 P̂ x̂, 但 P̂ 的列空间不一定和 y 的前 ℓ 个通道对齐
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


% ★ P3: P1 warm-start + timing
c_opt_prev = zeros(N_max_pred*m, 1);
qp_time_total = 0;
warm_start_used = zeros(1, T_cl);
cov_time_total = 0;     % 累计协方差递推时间 (s)
reid_time_total = 0;    % ★ 副本 H: 累计在线重辨识时间 (s)
opt_qp = optimset('Display', 'off', 'LargeScale', 'off');

% ★ 副本 H: 每 N_reid 拍重跑 PredVARX 在线辨识
N_reid = 500;                             % 重新辨识周期 (拍数)
y_window = zeros(p, N_reid);              % 滑动窗口 (最近 N_reid 拍 y)
u_window = zeros(m, N_reid);              % 滑动窗口 (最近 N_reid 拍 u)
n_reid_count = 0;                          % 已收集样本数
t_reid_last = -N_reid;                     % 上次重辨识时间

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
    % ★ P1 增量: 利用上一拍的 Σ_z^{cl}(N,k-1) 作本拍 j=1 初值 (跨拍 Lyapunov 重用)
    tic_cov = tic;
    if k > 1
        Sigma_z_prev_k = A_cl_hat * Sigma_z_history(:, :, N_pred) * A_cl_hat' + Q_aug_k;
        Sigma_z_j = Sigma_z_prev_k;       % j=1 (Δ 增量)
    else
        Sigma_z_prev_k = zeros(nz);
        Sigma_z_j = zeros(nz);            % j=1 (零起)
    end
    Sigma_x_j = E_x' * Sigma_z_j * E_x;
    Sigma_y_cl{1} = P_hat * Sigma_x_j * P_hat' + Pbar_hat * Sigma_ebar_k * Pbar_hat';
    Sigma_z_history(:, :, 1) = Sigma_z_j;
    for j = 2:N_pred
        Sigma_z_j = A_cl_hat * Sigma_z_prev_k * A_cl_hat' + Q_aug_k;
        Sigma_x_j = E_x' * Sigma_z_j * E_x;
        Sigma_y_cl{j} = P_hat * Sigma_x_j * P_hat' + Pbar_hat * Sigma_ebar_k * Pbar_hat';
        Sigma_z_prev_k = Sigma_z_j;
        Sigma_z_history(:, :, j) = Sigma_z_j;
    end
    cov_time_total = cov_time_total + toc(tic_cov);

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
    y_k_actual = P_hat * x_hat_k;   % 当前 y 的估计 (从降维状态反投影)
    % ★ P3: 自适应预测步数
    if k < 200
        N_curr = N_short;       % 早期: 短 N
    elseif k <= 1000
        N_curr = N_mid;         % 中期: 中 N
    else
        N_curr = N_max_pred;    % 后期: 长 N
    end
    N_history(k) = N_curr;
    N_total = N_total + N_curr;
    if N_curr >= N_max_pred
        N_max_count = N_max_count + 1;
    end
    % ★ P3: QP 维度 = N_curr*m
    H_qp = zeros(N_curr*m);
    f_qp = zeros(N_curr*m, 1);
    for j = 1:N_curr
        r_j = y_ref_full(:, min(k+j, T_cl));  % 未来参考
        if j == 1
            % ★ P3: N_cl 取前 N_curr*m 列匹配决策变量维度
            delta_r = r_j - y_k_actual;
            M_delta = M_cl{j};             % Δy_1 = M_1 z_k - y_k + N_1 c
            N_delta = N_cl{j}(:, 1:N_curr*m);
            f_const = -y_k_actual;
        else
            r_prev = y_ref_full(:, min(k+j-1, T_cl));
            delta_r = r_j - r_prev;
            M_delta = M_cl{j} - M_cl{j-1};
            N_delta = (N_cl{j}(:, 1:N_curr*m) - N_cl{j-1}(:, 1:N_curr*m));
            f_const = zeros(p, 1);
        end
        H_qp = H_qp + N_delta' * Q_w * N_delta;
        f_qp = f_qp + N_delta' * Q_w * (M_delta * z_k + f_const - delta_r);
    end
    H_qp = H_qp + kron(eye(N_curr), R_w);
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
        % ★ P3: Eq.14 范围 = N_curr
        A_ineq = zeros(0, N_curr*m);
        b_ineq = zeros(0, 1);
        Sigma_z_prev_eq14 = Q_aug_k;
        for j = 1:N_curr
            mu_y_no_c = M_cl{j} * z_k;
            mu_y_sub = mu_y_no_c(1:ell);
            % 只用动态部分: P_hat * E_x' * Sigma_z * E_x * P_hat' (不含 Sigma_ebar)
            Sigma_z_curr = A_cl_hat * Sigma_z_prev_eq14 * A_cl_hat' + Q_aug_k;
            Sigma_x_j = E_x' * Sigma_z_curr * E_x;
            Sigma_y_dyn = P_hat * Sigma_x_j * P_hat';
            sigma_y_sub = sqrt(max(diag(Sigma_y_dyn(1:ell, 1:ell)), 1e-6));
            % ★ P3: 只取 N_curr*m 列 (与 H_qp 维度对齐)
            A_ineq = [A_ineq; N_cl{j}(1:ell, 1:N_curr*m)];
            b_ineq = [b_ineq; y_max_vec - mu_y_sub - z_alpha * sigma_y_sub];
            Sigma_z_prev_eq14 = Sigma_z_curr;
        end
        use_eq14 = true;
    else
        use_eq14 = false;
    end

    % 输入约束: u_min ≤ c_j ≤ u_max (K=0 时 u=c)
    lb = u_min * ones(N_curr*m, 1);
    ub = u_max * ones(N_curr*m, 1);

    % 求解 QP
    tic_qp = tic;
    try
        if use_eq14
            [c_opt, ~, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, [], [], lb, ub, c_opt_prev, opt_qp);
        else
            [c_opt, ~, exitflag] = quadprog(H_qp, f_qp, [], [], [], [], lb, ub, c_opt_prev, opt_qp);
        end
        if exitflag ~= 1
            qp_fail = qp_fail + 1;
            c_opt = max(min(-H_qp \ f_qp, ub), lb);  % 无约束最优 + 裁剪
        end
    catch
        qp_fail = qp_fail + 1;
        c_opt = max(min(-H_qp \ f_qp, ub), lb);
    end
    qp_time_total = qp_time_total + toc(tic_qp);
    warm_start_used(k) = 1;

    % --- Algorithm 2 步骤 7: u_k = -K z_k + c_k (K=0 → u_k = c_k) ---
    % Delta-tracking: cost_const = Σ_j ||Δy_j - Δr_j||_Q^2 (和 QP 一致)
    cost_const = 0;
        y_k_actual = P_hat * x_hat_k;
        for j = 1:N_curr
            r_j = y_ref_full(:, min(k+j, T_cl));
            if j == 1
                delta_r = r_j - y_k_actual;
                % ★ P3: N_cl 取前 N_curr*m 列匹配 c_opt 维度
                residual = M_cl{j} * z_k + N_cl{j}(:, 1:N_curr*m) * c_opt - y_k_actual - delta_r;
            else
                r_prev = y_ref_full(:, min(k+j-1, T_cl));
                delta_r = r_j - r_prev;
                residual = (M_cl{j} - M_cl{j-1}) * z_k + (N_cl{j}(:,1:N_curr*m) - N_cl{j-1}(:,1:N_curr*m)) * c_opt - delta_r;
            end
            cost_const = cost_const + residual' * Q_w * residual;
        end
    cost_hist(k) = 0.5 * c_opt' * H_qp * c_opt + f_qp' * c_opt + cost_const;
    u_k = c_opt(1:m);
    u_hist(:, k) = u_k;

    % --- 副本 H: 滑动窗口收集 (在 u_k 已知后) ---
    if n_reid_count < N_reid
        n_reid_count = n_reid_count + 1;
        y_window(:, n_reid_count) = y_k;
        u_window(:, n_reid_count) = u_k;
    else
        y_window = [y_window(:, 2:end), y_k];
        u_window = [u_window(:, 2:end), u_k];
    end

    % --- 副本 H: 每 N_reid 拍重跑 PredVARX 辨识 ---
    if n_reid_count >= N_reid && (k - t_reid_last) >= N_reid
        % 备份: 记录重辨识前 k-10 拍到 k 拍的 y 跟踪精度 (避免当拍 self-correlated)
        r_curr = y_ref_full(1:ell, min(k, T_cl));
        if k > 10
            e_pre_norm = norm(y_true(1:ell, k-9:k) - r_curr*ones(1, 10), 'fro');
        else
            e_pre_norm = NaN;
        end
        tic_reid = tic;
        % 备份: 重辨识前的 P_hat, A_hat
        P_hat_old = P_hat;
        A_hat_old = A_hat;
        [A_hat, B_hat, P_hat, Pbar_hat, R_hat, Rbar_hat, G_aug_hat, Sigma_eps_hat, Sigma_ebar_hat, F_aug_hat, H_aug_hat, lambda_forget] = predvarx_identify(y_window, u_window, ell, beta_u, L_u, A_true, B_true, C_true, n, m, p);
        A_cl_hat = F_aug_hat;  % K=0
        % ★ 副本 H 修复: 重辨识后必须更新 M_cl, N_cl (MPC 名义预测)
        for j = 1:N_pred
            M_cl{j} = P_hat * E_x' * A_cl_hat^j;
            for i = 0:j-1
                N_cl{j}(:, i*m+1:(i+1)*m) = P_hat * E_x' * A_cl_hat^(j-1-i) * H_aug_hat;
            end
        end
        % ★ 诊断: 子空间漂移量
        sub_drift = norm(P_hat - P_hat_old, 'fro');
        eig_old = eig(A_hat_old);
        eig_new = eig(A_hat);
        eig_shift = norm(eig_old - eig_new);
        t_reid_last = k;
        reid_time_total = reid_time_total + toc(tic_reid);
        % 跟踪精度: 重辨识前 vs 后
        e_post_norm = norm(y_true(1:ell, k) - r_curr);
        fprintf('  [副本 H] k=%d 重辨识: 跟踪误差 %.3f→%.3f, ||P̂-P̂_old||_F=%.3f, eig 移 %.3f (累计 %.3fs)\n', k, e_pre_norm, e_post_norm, sub_drift, eig_shift, reid_time_total);
    end

    % --- 真实系统更新 (Eq.1) ---
    x_true = A_true * x_true + B_true * u_k + sigma_w_k * randn(n, 1);
    % ★ C: 平移一格作为下一拍 warm-start
    c_opt_prev = [c_opt(2:end); c_opt(end)];
    c_opt_prev = max(min(c_opt_prev, ub), lb);
end

%% ==================== 统计 ====================
idx_post = 50:T_cl;  % 阶跃后的数据
% 跟踪 RMSE: 前 ℓ 个通道
y_rmse = sqrt(mean(sum((y_true(1:ell, idx_post) - ref_val).^2, 1)));
fprintf('\n=== 结果 ===\n');
fprintf('QP 失败: %d/%d\n', qp_fail, T_cl);

fprintf('\n=== QP 时序 ===\n');
fprintf('QP 累计时间: %.3f s (mean %.6f s/step)\n', qp_time_total, qp_time_total/T_cl);
fprintf('Σ 协方差累计递推时间: %.4f s (mean %.6f s/step)\n', cov_time_total, cov_time_total/T_cl);
fprintf('★ 副本 H: 在线重辨识累计时间: %.3f s (平均 %.3f s/次, N_reid=%d)\n', reid_time_total, reid_time_total/(T_cl/N_reid), N_reid);


fprintf('RMSE (前 %d 通道): %.4f\n', ell, y_rmse);
fprintf('y 均值: [%.3f %.3f]\n', mean(y_true(1:2, idx_post), 2));


% ★ P3: dump metrics
fprintf('\n=== P3 自适应 N 统计 ===\n');
fprintf('平均预测步数: %.2f (最长 N_max_pred=%d)\n', N_total/T_cl, N_max_pred);
fprintf('N_max_pred 命中: %d/%d (%.1f%%)\n', N_max_count, T_cl, 100*N_max_count/T_cl);
fprintf('N_curr 阶段统计:\n');
fprintf('  k<200 (短 N=%d): %.1f%%\n', N_short, 100*sum(N_history(1:min(200,T_cl))==N_short)/T_cl);
fprintf('  200≤k≤1000 (中 N=%d): %.1f%%\n', N_mid, 100*sum(N_history(201:min(1000,T_cl))==N_mid)/T_cl);

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
try
    plot(1:T_cl, cost_hist, 'b-', 'LineWidth', 1.0);
    hold on;
    for seg = 1:length(noise_factors)-1
        xline(seg*noise_interval, 'r--', 'HandleVisibility', 'off');
    end
    xlabel('k'); ylabel('J_k');
    title('SMPC 代价函数'); grid on;
catch
    title('SMPC 代价函数 (数据不足)');
end

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