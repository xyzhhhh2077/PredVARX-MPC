function [A_hat, B_hat, P_hat, Pbar_hat, R_hat, Rbar_hat, G_aug_hat, Sigma_eps_hat, Sigma_ebar_hat, F_aug_hat, H_aug_hat, lambda_forget, prewhiten] = predvarx_identify(y_data, u_data, ell, beta_u, L_u, A_true, B_true, C_true, n, m, p)
% predvarx_identify: IVR (Mo Algorithm 1) + OLS u-stripping + VARX regression
%   严格对齐理论文档 §2.5: 阶段1(OLS剥离u) → 阶段2(IVR迭代精炼) → 阶段3(VARX) → 阶段4(噪声)
%
% Outputs:
%   A_hat      ℓ×ℓ  降维动力学 (used as Akh for state prediction)
%   B_hat      ℓ×m  降维输入矩阵 (used as Bk for state prediction)
%   P_hat      p×ℓ  DLV子空间 (used as Pk for prediction matrices)
%   Pbar_hat   p×(p-ℓ)  静态噪声基底 (null(P_hat'))
%   R_hat      ℓ×p  提取矩阵 (used as Rk: xk = Rk'*yk)
%   Rbar_hat   (p-ell)×p  静态提取
%   G_aug_hat  ℓ×ℓ  增广G (eye(ℓ) for 1st-order)
%   Sigma_eps_hat   ℓ×ℓ  过程噪声协方差
%   Sigma_ebar_hat  (p-ell)×(p-ell)  静态噪声协方差
%   F_aug_hat  ℓ×ℓ  增广动力学 (= A_hat for 1st-order)
%   H_aug_hat  ℓ×m  增广输入 (= B_hat for 1st-order)
%   lambda_forget  1×1  遗忘因子 (0.95)
%   prewhiten      struct  预白化/剥离参数

T = size(y_data, 2);
N = T - L_u;

%% ═══════════════════════════════════════════════════════════════
%% 阶段 1: OLS 剥离 u (理论文档 §2.1)
%%   理论维度: Φ_u ∈ R^{N×mL}, Y ∈ R^{N×p} → C_hat ∈ R^{mL×p}
%%   残差: yr ∈ R^{p×N} (每列是一个残差向量)
%% ═══════════════════════════════════════════════════════════════
Phi_u = zeros(N, m * L_u);    % N × mL  (每行 = [u_{k-1}', u_{k-2}', ...])
Y_strip = zeros(N, p);         % N × p   (每行 = y_k')

for i = 1:N
    t = L_u + i;
    phi = zeros(m * L_u, 1);
    for lag = 1:L_u
        phi((lag-1)*m + 1 : lag*m) = u_data(:, t - lag);
    end
    Phi_u(i, :) = phi';
    Y_strip(i, :) = y_data(:, t)';
end

% OLS: min ||Y - Φ_u * C_hat||^2 → C_hat = (Φ_u'Φ_u)^{-1} Φ_u' Y
C_hat_ols = (Phi_u' * Phi_u + 1e-6 * eye(m * L_u)) \ (Phi_u' * Y_strip);
% C_hat_ols ∈ R^{mL × p}

% 残差: yr = (Y - Φ_u * C_hat_ols)' → p × N
yr = (Y_strip - Phi_u * C_hat_ols)';

%% ═══════════════════════════════════════════════════════════════
%% 阶段 2: IVR 迭代精炼 (理论文档 §2.2, Mo Algorithm 1 步骤 1-6)
%%   输入: yr ∈ R^{p×N} (u 剥离后的残差)
%%   输出: Ĉ ∈ R^{p×ℓ}, Ĝ ∈ R^{p×ℓ} (Ĝ Ĉ = I_ℓ)
%% ═══════════════════════════════════════════════════════════════

% --- Step 1: 中心化 (Mo Input) ---
yr_mean = mean(yr, 2);               % p × 1
yr_c = yr - yr_mean;                 % p × N, 零均值

% --- Step 2: 预白化归一化 (Eq.27-28) ---
%   Σ̃_y = UDU^T, Y* = UD^{-1/2} (协方差白化为 ≈ I)
Sigma_yr = (yr_c * yr_c') / N;      % p × p
[U_d, D_d] = eig(Sigma_yr);
[D_sorted, idx] = sort(diag(D_d), 'descend');
U_d = U_d(:, idx);
D_sqrt_inv = diag(1 ./ sqrt(max(D_sorted, 1e-12)));  % p × p
Y_star = D_sqrt_inv * U_d' * yr_c;  % p × N, 预白化后

% --- Step 3: 工具变量初始化 (Eq.32) ---
Nm1 = N - 1;
Y_star_lag = Y_star(:, 1:Nm1);      % p × (N-1)  工具变量 (滞后)
Y_star_cur = Y_star(:, 2:N);        % p × (N-1)  目标

%   Π = Y_instr * inv(Y_instr'Y_instr + εI) * Y_instr'  (p×p)
%   这是投影到 Y_instr 列空间的 p×p 矩阵
Pi_instr = Y_star_lag / (Y_star_lag' * Y_star_lag + 1e-6 * eye(Nm1)) * Y_star_lag';  % p × p

%   M_init = Π * (Y_cur * Y_cur') / N  (p×p)
M_init = Pi_instr * (Y_star_cur * Y_star_cur') / N;  % p × p

%   Ĉ*_0 = Top-ℓ 特征向量 of M_init
[U_init, D_init] = eig(M_init);
[d_init, idx_init] = sort(diag(D_init), 'descend');
C_star = U_init(:, idx_init(1:ell));  % p × ℓ, 初始 Ĉ*

% 收敛监控
trace_prev = 0;
max_iter = 30;
converge_tol = 1e-4;  % 相对变化判据
ivr_iter = 1;

% --- Step 4-5: 交替迭代精炼 (Eq.33) ---
for ivr_iter = 1:max_iter
    % (3a) 提取降维状态: X̂ = Y*_cur^T * Ĉ*
    X_hat = Y_star_cur' * C_star;         % (N-1) × ℓ

    % (3b) VAR 回归: Â = (X̂_lag^T X̂_lag)^{-1} X̂_lag^T X̂
    X_lag = Y_star_lag' * C_star;         % (N-1) × ℓ
    A_var = (X_lag' * X_lag + 1e-8 * eye(ell)) \ (X_lag' * X_hat);  % ℓ × ℓ

    % (3c) 预测: X̂_pred = X̂_lag * Â
    X_pred = X_lag * A_var;               % (N-1) × ℓ

    % (3d) 收敛判据: trace(Σ̂_x̂) = trace(X̂_pred' X̂_pred / (N-1))
    Sigma_x_hat = (X_pred' * X_pred) / Nm1;  % ℓ × ℓ
    trace_curr = trace(Sigma_x_hat);

    if ivr_iter > 1 && abs(trace_curr - trace_prev) < converge_tol * max(abs(trace_curr), 1)
        break;
    end
    trace_prev = trace_curr;

    % (3e) 用 X_pred 作新工具变量, 重投影更新 Ĉ*
    %   内存优化: W = Y_cur * X_pred (p×ℓ), 避免构造 (N-1)×(N-1) Π
    %   M_IVR = W * inv(X_pred' X_pred + εI) * W'  (p×p)
    W = Y_star_cur * X_pred;                              % p × ℓ
    M_IVR = W / (X_pred' * X_pred + 1e-6 * eye(ell)) * W';  % p × p

    [U_ivr, D_ivr] = eig(M_IVR);
    [d_ivr, idx_ivr] = sort(diag(D_ivr), 'descend');
    C_star = U_ivr(:, idx_ivr(1:ell));  % 更新 Ĉ* (p × ℓ)
end

% --- Step 6: 反归一化 (Eq.18) ---
%   Ĉ = UD^{1/2} Ĉ*,   Ĝ = UD^{-1/2} Ĉ*
P_hat_raw = U_d * diag(sqrt(max(D_sorted, 1e-12))) * C_star;  % p × ℓ  (Ĉ)
R_hat_raw = U_d * D_sqrt_inv * C_star;                         % p × ℓ  (Ĝ)

% QR 列正交化: 确保 P_hat 列正交 → P_hat' P_hat = I_ℓ → R_hat = P_hat
[P_hat, ~] = qr(P_hat_raw, 0);
R_hat = P_hat;  % P_hat 列正交 → P_hat' P_hat = I → R_hat = P_hat

% 验证
proj_err = norm(R_hat' * P_hat - eye(ell), 'fro');
if proj_err > 1e-6
    warning('IVR: R_hat'' * P_hat ≠ I_ℓ (err=%.2e)', proj_err);
end

%% ═══════════════════════════════════════════════════════════════
%% 阶段 3: VARX 回归 (理论文档 §2.3, Mo 步骤 9-12)
%%   ★ 从原始 y 提取降维状态 (不是残差 yr!)
%%   x̂_k = Ĝ * y_k = R_hat' * y_data
%% ═══════════════════════════════════════════════════════════════
y_c = y_data - mean(y_data, 2);   % p × T, 中心化原始 y
xl = R_hat' * y_c;                % ℓ × T, 降维状态

% 回归: x_{k+1} = A_hat * x_k + B_hat * u_k
xn = xl(:, 2:end);                % ℓ × (T-1)  目标 (k+1)
xc = xl(:, 1:end-1);             % ℓ × (T-1)  预测器 x_k
u_c = u_data - mean(u_data, 2);  % m × T, 中心化输入 ★
ur = u_c(:, 1:end-1);            % m × (T-1)  预测器 u_k

Phi_ab = [xc; ur];                % (ℓ+m) × (T-1)
Theta = (Phi_ab * Phi_ab' + 1e-6 * eye(ell + m)) \ (Phi_ab * xn');  % (ℓ+m) × ℓ

A_hat = Theta(1:ell, :)';         % ℓ × ℓ
B_hat = Theta(ell+1:end, :)';     % ℓ × m

%% ═══════════════════════════════════════════════════════════════
%% 阶段 4: 噪声估计 (理论文档 §2.4)
%% ═══════════════════════════════════════════════════════════════

% 过程噪声 Σ_w: DLV 方向残差
Inn_w = xn - A_hat * xc - B_hat * ur;    % ℓ × (T-1)
Sigma_eps_hat = (Inn_w * Inn_w') / max(size(Inn_w, 2) - 1, 1);
Sigma_eps_hat = (Sigma_eps_hat + Sigma_eps_hat') / 2;  % 确保对称

% 静态噪声 Σ_ebar: P⊥ 方向观测残差
y_hat_reconstruct = P_hat * xl;            % p × T, DLV 重建
resid_obs = y_c - y_hat_reconstruct;      % p × T, 观测残差
Pbar_hat = null(P_hat');                  % p × (p-ℓ), 互补子空间
ebar_proj = Pbar_hat' * resid_obs;        % (p-ℓ) × T
Sigma_ebar_hat = (ebar_proj * ebar_proj') / max(T - 1, 1);
Sigma_ebar_hat = (Sigma_ebar_hat + Sigma_ebar_hat') / 2;

% 互补提取
Rbar_hat = null(R_hat');

%% ═══════════════════════════════════════════════════════════════
%% 增广矩阵 (1 阶, s=1, s_u=1)
%% ═══════════════════════════════════════════════════════════════
F_aug_hat = A_hat;    % ℓ×ℓ  (增广动力学, 1阶)
H_aug_hat = B_hat;    % ℓ×m  (增广输入, 1阶)
G_aug_hat = eye(ell); % ℓ×ℓ  (噪声映射, 1阶)
lambda_forget = 0.95;

%% 预白化/剥离参数 (供在线阶段使用)
prewhiten.C_hat_u = C_hat_ols;     % mL × p
prewhiten.y_mean = yr_mean;        % p × 1  (yr 的均值)
prewhiten.u_mean = mean(u_data, 2); % m × 1
prewhiten.Joint_U = U_d;           % p × p
prewhiten.Joint_D_sqrt_inv = D_sqrt_inv;  % p × p
prewhiten.L_u = L_u;

fprintf('  IVR: P=%dx%d A=%dx%d B=%dx%d iter=%d trace=%.4f proj_err=%.1e\n', ...
    size(P_hat,1), size(P_hat,2), size(A_hat,1), size(A_hat,2), ...
    size(B_hat,1), size(B_hat,2), ivr_iter, trace_curr, proj_err);
end
