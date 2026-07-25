function [A_hat, B_hat, P_hat, Pbar_hat, R_hat, Rbar_hat, G_aug_hat, Sigma_eps_hat, Sigma_ebar_hat, F_aug_hat, H_aug_hat, lambda_forget] = predvarx_identify(y_data, u_data, ell, beta_u, L_u, A_true, B_true, C_true, n, m, p)
% PREDVARX_IDENTIFY 离线/在线 PredVARX 辨识 (§3 Algorithm 1)
%   输入:
%     y_data: p × T 数据
%     u_data: m × T 数据
%     ell:   降维维度
%     beta_u: IVR 中 Z_u 对 u 的权重 (Route C.1 用)
%     L_u:   u 滞后阶数
%     A_true, B_true, C_true: 真值 (用于 cos 角度日志, 可传 [])
%     n, m, p: 维度
%   输出:
%     A_hat, B_hat: 降维动力学 (ℓ × ℓ, ℓ × m)
%     P_hat, Pbar_hat: 观测子空间与互补 (p × ℓ, p × (p-ℓ))
%     R_hat, Rbar_hat: 提取器与互补 (ℓ × p, (p-ℓ) × p)
%     G_aug_hat: 噪声映射 (默认 eye(ℓ))
%     Sigma_eps_hat, Sigma_ebar_hat: 动态 / 静态噪声协方差
%     F_aug_hat, H_aug_hat: 增广矩阵 (s=1 时 F_aug = Â, H_aug = B̂)
%     lambda_forget: 默认 0.95

T_off = size(y_data, 2);
fprintf('\n=== PredVARX 辨识 (T_off=%d) ===\n', T_off);

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
    u_lag = []; for lag = 1:L_u, u_lag = [u_lag; u_data(:, t-lag)]; end
    Phi_u(i, :) = u_lag';      % [u_{k-1}^T, u_{k-2}^T]
    Y_strip(:, i) = y_data(:, t);  % y_k
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
u_c = u_data(:, L_u+1:L_u+N_strip) - mean(u_data(:, L_u+1:L_u+N_strip), 2);  % m × N_strip (中心化 u)
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
u_past_for_iv = u_data(:, L_u+1:L_u+N_ivr)';   % u_{k-1}^T, N_ivr × m (转置!)
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
X_hat_all = R_hat' * y_data;  % ℓ × T_off (降维状态序列)
t_start = 2;
N_reg = T_off - t_start;

% 构造回归矩阵 Φ = [x̂_k, u_k]
Phi_ab = zeros(N_reg, ell + m);
X_target = zeros(N_reg, ell);
for i = 1:N_reg
    t = t_start + i - 1;
    Phi_ab(i, :) = [X_hat_all(:, t); u_data(:, t)]';  % [x̂_k^T, u_k^T]
    X_target(i, :) = X_hat_all(:, t+1)';              % x̂_{k+1}^T
end

% Ridge 回归: Θ̂ = (Φ^T Φ + λI)^{-1} Φ^T X_target
Theta_ab = (Phi_ab'*Phi_ab + 1e-4*eye(ell+m)) \ (Phi_ab' * X_target);
A_hat = Theta_ab(1:ell, :)';      % ℓ × ℓ (降维动力学矩阵 Â)
B_hat = Theta_ab(ell+1:end, :)';  % ℓ × m (降维控制矩阵 B̂)

% --- 阶段 4: 噪声估计 (§2.4 文档) ---
X_pred_off = Phi_ab * Theta_ab;
Sigma_eps_hat = ((X_target - X_pred_off)' * (X_target - X_pred_off)) / N_reg;

resid_y = y_data - P_hat * X_hat_all;
Sigma_resid = (resid_y * resid_y') / T_off;
Sigma_ebar_hat = Rbar_hat' * Sigma_resid * Rbar_hat;

fprintf('  Â eig: [%.3f, %.3f]\n', eig(A_hat));
fprintf('  Σ̂_ε  diag: [%.4f %.4f]\n', diag(Sigma_eps_hat));

G_aug_hat = eye(ell);  % 默认
F_aug_hat = A_hat;
H_aug_hat = B_hat;

lambda_forget = 0.95;

end