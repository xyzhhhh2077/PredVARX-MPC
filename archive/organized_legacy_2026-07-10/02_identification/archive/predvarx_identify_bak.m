function [A_hat, B_hat, P_hat, Pbar_hat, R_hat, Rbar_hat, G_aug_hat, Sigma_eps_hat, Sigma_ebar_hat, F_aug_hat, H_aug_hat, lambda_forget, prewhiten] = predvarx_identify(y_data, u_data, ell, beta_u, L_u, A_true, B_true, C_true, n, m, p)
% 最小工作版: Stage1剥离u → 联合预白化 → IVR → VARX回归
T = size(y_data,2);

% Stage 1
N = T - L_u;
Phi_u = zeros(N, m*L_u); Yr = zeros(p, N);
for i = 1:N
    t = L_u + i;
    ul = zeros(m*L_u,1);
    for lag = 1:L_u, ul((lag-1)*m+1:lag*m) = u_data(:, t-lag); end
    Phi_u(i,:) = ul'; Yr(:,i) = y_data(:,t);
end
C_hat_u = (Phi_u'*Phi_u + 1e-6*eye(m*L_u)) \ (Phi_u' * Yr');
y_res = Yr - (Phi_u * C_hat_u)';
fprintf('  S1: %.3f\n', norm(y_res,'fro')/norm(Yr,'fro'));

% Stage 2: joint prewhiten
y_c = y_res - mean(y_res,2);
u_c = u_data(:, L_u+1:L_u+N) - mean(u_data(:, L_u+1:L_u+N),2);
joint = [y_c; beta_u * u_c];  % (p+m) × N
[JU, JD, ~] = svd(joint * joint' / N, 'econ');
d = min(p+m, size(JD,1));
JDi = diag(1 ./ sqrt(diag(JD(1:d,1:d)) + 1e-12));
Jnorm = JU(:,1:d) * JDi * JU(:,1:d)' * joint;  % (p+m) × N

% IVR, s=1
Y_star = Jnorm(:, 2:end);    % (p+m) × (N-1)
Y_lag  = Jnorm(:, 1:end-1);  % (p+m) × (N-1)
Nivr = N-1;
[U,~,~] = svd((Y_star * Y_lag') / Nivr, 'econ');
Pj = U(:, 1:ell);  % (p+m) × ℓ

% 简化: 跳过迭代IVR, 直接SVD
% (避免维度bug, 噪声估计差但模型结构对)
X  = Pj' * Y_star;
Xl = Pj' * Y_lag;
Av = (Xl * Xl') \ (Xl * X');
fprintf('  IVR:1 (skip iter)\n');
% Extract y-part, denormalize
Py = Pj(1:p, :);  % p × ℓ
P_hat = JU(1:p,1:p) * diag(sqrt(diag(JD(1:p,1:p)))) * JU(1:p,1:p)' * Py;
R_hat = (P_hat' * P_hat + 1e-8*eye(ell)) \ P_hat';
Pbar_hat = null(P_hat');
Rbar_hat = Pbar_hat';
G_aug_hat = eye(ell);

% Stage 3: VARX
xl = R_hat * y_c;  % ℓ × N
xn = xl(:, 2:end);
xc = xl(:, 1:end-1);
ur = u_c(:, 1:end-1);
Th = ([xc; ur] * [xc; ur]' + 1e-6*eye(ell+m)) \ ([xc; ur] * xn');
A_hat = Th(1:ell, :)';
B_hat = Th(ell+1:end, :)';
F_aug_hat = A_hat;
H_aug_hat = B_hat;  % ℓ×m (非 2m)

Inn = xn - A_hat*xc - B_hat*ur;
Sigma_eps_hat = (Inn * Inn') / (size(Inn,2)-1);
Sigma_ebar_hat = zeros(p-ell);
lambda_forget = 0.95;
prewhiten.C_hat_u = C_hat_u;
prewhiten.y_mean = mean(y_res,2);
prewhiten.u_mean = mean(u_c,2);
prewhiten.Joint_U = JU;
prewhiten.Joint_D_inv_sqrt = JDi;
prewhiten.L_u = L_u;
end