function [A_hat, B_hat, P_hat, Pbar_hat, R_hat, Rbar_hat, G_aug_hat, Sigma_eps_hat, Sigma_ebar_hat, F_aug_hat, H_aug_hat, lambda_forget, prewhiten] = predvarx_identify_v2(y_data, u_data, ell, beta_u, L_u, A_true, B_true, C_true, n, m, p)
% predvarx_identify_V2: 剥离u → SVD降维 → VARX回归 (最小可工作版)
T = size(y_data,2);
N = T - L_u;

% Stage 1: strip u
Ph = zeros(N, m*L_u); Yr = zeros(p, N);
for i = 1:N
    t = L_u + i; ul = zeros(m*L_u,1);
    for lag = 1:L_u, ul((lag-1)*m+1:lag*m) = u_data(:,t-lag); end
    Ph(i,:) = ul'; Yr(:,i) = y_data(:,t);
end
Cu = (Ph'*Ph + 1e-6*eye(m*L_u)) \ (Ph' * Yr');
yr = Yr - (Ph * Cu)';  % p × N, u的线性影响已剥离

% Stage 2: joint prewhiten + SVD
y_c = yr - mean(yr,2);
u_c = u_data(:,L_u+1:L_u+N) - mean(u_data(:,L_u+1:L_u+N),2);
[Uy,Sy,~] = svd(y_c*y_c'/N, 'econ');
P_hat = Uy(:, 1:ell);  % p × ℓ
R_hat = (P_hat'*P_hat + 1e-8*eye(ell)) \ P_hat';  % ℓ × p

% Stage 3: VARX regression
xl = R_hat * y_c;  % ℓ × N
xn = xl(:, 2:end);  % ℓ × (N-1)
xc = xl(:, 1:end-1);
ur = u_c(:, 1:end-1);
Phi = [xc; ur];
Th = (Phi*Phi' + 1e-6*eye(ell+m)) \ (Phi * xn');
A_hat = Th(1:ell, :)';
B_hat = Th(ell+1:end, :)';

% Noise
Inn = xn - A_hat*xc - B_hat*ur;
Sigma_eps_hat = (Inn * Inn') / max(size(Inn,2)-1, 1);
Sigma_ebar_hat = zeros(p-ell);

Pbar_hat = null(P_hat');
Rbar_hat = Pbar_hat';
G_aug_hat = eye(ell);
F_aug_hat = A_hat;
H_aug_hat = B_hat;  % ℓ×m
lambda_forget = 0.95;

prewhiten.C_hat_u = Cu;
prewhiten.y_mean = mean(yr,2);
prewhiten.u_mean = mean(u_c,2);
prewhiten.Joint_U = Uy;
prewhiten.Joint_D_inv_sqrt = diag(1./sqrt(diag(Sy)+1e-12));
prewhiten.L_u = L_u;
fprintf('  V2: P=%dx%d A=%dx%d B=%dx%d\\n',size(P_hat,1),size(P_hat,2),size(A_hat,1),size(A_hat,2),size(B_hat,1),size(B_hat,2));
end