function [A_hat, B_hat, P_hat, Pbar_hat, R_hat, Rbar_hat, G_aug_hat, Sigma_eps_hat, Sigma_ebar_hat, F_aug_hat, H_aug_hat, lambda_forget, prewhiten] = predvarx_identify_oblique(y_data, u_data, ell, beta_u, L_u, A_true, B_true, C_true, n, m, p)
% PredVARX with oblique dual bases.  The IVR-estimated loading/extraction
% pair is aligned so R_hat' * P_hat = I instead of forcing R_hat=P_hat.
% Unused compatibility inputs are kept to match predvarx_identify().

T = size(y_data, 2);
N = T - L_u;

% Stage 1: regress out delayed input effects before IVR.
Phi_u = zeros(N, m * L_u);
Y_strip = zeros(N, p);
for i = 1:N
    t = L_u + i;
    phi = zeros(m * L_u, 1);
    for lag = 1:L_u
        phi((lag-1)*m+1:lag*m) = u_data(:, t-lag);
    end
    Phi_u(i,:) = phi';
    Y_strip(i,:) = y_data(:,t)';
end
C_hat_ols = (Phi_u' * Phi_u + 1e-6 * eye(m*L_u)) \ (Phi_u' * Y_strip);
yr = (Y_strip - Phi_u * C_hat_ols)';

% Stage 2: center and prewhiten the residual observation sequence.
yr_mean = mean(yr, 2);
yr_c = yr - yr_mean;
Sigma_yr = (yr_c * yr_c') / N;
[U_d, D_d] = eig((Sigma_yr + Sigma_yr')/2);
[D_sorted, idx] = sort(real(diag(D_d)), 'descend');
U_d = real(U_d(:,idx));
D_sqrt = diag(sqrt(max(D_sorted, 1e-12)));
D_sqrt_inv = diag(1 ./ sqrt(max(D_sorted, 1e-12)));
Y_star = D_sqrt_inv * U_d' * yr_c;

% Stage 3: IVR initialization and alternating predictable-subspace update.
Nm1 = N - 1;
Y_star_lag = Y_star(:,1:Nm1);
Y_star_cur = Y_star(:,2:N);
Pi_instr = Y_star_lag / (Y_star_lag' * Y_star_lag + 1e-6 * eye(Nm1)) * Y_star_lag';
M_init = Pi_instr * (Y_star_cur * Y_star_cur') / N;
C_star = top_eigvecs(M_init, ell);
trace_prev = 0;
trace_curr = 0;
max_iter = 30;
for ivr_iter = 1:max_iter
    X_hat = Y_star_cur' * C_star;
    X_lag = Y_star_lag' * C_star;
    A_var = (X_lag' * X_lag + 1e-8 * eye(ell)) \ (X_lag' * X_hat);
    X_pred = X_lag * A_var;
    Sigma_x_hat = (X_pred' * X_pred) / Nm1;
    trace_curr = real(trace(Sigma_x_hat));
    if ivr_iter > 1 && abs(trace_curr-trace_prev) < 1e-4 * max(abs(trace_curr),1)
        break;
    end
    trace_prev = trace_curr;
    W = Y_star_cur * X_pred;
    M_IVR = W / (X_pred' * X_pred + 1e-6 * eye(ell)) * W';
    C_star = top_eigvecs(M_IVR, ell);
end

% Stage 4: recover an oblique dual pair and align it by SVD.
P_raw = U_d * D_sqrt * C_star;
R_raw = U_d * D_sqrt_inv * C_star;
S = R_raw' * P_raw;
[Us, Ss, Vs] = svd(S, 'econ');
P_hat = P_raw * Vs / Ss;
R_hat = R_raw * Us;

% Complete the dual bases.  null(R') guarantees R'Pbar=0.  Since P is
% generally oblique to Pbar, construct the remaining dual rows from the
% inverse of the full square basis [P Pbar].
Pbar_hat = null(R_hat');
Qfull = [P_hat, Pbar_hat];
if rcond(Qfull) < 1e-12
    error('Oblique basis [P Pbar] is numerically singular.');
end
Dual = inv(Qfull)';
R_hat = Dual(:,1:ell);
Rbar_hat = Dual(:,ell+1:end);

% Stage 5: refit a first-order VARX model in the oblique coordinates.
y_c = y_data - mean(y_data, 2);
xl = R_hat' * y_c;
xn = xl(:,2:end);
xc = xl(:,1:end-1);
u_c = u_data - mean(u_data, 2);
ur = u_c(:,1:end-1);
Phi_ab = [xc; ur];
Theta = (Phi_ab * Phi_ab' + 1e-6 * eye(ell+m)) \ (Phi_ab * xn');
A_hat = Theta(1:ell,:)';
B_hat = Theta(ell+1:end,:)';

% Stage 6: dynamic and complementary observation residual covariances.
Inn_w = xn - A_hat * xc - B_hat * ur;
Sigma_eps_hat = (Inn_w * Inn_w') / max(size(Inn_w,2)-1,1);
Sigma_eps_hat = (Sigma_eps_hat + Sigma_eps_hat')/2;
y_hat_reconstruct = P_hat * xl;
resid_obs = y_c - y_hat_reconstruct;
ebar_proj = Rbar_hat' * resid_obs;
Sigma_ebar_hat = (ebar_proj * ebar_proj') / max(T-1,1);
Sigma_ebar_hat = (Sigma_ebar_hat + Sigma_ebar_hat')/2;

F_aug_hat = A_hat;
H_aug_hat = B_hat;
G_aug_hat = eye(ell);
lambda_forget = 0.95;
prewhiten.C_hat_u = C_hat_ols;
prewhiten.y_mean = yr_mean;
prewhiten.u_mean = mean(u_data,2);
prewhiten.Joint_U = U_d;
prewhiten.Joint_D_sqrt_inv = D_sqrt_inv;
prewhiten.L_u = L_u;

fprintf('  Oblique IVR: P=%dx%d A=%dx%d B=%dx%d iter=%d trace=%.4f dual_err=%.1e\n', ...
    size(P_hat,1),size(P_hat,2),size(A_hat,1),size(A_hat,2),size(B_hat,1),size(B_hat,2), ...
    ivr_iter,trace_curr,norm(R_hat'*P_hat-eye(ell),'fro'));
end

function V = top_eigvecs(M, ell)
M = (M + M')/2;
[U,D] = eig(M);
[~,idx] = sort(real(diag(D)), 'descend');
V = real(U(:,idx(1:ell)));
end
