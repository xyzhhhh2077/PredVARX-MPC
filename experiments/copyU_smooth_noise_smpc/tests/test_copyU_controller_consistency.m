%% test_copyU_controller_consistency
% Regression test: all N*nu inputs must be bounded, mean prediction must use
% centered inputs, and out.cost must equal the raw MPC objective.
clear; clc; addpath(fileparts(fileparts(mfilename('fullpath'))));

ny = 2; ell = 2; nu = 2; N = 3;
model.A = [0.82 0.08; 0 0.74];
model.B = [0.45 -0.10; 0.15 0.35];
model.P = eye(ny); model.R = eye(ny);
model.y_mean = [0.30; -0.20];
model.u_mean = [0.55; -0.40];
model.Sigma_eps = 1e-3*eye(ell);
model.Sigma_obs = 1e-3*eye(ny);

opt.N = N;
opt.Q = diag([4 3]);
opt.Ru = 0.2*eye(nu);
opt.u_min = -0.25;
opt.u_max = 0.25;
opt.H = eye(ny);
opt.h = 10*ones(ny,1);
opt.alpha_joint = 0.10;

y = [1.1; -0.7];
r = [2.0; 1.5];
[z, y_pred, U, out] = centered_smpc_step(y,r,model,opt);

assert(numel(U) == N*nu, 'U dimension mismatch');
assert(all(U >= opt.u_min-1e-9 & U <= opt.u_max+1e-9), ...
    'Every element of the N*nu control plan must satisfy input bounds');

U0 = repmat(model.u_mean,N,1);
M1 = model.P*model.A;
G1 = zeros(ny,N*nu);
G1(:,1:nu) = model.P*model.B;
y_pred_expected = model.y_mean + M1*z + G1*(U-U0);
assert(norm(y_pred-y_pred_expected,inf) < 1e-10, ...
    'Reported y_pred must use centered input U-U0');

J_raw = 0;
for j = 1:N
    Mj = model.P*(model.A^j);
    Gj = zeros(ny,N*nu);
    for i = 0:j-1
        Gj(:,i*nu+1:(i+1)*nu) = model.P*(model.A^(j-1-i))*model.B;
    end
    ej = model.y_mean + Mj*z + Gj*(U-U0) - r;
    J_raw = J_raw + ej'*opt.Q*ej;
end
DU = U-U0;
J_raw = J_raw + DU'*kron(eye(N),opt.Ru)*DU;
assert(abs(out.cost-J_raw) < 1e-8*max(1,abs(J_raw)), ...
    'out.cost must equal the raw centered MPC objective');

assert(isfield(out,'lb') && isfield(out,'ub'), 'Controller must expose bounds for audit');
assert(numel(out.lb)==N*nu && numel(out.ub)==N*nu, ...
    'Audit bounds must cover all N*nu decision variables');

fprintf('PASS test_copyU_controller_consistency: max|U|=%.4f J=%.6f\n',max(abs(U)),out.cost);
