%% test_copyT_process_lv_smpc_static
% Static consistency checks for copyT process LV-SMPC setup.
clear; clc; addpath(fileparts(fileparts(mfilename('fullpath'))));
rng(7,'twister');
p = 30; n = 6; m = 3; ell = 5; tracked = [1 2]; T = 400;
A = diag([0.94 0.88 0.78 0.64 0.50 0.35]);
B = randn(n,m)*0.2;
C = randn(p,n)*0.3; C(1,:) = [1 0 0.1 0 0 0]; C(2,:) = [0 1 0 -0.1 0 0];
x = zeros(n,T+1); y = zeros(p,T); u = randn(m,T);
for k = 1:T
    y(:,k) = C*x(:,k) + 0.03*randn(p,1);
    x(:,k+1) = A*x(:,k) + B*u(:,k) + 0.02*randn(n,1);
end
[Ahat,Bhat,P,R,Sigma_eps,stats] = control_aware_subspace_varx(y,u,ell,tracked);
assert(all(size(Ahat)==[ell ell]), 'Ahat size wrong');
assert(all(size(Bhat)==[ell m]), 'Bhat size wrong');
assert(all(size(P)==[p ell]), 'P size wrong');
assert(norm(R-P,'fro') < 1e-12, 'copyT uses orthogonal control-aware variant, so R=P expected');
assert(stats.tracked_projection_error < 1e-10, 'tracked axes not exactly covered');
assert(min(eig((Sigma_eps+Sigma_eps')/2)) > -1e-8, 'Sigma_eps not PSD within tolerance');
fprintf('PASS test_copyT_process_lv_smpc_static: projErr=%.3e recon=%.3f\n', stats.tracked_projection_error, stats.reconstruction_residual);
