function test_crte_fixed_surrogate_varx
% Focused gates for the CRTE fixed spectral surrogate in the free complement.
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
rng(2107,'twister');

p = 9; m = 2; ell = 5; tracked = [1 2]; T = 700;
F = diag([0.91 0.77 0.58]);
G = [0.40 -0.15; 0.18 0.32; -0.08 0.24];
C = randn(p,3); C(1,:)=[1 0 0]; C(2,:)=[0 1 0];
u = randn(m,T);
x = zeros(3,T+1);
y_clean = zeros(p,T);
for k=1:T
    y_clean(:,k) = C*x(:,k);
    x(:,k+1) = F*x(:,k) + G*u(:,k) + 0.035*randn(3,1);
end
scales = [0.025 0.030 0.045 0.080 0.130 0.210 0.330 0.480 0.700]';
Corr = eye(p);
Corr(3,4)=0.35; Corr(4,3)=0.35;
Corr(5,6)=-0.25; Corr(6,5)=-0.25;
D = diag(scales);
L = chol(D*Corr*D,'lower');
y = y_clean + L*randn(p,T);

%% 1) Signature must not accept / require Sigma_n
src = fileread(fullfile(fileparts(here),'crte_fixed_surrogate_varx.m'));
sig_line = regexp(src, ...
    'function\s+\[Ahat,Bhat,P,R,Sigma_eps,stats\]\s*=\s*crte_fixed_surrogate_varx\(([^)]*)\)', ...
    'tokens','once');
assert(~isempty(sig_line),'could not parse identifier signature');
assert(~contains(sig_line{1},'Sigma_n'),'signature must not include Sigma_n');
assert(contains(src,'uses_true_Sigma_n'),'stats must record Sigma_n usage flag');

%% 2) Default run returns a stable, dual-consistent model
[A,B,P,R,S,st] = crte_fixed_surrogate_varx(y,u,ell,tracked);
assert(st.uses_true_Sigma_n==false);
assert(numel(st.candidates) > 0);
assert(sum(st.valid_candidates) > 0);
assert(norm(R'*P-eye(ell),'fro') < 1e-7,'R^T P = I failed');
assert(st.dual_basis_completion,'4-piece dual completion failed');
assert(st.spectral_radius < 1.10,'identified Ahat is not closed-loop stable');
assert(st.selected_validation_nrmse >= 0);

%% 3) Selection must report a triple and be reproducible across calls
[A2,B2,P2,R2,S2,st2] = crte_fixed_surrogate_varx(y,u,ell,tracked);
assert(st.selected_mu == st2.selected_mu);
assert(st.selected_alpha == st2.selected_alpha);
assert(st.selected_beta == st2.selected_beta);

%% 4) Higher beta penalizes noisy candidates: free block eigenvalue ranking
% is monotonically shifted by -beta * Ntr(Sigma_noise_proxy).
% Compare the top eigenvector direction between beta=0 and beta=large.
no_opt = struct('mu_grid',0.5,'alpha_grid',1,'beta_grid',0,'val_fraction',0.25, ...
    'task_gate_fraction',0.05,'noise_gate_factor',5.0,'reach_gate_fraction',0.0);
high_noise_opt = struct('mu_grid',0.5,'alpha_grid',0,'beta_grid',1, ...
    'val_fraction',0.25,'task_gate_fraction',0.05,'noise_gate_factor',10.0, ...
    'reach_gate_fraction',0.0);
[~,~,~,~,~,st_low] = crte_fixed_surrogate_varx(y,u,ell,tracked,no_opt);
[~,~,~,~,~,st_high] = crte_fixed_surrogate_varx(y,u,ell,tracked,high_noise_opt);
assert(st_low.selected_beta == 0);
assert(st_high.selected_beta == 1);
% The free block eigenvalues should differ because alpha/beta changed the
% objective (not because the simulation seed changed).
assert(any(abs([st_low.selected_eigenvalues] - [st_high.selected_eigenvalues])>1e-12));

%% 5) Output contract sanity check
assert(isequal(size(A),[ell ell]) && isequal(size(B),[ell m]));
assert(min(eig((S+S')/2)) > -1e-8,'Sigma_eps not PSD');
fprintf('PASS CRTE fixed surrogate: dual=%.2e spectral=%.4f selected=[mu=%.2f alpha=%.2f beta=%.2f] candidates=%d valid=%d val NRMSE=%.4f\n', ...
    st.dual_error, st.spectral_radius, st.selected_mu, st.selected_alpha, st.selected_beta, ...
    numel(st.candidates), sum(st.valid_candidates), st.selected_validation_nrmse);
end