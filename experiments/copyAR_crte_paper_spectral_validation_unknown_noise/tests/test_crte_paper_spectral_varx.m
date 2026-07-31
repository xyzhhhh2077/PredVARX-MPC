function test_crte_paper_spectral_varx
% Focused gates: paper spectral surrogate + residual unknown-noise + validation.
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
Sigma_n = D*Corr*D; Sigma_n=(Sigma_n+Sigma_n')/2;
L = chol(Sigma_n,'lower');
y = y_clean + L*randn(p,T);

%% 1) Signature must NOT accept Sigma_n
src = fileread(fullfile(fileparts(here),'crte_paper_spectral_varx.m'));
sig_line = regexp(src, ...
    'function\s+\[Ahat,Bhat,P,R,Sigma_eps,stats\]\s*=\s*crte_paper_spectral_varx\(([^)]*)\)', ...
    'tokens','once');
assert(~isempty(sig_line),'could not parse identifier signature');
assert(~contains(sig_line{1},'Sigma_n'),'signature must not include Sigma_n');
assert(contains(src,'paper_ntr'),'must implement paper Ntr');
assert(contains(src,'trv/d'),'paper Ntr must normalize trace by matrix dimension');
assert(contains(src,'epsilon_ntr = 10e-6'),'paper Ntr epsilon must be 10e-6');
assert(contains(src,'residual'),'must document residual proxy');
assert(contains(src,'validation_nrmse'),'must use validation selection');
assert(~contains(lower(src),'j_teacher'),'must NOT use profiled min-teacher');

runner_src = fileread(fullfile(fileparts(here), ...
    'copyAR_crte_paper_spectral_validation_unknown_noise.m'));
assert(contains(runner_src,'mu_grid = [0.10 0.25 0.50 0.75]'), ...
    'final copyAR grid must contain interior mu candidates only');
assert(contains(runner_src,'stats.selected_mu>0 && stats.selected_mu<1'), ...
    'final copyAR runner must enforce an interior selected mu');

%% 2) Default run
[A,B,P,R,S,st] = crte_paper_spectral_varx(y,u,ell,tracked);
assert(st.uses_true_Sigma_n==false);
assert(strcmp(st.ntr_mode,'paper_trace_normalize'));
assert(numel(st.candidates) > 0);
assert(sum(st.valid_candidates) > 0);
assert(norm(R'*P-eye(ell),'fro') < 1e-7,'R^T P = I failed');
assert(st.dual_basis_completion,'4-piece dual completion failed');
assert(st.spectral_radius < 1.10,'identified Ahat not stable enough');
assert(st.selected_validation_nrmse >= 0);
assert(st.selected_alpha==1 && st.selected_beta==1);

%% 3) Reproducible selection
[~,~,~,~,~,st2] = crte_paper_spectral_varx(y,u,ell,tracked);
assert(st.selected_mu == st2.selected_mu);
assert(abs(st.selected_validation_nrmse - st2.selected_validation_nrmse) < 1e-12);

%% 4) Output contract
assert(isequal(size(A),[ell ell]) && isequal(size(B),[ell m]));
assert(min(eig((S+S')/2)) > -1e-8,'Sigma_eps not PSD');
fprintf(['PASS paper CRTE spectral+val+unknownSn: dual=%.2e spectral=%.4f ' ...
    'selected=[mu=%.2f alpha=%.2f beta=%.2f] candidates=%d valid=%d val NRMSE=%.4f\n'], ...
    st.dual_error, st.spectral_radius, st.selected_mu, st.selected_alpha, st.selected_beta, ...
    numel(st.candidates), sum(st.valid_candidates), st.selected_validation_nrmse);
end
