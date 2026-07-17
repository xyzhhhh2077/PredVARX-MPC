function test_opinion02_noise_opt
% Opinion 2: declared Sigma_n free-block noise objective + isotropic collapse.
% - free_noise_objective / baseline / improvement must be present
% - free_noise_is_declared_Sigma_n must be true
% - isotropic Sigma_n => ||R-P|| < 1e-8 and isotropic_collapsed_to_R_eq_P
% - heteroscedastic Sigma_n => free_noise_improvement > 0
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
rng(20260718,'twister');

p = 9; m = 2; ell = 5; tracked = [1 2]; T = 650;
F = diag([0.91 0.77 0.58]);
G = [0.40 -0.15; 0.18 0.32; -0.08 0.24];
C = randn(p,3);
C(1,:) = [1 0 0];
C(2,:) = [0 1 0];
u = randn(m,T);
x = zeros(3,T+1);
y_clean = zeros(p,T);
for k = 1:T
    y_clean(:,k) = C*x(:,k);
    x(:,k+1) = F*x(:,k) + G*u(:,k) + 0.035*randn(3,1);
end

%% Case 1: isotropic declared Sigma_n -> R collapses to P
sigma = 0.06;
Sigma_iso = sigma^2*eye(p);
y_iso = y_clean + sigma*randn(p,T);
[~,~,Pi,Ri,~,sti] = split_control_free_ivr_varx(y_iso,u,ell,tracked,Sigma_iso);

assert(isfield(sti,'free_noise_baseline'),'missing free_noise_baseline');
assert(isfield(sti,'free_noise_objective'),'missing free_noise_objective');
assert(isfield(sti,'free_noise_improvement'),'missing free_noise_improvement');
assert(isfield(sti,'free_noise_is_declared_Sigma_n'), ...
    'missing free_noise_is_declared_Sigma_n');
assert(isfield(sti,'isotropic_collapsed_to_R_eq_P'), ...
    'missing isotropic_collapsed_to_R_eq_P');
assert(sti.free_noise_is_declared_Sigma_n == true, ...
    'free_noise objective must be declared as Sigma_n metric');
assert(norm(Ri-Pi,'fro') < 1e-8, ...
    sprintf('isotropic must ||R-P||<1e-8, got %.3e', norm(Ri-Pi,'fro')));
assert(sti.isotropic_collapsed_to_R_eq_P == true, ...
    'isotropic_collapsed_to_R_eq_P must be true for Sigma_n ~ I');
assert(sti.free_noise_improvement >= -1e-10, ...
    'isotropic free-noise objective regressed vs R=P baseline');

%% Case 2: heteroscedastic/correlated declared Sigma_n -> strict improvement
scales = [0.025 0.030 0.045 0.080 0.130 0.210 0.330 0.480 0.700]';
Corr = eye(p);
Corr(3,4)=0.35; Corr(4,3)=0.35;
Corr(5,6)=-0.25; Corr(6,5)=-0.25;
D = diag(scales);
Sigma_het = D*Corr*D;
L = chol(Sigma_het,'lower');
y_het = y_clean + L*randn(p,T);
[~,~,Ph,Rh,~,sth] = split_control_free_ivr_varx(y_het,u,ell,tracked,Sigma_het);

assert(sth.free_noise_is_declared_Sigma_n == true, ...
    'heteroscedastic free_noise must still be declared Sigma_n metric');
assert(sth.free_noise_improvement > 0, ...
    sprintf('heteroscedastic must free_noise_improvement>0, got %.3e', ...
    sth.free_noise_improvement));
assert(sth.free_noise_objective <= sth.free_noise_baseline + 1e-10, ...
    'optimized free dual worse than orthogonal free baseline');
assert(norm(Rh-Ph,'fro') > 1e-4, ...
    'heteroscedastic free dual should be genuinely oblique (R~=P)');
assert(sth.isotropic_collapsed_to_R_eq_P == false, ...
    'heteroscedastic case must not claim isotropic R=P collapse');

fprintf(['PASS opinion02: iso ||R-P||=%.3e collapsed=%d; ' ...
    'het improvement=%.3e collapsed=%d; declared_Sigma_n=%d/%d\n'], ...
    norm(Ri-Pi,'fro'), sti.isotropic_collapsed_to_R_eq_P, ...
    sth.free_noise_improvement, sth.isotropic_collapsed_to_R_eq_P, ...
    sti.free_noise_is_declared_Sigma_n, sth.free_noise_is_declared_Sigma_n);
end
