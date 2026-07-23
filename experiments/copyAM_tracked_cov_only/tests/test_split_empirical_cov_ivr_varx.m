function test_split_empirical_cov_ivr_varx
% Focused gates for empirical-total-covariance free dual (no Sigma_n).
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
rng(2107,'twister');

p = 9; m = 2; ell = 5; tracked = [1 2]; T = 700;
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
% Plant may have heteroscedastic sensor noise, but identifier never sees Sigma_n.
scales = [0.025 0.030 0.045 0.080 0.130 0.210 0.330 0.480 0.700]';
Corr = eye(p);
Corr(3,4)=0.35; Corr(4,3)=0.35;
Corr(5,6)=-0.25; Corr(6,5)=-0.25;
D = diag(scales);
Sigma_true = D*Corr*D;  %#ok<NASGU> plant only; not passed to identifier
L = chol(Sigma_true,'lower');
y = y_clean + L*randn(p,T);

E = zeros(p,numel(tracked)); E(tracked,:) = eye(numel(tracked));

%% 1) Signature must not accept / require Sigma_n as an input argument
nargin_id = nargin('split_empirical_cov_ivr_varx');
assert(nargin_id == 5, 'identifier should have arity 5 (y,u,ell,tracked,oblique_alpha)');
src = fileread(fullfile(fileparts(here),'split_empirical_cov_ivr_varx.m'));
sig_line = regexp(src, ...
    'function\s+\[Ahat,Bhat,P,R,Sigma_eps,stats\]\s*=\s*split_empirical_cov_ivr_varx\(([^)]*)\)', ...
    'tokens','once');
assert(~isempty(sig_line), 'could not parse identifier signature');
assert(~contains(sig_line{1},'Sigma_n'), ...
    'function signature must not include Sigma_n input');
assert(contains(src,'empirical-total-covariance metric'), ...
    'stats metric name must be empirical-total-covariance metric');
assert(contains(src,'NOT a sensor-noise optimum') || contains(src,'not a sensor-noise'), ...
    'source must explicitly deny sensor-noise optimum claim');

%% 2) alpha=1: dual, left/right invariance, tracked columns, R'*P=I
[A1,B1,P1,R1,S1,st1] = split_empirical_cov_ivr_varx(y,u,ell,tracked,1);
assert(st1.oblique_alpha == 1);
assert(st1.uses_true_Sigma_n == false);
assert(strcmp(st1.metric_name,'empirical-total-covariance metric'));
assert(norm(R1'*P1-eye(ell),'fro') < 1e-8,'dual R''*P=I failed');
assert(norm(P1*R1'*E-E,'fro') < 1e-8,'right coverage failed');
assert(norm(E'*P1*R1'-E','fro') < 1e-8,'left value preservation failed');
assert(norm(R1(:,1:numel(tracked))-E,'fro') < 1e-12,'tracked extractor columns changed');
assert(isequal(size(A1),[ell ell]) && isequal(size(B1),[ell m]));
assert(min(eig((S1+S1')/2)) > -1e-8,'Sigma_eps not PSD');

%% 3) Empirical free-block objective: alpha=1 not worse than orthogonal baseline
assert(st1.free_emp_cov_objective <= st1.free_emp_cov_baseline + 1e-8, ...
    'empirical free objective worse than orthogonal baseline');
assert(st1.free_emp_cov_improvement >= -1e-10);
% On generic data free dual should typically be oblique under total cov.
assert(st1.free_oblique_norm > 1e-6 || st1.free_emp_cov_improvement < 1e-12, ...
    'unexpected: free dual neither oblique nor metric-tied to orthogonal');

%% 4) Closed-form / KKT optimality under free empirical metric + feasible perturbation
% Reconstruct free dual in free coordinates.
Nperp = null(E');
r = ell-numel(tracked);
q = numel(tracked);
V = Nperp'*P1(:,q+1:end);
W = Nperp'*R1(:,q+1:end);
assert(norm(W'*V-eye(r),'fro') < 1e-8,'free dual W''V=I failed');
C_reg = st1.C_emp + st1.C_emp_ridge*eye(size(st1.C_emp));
objW = trace(W'*C_reg*W);

% KKT: grad tr(W'C W)=2 C W must lie in range(V) when W'V=I is active.
% Equiv: (I - V*pinv(V)) * C_reg * W ~ 0  (or residual after least-squares on V).
CW = C_reg * W;
Lambda = V \ CW;  % least-squares multiplier map
kkt_res = norm(CW - V*Lambda, 'fro') / max(norm(CW,'fro'), eps);
assert(kkt_res < 1e-6, sprintf('KKT residual too large: %.3e', kkt_res));

% Feasible free directions: Wpert = W + Uker*S keeps Wpert'V = I if V'*Uker=0.
Uker = null(V');
assert(size(Uker,2) >= 1,'need nontrivial free ambient dim for perturbation');
n_trials = 8;
max_obj_increase = -inf;
eps_pert = 1e-4;
for t = 1:n_trials
    S = eps_pert*randn(size(Uker,2),r);
    Wpert = W + Uker*S;
    assert(norm(Wpert'*V-eye(r),'fro') < 1e-9);
    objp = trace(Wpert'*C_reg*Wpert);
    max_obj_increase = max(max_obj_increase, objp - objW);
end
% Second-order: objective increase is O(eps^2); first-order decrease forbidden.
assert(max_obj_increase >= -1e-10, 'feasible perturbation decreased objective (unexpected)');
assert(max_obj_increase < 50*(eps_pert^2)*max(objW,1), ...
    sprintf('feasible free perturbation increased empirical objective by %.3e',max_obj_increase));

% Closed-form recompute must match returned free dual at alpha=1
SinvV = C_reg \ V;
Wcf = SinvV / ((V'*SinvV + (V'*SinvV)')/2);
assert(norm(W - Wcf,'fro') < 1e-7,'returned free dual mismatches closed form');

%% 5) alpha=0 collapses free dual toward orthogonal loading (interpolation only)
[~,~,P0,R0,~,st0] = split_empirical_cov_ivr_varx(y,u,ell,tracked,0);
assert(st0.oblique_alpha == 0);
assert(norm(R0(:,q+1:end)-P0(:,q+1:end),'fro') < 1e-10, ...
    'alpha=0 free dual must equal free loading');
assert(norm(R0'*P0-eye(ell),'fro') < 1e-8);

%% 6) alpha=0.5 is convex interpolation of free duals (geometry sanity)
[~,~,~,R05,~,st05] = split_empirical_cov_ivr_varx(y,u,ell,tracked,0.5);
W05 = Nperp'*R05(:,q+1:end);
W0 = Nperp'*R0(:,q+1:end);
assert(norm(W05 - (0.5*W0 + 0.5*W),'fro') < 1e-7, ...
    'alpha=0.5 free dual is not midpoint interpolation');
assert(st05.free_emp_cov_objective + 1e-8 >= min(st0.free_emp_cov_objective, st1.free_emp_cov_objective) - 1e-6);

fprintf('PASS empirical free dual: dual=%.2e left=%.2e right=%.2e free-oblique=%.3e; emp-obj %.6g -> %.6g (gain %.3e); pert max dJ=%.2e; no Sigma_n\n', ...
    st1.dual_error, st1.tracked_left_error, st1.tracked_right_error, st1.free_oblique_norm, ...
    st1.free_emp_cov_baseline, st1.free_emp_cov_objective, st1.free_emp_cov_improvement, ...
    max_obj_increase);
end
