function test_opinion05_sigma_denoms
% Opinion 5: primary Sigma_eps uses conditional zero-mean Gaussian ML.
% OLS residual-DOF remains a diagnostic; T-2 is not a required field.

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
rng(1705,'twister');

p = 9; m = 2; ell = 5; tracked = [1 2]; T = 400;
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
sigma = 0.06;
Sigma_n = sigma^2*eye(p);
y = y_clean + sigma*randn(p,T);

[Ahat,Bhat,P,R,Sigma_eps,stats] = split_control_free_ivr_varx( ...
    y,u,ell,tracked,Sigma_n);

% --- required multi-denom fields ---
req = {'Sigma_eps','Sigma_eps_ml','Sigma_eps_ols', ...
       'N_residual','Sigma_eps_denom_primary', ...
       'Sigma_eps_denom_ml','Sigma_eps_denom_ols'};
for i = 1:numel(req)
    assert(isfield(stats,req{i}),'missing stats.%s',req{i});
end

N = stats.N_residual;
assert(N == T-1,'N_residual should be T-1 (one lag lost in VARX)');
assert(stats.Sigma_eps_denom_primary == max(N,1),'primary ML denom mismatch');
assert(stats.Sigma_eps_denom_ml == max(N,1),'ML denom mismatch');
assert(stats.Sigma_eps_denom_ols == max(N-(ell+m),1),'OLS denom mismatch');

% --- shapes ---
assert(isequal(size(Sigma_eps),[ell ell]));
assert(isequal(size(stats.Sigma_eps),[ell ell]));
assert(isequal(size(stats.Sigma_eps_ml),[ell ell]));
assert(isequal(size(stats.Sigma_eps_ols),[ell ell]));

% --- primary return is the theory-backed ML scale ---
assert(norm(Sigma_eps-stats.Sigma_eps,'fro') < 1e-14, ...
    'return Sigma_eps must equal stats.Sigma_eps (primary ML scale)');

% Reconstruct Gram from each denom and check consistency
G_primary = stats.Sigma_eps * stats.Sigma_eps_denom_primary;
G_ml  = stats.Sigma_eps_ml  * stats.Sigma_eps_denom_ml;
G_ols = stats.Sigma_eps_ols * stats.Sigma_eps_denom_ols;
assert(norm(G_primary-G_ml,'fro') < 1e-10*max(norm(G_primary,'fro'),1), ...
    'primary and ML Gram mismatch');
assert(norm(G_primary-G_ols,'fro') < 1e-10*max(norm(G_primary,'fro'),1), ...
    'primary and OLS Gram mismatch');

% Explicit scale relations (positive denoms)
assert(norm(stats.Sigma_eps_ml - stats.Sigma_eps,'fro') < 1e-12, ...
    'primary vs ML scale relation failed');
assert(norm(stats.Sigma_eps_ols - stats.Sigma_eps*N/max(N-(ell+m),1),'fro') < 1e-12, ...
    'OLS vs primary scale relation failed');

% --- PSD (symmetric + eigenvalues >= -tol) ---
check_psd(Sigma_eps,'return Sigma_eps');
check_psd(stats.Sigma_eps,'stats.Sigma_eps');
check_psd(stats.Sigma_eps_ml,'stats.Sigma_eps_ml');
check_psd(stats.Sigma_eps_ols,'stats.Sigma_eps_ols');

% Downstream contract: primary return remains usable as model.Sigma_eps
assert(min(eig((Sigma_eps+Sigma_eps')/2)) > -1e-8, ...
    'default Sigma_eps not usable as PSD covariance');
assert(isequal(size(Ahat),[ell ell]) && isequal(size(Bhat),[ell m]));
assert(norm(R'*P-eye(ell),'fro') < 1e-8,'dual identity broken by Sigma change');

fprintf(['PASS opinion05 sigma denoms: N=%d; denoms primary/ml/ols=%d/%d/%d; ' ...
    'eigmin primary/ml/ols=%.3e/%.3e/%.3e\n'], ...
    N, stats.Sigma_eps_denom_primary, stats.Sigma_eps_denom_ml, stats.Sigma_eps_denom_ols, ...
    min(eig((stats.Sigma_eps+stats.Sigma_eps')/2)), ...
    min(eig((stats.Sigma_eps_ml+stats.Sigma_eps_ml')/2)), ...
    min(eig((stats.Sigma_eps_ols+stats.Sigma_eps_ols')/2)));
end

function check_psd(S,name)
Ssym = (S+S')/2;
assert(norm(S-Ssym,'fro') < 1e-12, '%s not symmetric', name);
assert(min(eig(Ssym)) > -1e-8, '%s is not PSD (min eig)', name);
end
