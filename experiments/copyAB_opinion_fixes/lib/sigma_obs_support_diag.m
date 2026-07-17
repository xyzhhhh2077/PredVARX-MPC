function [rank_o, rank_n, trace_o, trace_on, note] = sigma_obs_support_diag(o_window, Sigma_n_declared, p, ell)
%SIGMA_OBS_SUPPORT_DIAG  Rank/trace diagnostics for Sigma_obs support (opinion 7).
%
%   [rank_o, rank_n, trace_o, trace_on, note] = sigma_obs_support_diag( ...
%       o_window, Sigma_n_declared, p, ell)
%
% Inputs
%   o_window         p x L matrix of projected residuals
%                    o = (I - P*R')*(y - y_mean). Support is at most p-ell.
%   Sigma_n_declared p x p declared full-space sensor-noise covariance
%                    (or a scaled shape of it). Typically PD with rank p.
%   p                output dimension
%   ell              latent dimension
%
% Outputs
%   rank_o   numerical rank of empirical Cov(o) from the residual window
%   rank_n   numerical rank of declared Sigma_n
%   trace_o  trace of empirical Cov(o)
%   trace_on trace of the online SMPC proxy
%            Sigma_on = max(sigma_hat^2, 1e-8) * shape(Sigma_n),
%            with sigma_hat from residual Frobenius scale using denom (p-ell)
%   note     short string: objects are NOT identical; proxy is not Cov(o)
%
% Theory note (opinion 7)
%   Projected residual covariance should live on the unresolved subspace
%   (rank ~ p-ell), e.g. Sigma_o = Q_o * Sigma_xi * Q_o'. Declared sensor
%   noise Sigma_n is a full-space object (rank ~ p). The online agent
%   Sigma_obs used in SMPC is a scaled declared shape: a sensor-noise-floor
%   proxy, NOT the identity Cov(o). Rank/trace comparison is diagnostic only.

if nargin < 4
    error('sigma_obs_support_diag:Args', ...
        'Usage: sigma_obs_support_diag(o_window, Sigma_n_declared, p, ell)');
end

o_window = double(o_window);
Sigma_n_declared = double(Sigma_n_declared);
p = double(p);
ell = double(ell);

if isempty(o_window)
    error('sigma_obs_support_diag:EmptyO', 'o_window must be nonempty.');
end
if size(o_window, 1) ~= p
    error('sigma_obs_support_diag:DimO', ...
        'o_window must be p x L (got %d x %d, p=%g).', size(o_window,1), size(o_window,2), p);
end
if ~isequal(size(Sigma_n_declared), [p p])
    error('sigma_obs_support_diag:DimN', ...
        'Sigma_n_declared must be p x p.');
end
if ~(isscalar(ell) && ell >= 0 && ell < p)
    error('sigma_obs_support_diag:Ell', ...
        'ell must satisfy 0 <= ell < p.');
end

L = size(o_window, 2);
O = o_window - mean(o_window, 2);
den_cov = max(L - 1, 1);
Sigma_o = (O * O') / den_cov;
Sigma_o = (Sigma_o + Sigma_o') / 2;

Sigma_n = (Sigma_n_declared + Sigma_n_declared') / 2;

% Numerical ranks (relative eigenvalue tolerance)
tol_o = max(size(Sigma_o)) * eps(norm(Sigma_o, 2));
tol_n = max(size(Sigma_n)) * eps(norm(Sigma_n, 2));
if ~isfinite(tol_o) || tol_o <= 0
    tol_o = eps;
end
if ~isfinite(tol_n) || tol_n <= 0
    tol_n = eps;
end
eo = sort(max(real(eig(Sigma_o)), 0), 'descend');
en = sort(max(real(eig(Sigma_n)), 0), 'descend');
rank_o = sum(eo > tol_o);
rank_n = sum(en > tol_n);

trace_o = real(trace(Sigma_o));

% Online proxy matching copyAB/copyAA runner convention:
%   sigma_hat^2 from residual energy with residual DOF (p-ell),
%   shape from full-space declared Sigma_n (trace-normalized to mean 1).
den_scale = max((p - ell) * (L - 1), 1);
sigma_hat = norm(O, 'fro') / sqrt(den_scale);
tr_n = real(trace(Sigma_n));
if tr_n <= 0 || ~isfinite(tr_n)
    shape_n = eye(p);
else
    shape_n = Sigma_n / (tr_n / p);
end
Sigma_on = max(sigma_hat^2, 1e-8) * shape_n;
Sigma_on = (Sigma_on + Sigma_on') / 2;
trace_on = real(trace(Sigma_on));

note = sprintf([ ...
    'Sigma_obs online proxy is scaled declared Sigma_n (rank_n=%d, tr_on=%.4g); ' ...
    'NOT Cov(o) (rank_o=%d <= p-ell=%d, tr_o=%.4g). ' ...
    'Objects differ: residual support vs full-space sensor floor.'], ...
    rank_n, trace_on, rank_o, p - ell, trace_o);
end
