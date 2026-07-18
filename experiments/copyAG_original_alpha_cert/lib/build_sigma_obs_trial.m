function [Sigma_obs, meta] = build_sigma_obs_trial(mode, o_window, Sigma_n, p, ell, scale_floor)
% BUILD_SIGMA_OBS_TRIAL Opinion 7 experimental modes for online Sigma_obs.
%
% mode:
%   'declared_shape'   - scale * shape(Sigma_n)   [copyAB proxy]
%   'residual_support' - scale * Q_o Q_o'         [projection support]
%   'additive'         - Sigma_n + scale * Q_o Q_o'  [sensor + residual]
%
% scale is estimated from residual Frobenius energy with DOF (p-ell).
% These are engineering trial laws, not proved identities for Cov(o).

if nargin < 6 || isempty(scale_floor), scale_floor = 1e-8; end
mode = lower(char(mode));
meta = struct();
meta.mode = mode;

O = o_window;
if size(O,1) ~= p
    error('build_sigma_obs_trial: o_window must have p rows');
end
L = size(O,2);
Oc = O - mean(O,2);
dof = max((p-ell)*max(L-1,1), 1);
sigma2 = (norm(Oc,'fro')^2) / dof;
sigma2 = max(sigma2, scale_floor);
meta.sigma2 = sigma2;

Sn = (Sigma_n + Sigma_n')/2;
shape_n = Sn / max(trace(Sn)/p, eps);
meta.rank_n = rank(Sn, 1e-10*norm(Sn,2));

% Residual support projector from sample residual span (thin QR)
if L >= 1
    [Q,~] = qr(Oc, 0);
    % Keep at most p-ell columns with nontrivial energy
    s = sqrt(sum(Q.^2,1));
    keep = s > 1e-12;
    Q = Q(:, keep);
    if size(Q,2) > (p-ell)
        Q = Q(:, 1:(p-ell));
    end
    if isempty(Q)
        Q = zeros(p,0);
    end
else
    Q = zeros(p,0);
end
meta.rank_o = size(Q,2);
QoQo = Q*Q';

switch mode
    case 'declared_shape'
        Sigma_obs = sigma2 * shape_n;
        meta.note = 'proxy: residual-scale * declared shape; NOT Cov(o)';
    case 'residual_support'
        Sigma_obs = sigma2 * QoQo + scale_floor * eye(p);
        meta.note = 'trial: residual-support sigma2*Q_oQ_o + floor I';
    case 'additive'
        Sigma_obs = Sn + sigma2 * QoQo;
        meta.note = 'trial: declared Sigma_n + residual-support term';
    otherwise
        error('build_sigma_obs_trial: unknown mode %s', mode);
end
Sigma_obs = (Sigma_obs + Sigma_obs')/2;
meta.trace = trace(Sigma_obs);
meta.rank = rank(Sigma_obs, 1e-10*max(meta.trace,1));
end
