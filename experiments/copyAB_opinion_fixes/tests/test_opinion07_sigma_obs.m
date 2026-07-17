function test_opinion07_sigma_obs
%TEST_OPINION07_SIGMA_OBS  Rank/trace support diagnostics for Sigma_obs.
%
% Opinion 7: projected residual Cov(o) lives on a (p-ell)-dim support,
% while declared Sigma_n and the online Sigma_obs proxy are full-space.
% The online agent is a scaled sensor-noise floor, NOT Cov(o) itself.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(fullfile(root, 'lib'));
rng(707,'twister');

p = 8;
ell = 3;
L = 60;
p_minus_ell = p - ell;

% --- Synthetic residual support: o = Q_o * xi, rank(Cov(o)) = p-ell ---
[Qfull, ~] = qr(randn(p, p), 0);
Q_o = Qfull(:, ell+1:end);           % p x (p-ell), orthonormal complement
assert(size(Q_o, 2) == p_minus_ell);
xi = 0.4 * randn(p_minus_ell, L);
o_window = Q_o * xi;                 % exact support in col(Q_o)

% Declared full-space heteroscedastic sensor noise (rank p)
scales = linspace(0.05, 0.25, p)';
Corr = eye(p);
Corr(1,2) = 0.2; Corr(2,1) = 0.2;
D = diag(scales);
Sigma_n = D * Corr * D;
Sigma_n = (Sigma_n + Sigma_n') / 2;

[rank_o, rank_n, trace_o, trace_on, note] = sigma_obs_support_diag( ...
    o_window, Sigma_n, p, ell);

% Residual rank should match unresolved DOF; declared is full rank.
assert(rank_o == p_minus_ell, ...
    sprintf('expected rank_o=%d, got %d', p_minus_ell, rank_o));
assert(rank_n == p, ...
    sprintf('expected rank_n=%d, got %d', p, rank_n));
assert(trace_o > 0 && isfinite(trace_o), 'trace_o invalid');
assert(trace_on > 0 && isfinite(trace_on), 'trace_on invalid');
assert(ischar(note) || isstring(note), 'note must be a string');
assert(contains(lower(char(note)), 'not'), ...
    'note must state proxy is NOT Cov(o)');

% Manual online proxy check (same convention as runner)
O = o_window - mean(o_window, 2);
sigma_hat = norm(O, 'fro') / sqrt(max((p - ell) * (L - 1), 1));
shape_n = Sigma_n / (trace(Sigma_n) / p);
Sigma_on = max(sigma_hat^2, 1e-8) * shape_n;
assert(abs(trace_on - trace(Sigma_on)) < 1e-9 * max(1, abs(trace_on)), ...
    'trace_on mismatch vs manual Sigma_on');

% Empirical Cov(o) rank cannot exceed p-ell, and proxy is full-rank
eo = sort(real(eig((O*O')/max(L-1,1))), 'descend');
assert(sum(eo > 1e-10) <= p_minus_ell, 'empirical Cov(o) rank exceeds p-ell');
assert(rank(Sigma_on) == p, 'online proxy should be full rank p');

% Isotropic declared noise still full-rank proxy
[rank_o2, rank_n2, ~, ~, note2] = sigma_obs_support_diag( ...
    o_window, 0.1^2 * eye(p), p, ell);
assert(rank_o2 == p_minus_ell);
assert(rank_n2 == p);
assert(contains(lower(char(note2)), 'not'));

% Dimension guard
try
    sigma_obs_support_diag(randn(p-1, L), Sigma_n, p, ell); %#ok<NASGU>
    error('test_opinion07_sigma_obs:ExpectedDimError', 'should have failed');
catch ME
    if strcmp(ME.identifier, 'test_opinion07_sigma_obs:ExpectedDimError')
        rethrow(ME);
    end
    % expected dimension error path
end

fprintf(['PASS opinion07 Sigma_obs support: rank_o=%d (=p-ell), ' ...
    'rank_n=%d (=p), trace_o=%.4g, trace_on=%.4g\\n  note: %s\\n'], ...
    rank_o, rank_n, trace_o, trace_on, char(note));
end
