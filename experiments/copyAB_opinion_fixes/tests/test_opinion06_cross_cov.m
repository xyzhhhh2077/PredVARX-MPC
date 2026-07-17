function test_opinion06_cross_cov
% TEST_OPINION06_CROSS_COV Opinion 6: R'*o=0 does not imply Cov(z,o)=0.
% Random data: full reconstruction relative error ~0, drop_cross can be >0.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(fullfile(root, 'lib'));
rng(20260718, 'twister');

p = 8; ell = 3; T = 800;

% --- Case A: random oblique dual pair + generic Gaussian data ------------
[P, R] = local_random_dual_pair(p, ell);
y = randn(p, T);
y = y + 0.4 * [1; zeros(p-1,1)] * randn(1, T);  % mild cross-channel structure

[z, o, Sz, So, Szo, drop_rel] = cross_cov_diagnostics(y, P, R);

assert(isequal(size(z), [ell T]), 'z size mismatch');
assert(isequal(size(o), [p T]), 'o size mismatch');
assert(isequal(size(Sz), [ell ell]), 'Sigma_z size mismatch');
assert(isequal(size(So), [p p]), 'Sigma_o size mismatch');
assert(isequal(size(Szo), [ell p]), 'Sigma_zo size mismatch');

% Samplewise dual residual identity (R'*o = 0 when R'*P = I)
assert(norm(R' * o, 'fro') < 1e-10, 'R''*o should be ~0 samplewise');
assert(norm(R' * P - eye(ell), 'fro') < 1e-10, 'dual pair broken');

% Exact reconstruction of each sample
yc = y - mean(y, 2);
assert(norm(yc - (P * z + o), 'fro') < 1e-10, 'y_c ~= Pz + o');

% Full block reconstruction of Cov(y) should match empirical Cov(y)
den = max(T - 1, 1);
Sy = (yc * yc') / den;
Sy = (Sy + Sy') / 2;
Sy_full = P * Sz * P' + P * Szo + Szo' * P' + So;
Sy_full = (Sy_full + Sy_full') / 2;
full_rel = norm(Sy_full - Sy, 'fro') / max(norm(Sy, 'fro'), eps);
assert(full_rel < 1e-10, sprintf('full recon rel err too large: %.3e', full_rel));

% Drop-cross relative error must be reported and can be strictly positive
assert(isfinite(drop_rel) && drop_rel >= 0, 'drop_cross_rel_err invalid');
assert(drop_rel > 1e-4, ...
    sprintf('expected drop_cross_rel_err > 0 on generic data, got %.3e', drop_rel));
assert(norm(Szo, 'fro') > 1e-6, 'Sigma_zo should be nonzero on generic data');

% --- Case B: orthogonal extractor on isotropic data tends to small drop --
% (not required zero: P principal subspace of identity is any ONB, but
%  random orthogonal P with isotropic y has E[z o']=0 asymptotically)
[Q,~] = qr(randn(p, ell), 0);
Po = Q; Ro = Q;
y_iso = randn(p, T);
[~, ~, ~, ~, Szo_iso, drop_iso] = cross_cov_diagnostics(y_iso, Po, Ro);
% With isotropic samples and orthogonal dual, cross term is O(1/sqrt(T))
assert(drop_iso < 0.15, ...
    sprintf('isotropic orthogonal drop should be modest, got %.3e', drop_iso));
assert(norm(Szo_iso, 'fro') < 0.5, 'isotropic Sigma_zo unexpectedly large');

% --- Case C: oblique free dual as in copyAB (noise-optimal) --------------
% Reproduce a small split geometry: tracked E + free oblique dual.
tracked = [1 2]; q = numel(tracked);
E = zeros(p, q); E(tracked, :) = eye(q);
Nperp = null(E');
r = ell - q;
V = orth(randn(size(Nperp, 2), r));
Sigma_n = diag(linspace(0.04, 0.55, p).^2);
Sigma_n(3,4) = 0.02; Sigma_n(4,3) = 0.02;
Sigma_n = (Sigma_n + Sigma_n') / 2;
Sperp = Nperp' * Sigma_n * Nperp;
SinvV = Sperp \ V;
W = SinvV / (V' * SinvV);
Psplit = [E, Nperp * V];
Rsplit = [E, Nperp * W];
y_het = chol(Sigma_n, 'lower') * randn(p, T) + 0.7 * Psplit * randn(ell, T);
[~, ~, ~, ~, Szo_s, drop_s] = cross_cov_diagnostics(y_het, Psplit, Rsplit);
assert(norm(Rsplit' * Psplit - eye(ell), 'fro') < 1e-10, 'split dual failed');
assert(drop_s > 1e-4, sprintf('split oblique drop should be >0, got %.3e', drop_s));

fprintf(['PASS opinion06 cross-cov: full_rel=%.3e, drop_generic=%.3e, ' ...
    '||Szo||_F=%.3e; iso drop=%.3e; split drop=%.3e, ||Szo_split||_F=%.3e\n'], ...
    full_rel, drop_rel, norm(Szo, 'fro'), drop_iso, drop_s, norm(Szo_s, 'fro'));
end

function [P, R] = local_random_dual_pair(p, ell)
% Random full-rank P and a dual R with R'*P = I (not necessarily R=P).
P = randn(p, ell);
% Random positive definite metric for a weighted dual (oblique)
M = randn(p); M = M * M' + 0.3 * eye(p);
R = M \ P / (P' * (M \ P));
% Optional: slightly mix R away from pure metric dual while keeping dual
% identity — already dual; done.
end
