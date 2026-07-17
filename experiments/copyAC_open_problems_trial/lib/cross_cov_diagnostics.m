function [z, o, Sigma_z, Sigma_o, Sigma_zo, drop_cross_rel_err] = cross_cov_diagnostics(y, P, R)
% CROSS_COV_DIAGNOSTICS Oblique latent/residual cross-covariance diagnostics.
%
%   [z,o,Sigma_z,Sigma_o,Sigma_zo,drop_cross_rel_err] = cross_cov_diagnostics(y,P,R)
%
%   Inputs
%     y : p-by-T measurement samples (raw; mean is removed inside)
%     P : p-by-ell loading
%     R : p-by-ell extractor (dual to P when R'*P = I)
%
%   Outputs
%     z                  : ell-by-T latent coordinates  z = R' * y_c
%     o                  : p-by-T residual o = (I - P*R') * y_c
%     Sigma_z            : ell-by-ell sample Cov(z)
%     Sigma_o            : p-by-p   sample Cov(o)
%     Sigma_zo           : ell-by-p  sample Cov(z,o) = E[z o']
%     drop_cross_rel_err : relative Frobenius error of the drop-cross
%                          SMPC-style approximation of Cov(y):
%                            Sigma_y_drop = P*Sigma_z*P' + Sigma_o
%                          versus the empirical Cov(y).
%
%   Algebra (exact when R'*P = I):
%     y_c = P*z + o
%     Cov(y) = P*Sigma_z*P' + P*Sigma_zo + Sigma_zo'*P' + Sigma_o
%
%   Opinion 6: samplewise R'*o = 0 does NOT imply Sigma_zo = 0.
%   Therefore the SMPC construction
%     Sigma_y ≈ P*Sigma_z*P' + Sigma_obs
%   silently drops the cross blocks. This diagnostic quantifies that gap.
%
%   Full block reconstruction relative error is ~0 by construction
%   (up to floating-point). drop_cross_rel_err can be strictly positive
%   for generic data / oblique extractors.

[p, T] = size(y);
assert(isequal(size(P,1), p) && isequal(size(R,1), p), ...
    'cross_cov_diagnostics: P and R must have p rows matching y.');
assert(isequal(size(P,2), size(R,2)), ...
    'cross_cov_diagnostics: P and R must share latent dimension ell.');
assert(T >= 2, 'cross_cov_diagnostics: need at least two samples.');

yc = y - mean(y, 2);
z = R' * yc;
o = yc - P * z;   % == (I - P*R')*yc

den = max(T - 1, 1);
Sigma_z  = (z * z') / den;
Sigma_z  = (Sigma_z + Sigma_z') / 2;
Sigma_o  = (o * o') / den;
Sigma_o  = (Sigma_o + Sigma_o') / 2;
Sigma_zo = (z * o') / den;          % ell-by-p; not necessarily zero
Sigma_y  = (yc * yc') / den;
Sigma_y  = (Sigma_y + Sigma_y') / 2;

% Drop-cross approximation used by centered_smpc_step-style Sigma_y.
Sigma_y_drop = P * Sigma_z * P' + Sigma_o;
Sigma_y_drop = (Sigma_y_drop + Sigma_y_drop') / 2;

ny = max(norm(Sigma_y, 'fro'), eps);
drop_cross_rel_err = norm(Sigma_y_drop - Sigma_y, 'fro') / ny;
end
