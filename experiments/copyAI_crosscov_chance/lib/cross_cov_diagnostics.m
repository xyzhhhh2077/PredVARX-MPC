function [z, o, Sigma_z, Sigma_o, Sigma_zo, drop_cross_rel_err, extras] = cross_cov_diagnostics(y, P, R, H)
% CROSS_COV_DIAGNOSTICS Oblique latent/residual cross-covariance diagnostics.
%
%   [z,o,Sigma_z,Sigma_o,Sigma_zo,drop_cross_rel_err] = cross_cov_diagnostics(y,P,R)
%   [..., extras] = cross_cov_diagnostics(y,P,R)
%   [..., extras] = cross_cov_diagnostics(y,P,R,H)
%
%   Inputs
%     y : p-by-T measurement samples (raw; mean is removed inside)
%     P : p-by-ell loading
%     R : p-by-ell extractor (dual to P when R'*P = I)
%     H : (optional) nq-by-p chance-constraint directions (rows = hq')
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
%     extras (optional)  : struct with
%                            Sigma_y_emp, Sigma_y_full, Sigma_y_drop
%                            full_recon_rel_err
%                            ||Sigma_zo||_F as Sigma_zo_fro
%                            if H given: var_full, var_drop, var_ratio
%                              (var_* = hq' * Sigma_y_* * hq)
%
%   Algebra (exact when R'*P = I):
%     y_c = P*z + o
%     Cov(y) = P*Sigma_z*P' + P*Sigma_zo + Sigma_zo'*P' + Sigma_o
%            =: Sigma_y_full
%     Sigma_y_drop = P*Sigma_z*P' + Sigma_o   % cross blocks omitted
%
%   Opinion 6: samplewise R'*o = 0 does NOT imply Sigma_zo = 0.
%   Therefore the SMPC construction
%     Sigma_y ~ P*Sigma_z*P' + Sigma_obs
%   silently drops the cross blocks. This diagnostic quantifies that gap.
%
%   Honest boundary: quantifying drop-cross is NOT a re-proof of Boole
%   risk allocation under the full cross-covariance law.

if nargin < 4
    H = [];
end

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

% Drop-cross approximation used by centered_smpc_step default Sigma_y.
Sigma_y_drop = P * Sigma_z * P' + Sigma_o;
Sigma_y_drop = (Sigma_y_drop + Sigma_y_drop') / 2;

% Full block with cross terms.
Sigma_y_full = P * Sigma_z * P' + P * Sigma_zo + Sigma_zo' * P' + Sigma_o;
Sigma_y_full = (Sigma_y_full + Sigma_y_full') / 2;

ny = max(norm(Sigma_y, 'fro'), eps);
drop_cross_rel_err = norm(Sigma_y_drop - Sigma_y, 'fro') / ny;
full_recon_rel_err = norm(Sigma_y_full - Sigma_y, 'fro') / ny;

if nargout >= 7
    extras = struct();
    extras.Sigma_y_emp = Sigma_y;
    extras.Sigma_y_full = Sigma_y_full;
    extras.Sigma_y_drop = Sigma_y_drop;
    extras.full_recon_rel_err = full_recon_rel_err;
    extras.Sigma_zo_fro = norm(Sigma_zo, 'fro');
    extras.drop_cross_rel_err = drop_cross_rel_err;
    if ~isempty(H)
        assert(size(H, 2) == p, 'H must be nq-by-p');
        nq = size(H, 1);
        var_full = zeros(nq, 1);
        var_drop = zeros(nq, 1);
        for q = 1:nq
            hq = H(q, :)';
            var_full(q) = real(hq' * Sigma_y_full * hq);
            var_drop(q) = real(hq' * Sigma_y_drop * hq);
        end
        extras.var_full = var_full;
        extras.var_drop = var_drop;
        extras.var_ratio = var_full ./ max(var_drop, eps);  % full/drop
        extras.H = H;
    end
end
end
