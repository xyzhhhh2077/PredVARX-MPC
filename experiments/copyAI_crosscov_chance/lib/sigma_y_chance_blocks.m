function out = sigma_y_chance_blocks(P, Sigma_z, Sigma_zo, Sigma_o_or_obs, H)
% SIGMA_Y_CHANCE_BLOCKS Build drop-cross vs full Sigma_y and hq variances.
%
%   out = sigma_y_chance_blocks(P, Sigma_z, Sigma_zo, Sigma_o, H)
%
%   Sigma_y_full = P*Sz*P' + P*Szo + Szo'*P' + So
%   Sigma_y_drop = P*Sz*P' + So
%
%   For each row hq' of H (nq-by-p):
%     var_full(q) = hq' * Sigma_y_full * hq
%     var_drop(q) = hq' * Sigma_y_drop * hq
%     var_ratio(q) = var_full / max(var_drop, eps)
%
%   So may be empirical Sigma_o or the online Sigma_obs proxy.
%   This helper does NOT claim Boole risk allocation under full cross terms.

ny = size(P, 1);
nz = size(P, 2);
assert(isequal(size(Sigma_z), [nz nz]), 'Sigma_z size');
assert(isequal(size(Sigma_o_or_obs), [ny ny]), 'Sigma_o size');
if isempty(Sigma_zo)
    Sigma_zo = zeros(nz, ny);
end
assert(isequal(size(Sigma_zo), [nz ny]) || isequal(size(Sigma_zo), [ny nz]), ...
    'Sigma_zo must be ell-by-p or p-by-ell');
if size(Sigma_zo, 1) == ny && size(Sigma_zo, 2) == nz
    Sigma_zo = Sigma_zo';
end

So = (Sigma_o_or_obs + Sigma_o_or_obs') / 2;
Sz = (Sigma_z + Sigma_z') / 2;
Szo = Sigma_zo;

Sigma_y_drop = P * Sz * P' + So;
Sigma_y_drop = (Sigma_y_drop + Sigma_y_drop') / 2;
Sigma_y_full = Sigma_y_drop + P * Szo + Szo' * P';
Sigma_y_full = (Sigma_y_full + Sigma_y_full') / 2;

out = struct();
out.Sigma_y_drop = Sigma_y_drop;
out.Sigma_y_full = Sigma_y_full;
out.Sigma_zo = Szo;
out.Sigma_zo_fro = norm(Szo, 'fro');
out.cross_block_fro = norm(P * Szo + Szo' * P', 'fro');
out.drop_vs_full_rel = norm(Sigma_y_drop - Sigma_y_full, 'fro') / ...
    max(norm(Sigma_y_full, 'fro'), eps);

if nargin >= 5 && ~isempty(H)
    assert(size(H, 2) == ny, 'H must be nq-by-p');
    nq = size(H, 1);
    var_full = zeros(nq, 1);
    var_drop = zeros(nq, 1);
    for q = 1:nq
        hq = H(q, :)';
        var_full(q) = max(real(hq' * Sigma_y_full * hq), 0);
        var_drop(q) = max(real(hq' * Sigma_y_drop * hq), 0);
    end
    out.var_full = var_full;
    out.var_drop = var_drop;
    out.var_ratio = var_full ./ max(var_drop, eps);
    out.H = H;
end
end
