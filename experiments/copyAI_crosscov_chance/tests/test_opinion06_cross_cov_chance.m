function test_opinion06_cross_cov_chance
% TEST_OPINION06_CROSS_COV_CHANCE
% R'*o ~ 0 does not force Sigma_zo = 0; full recon err ~ 0;
% centered_smpc_step default use_cross_cov=false; ON path uses full blocks.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(fullfile(root, 'lib'));
rng(20260719, 'twister');

p = 8; ell = 3; T = 900;

% --- Case A: oblique dual + generic data ---------------------------------
[P, R] = local_random_dual_pair(p, ell);
y = randn(p, T);
y = y + 0.5 * [1; 0.4; zeros(p - 2, 1)] * randn(1, T);

H = zeros(2, p); H(1, 1) = 1; H(2, 2) = 1;
[z, o, Sz, So, Szo, drop_rel, xtra] = cross_cov_diagnostics(y, P, R, H);

assert(norm(R' * P - eye(ell), 'fro') < 1e-10, 'dual broken');
assert(norm(R' * o, 'fro') < 1e-10, 'R''*o should be ~0');
assert(xtra.full_recon_rel_err < 1e-10, ...
    sprintf('full recon too large: %.3e', xtra.full_recon_rel_err));
assert(drop_rel > 1e-4, sprintf('drop should be >0, got %.3e', drop_rel));
assert(norm(Szo, 'fro') > 1e-6, 'Sigma_zo should be nonzero');
assert(all(isfinite(xtra.var_ratio)) && numel(xtra.var_ratio) == 2, 'var_ratio');

blocks = sigma_y_chance_blocks(P, Sz, Szo, So, H);
assert(norm(blocks.Sigma_y_full - xtra.Sigma_y_full, 'fro') < 1e-12, 'block full');
assert(norm(blocks.Sigma_y_drop - xtra.Sigma_y_drop, 'fro') < 1e-12, 'block drop');
assert(max(abs(blocks.var_ratio(:) - xtra.var_ratio(:))) < 1e-12, 'ratio match');

% --- Case B: orthogonal isotropic -> modest drop -------------------------
[Q, ~] = qr(randn(p, ell), 0);
[~, ~, ~, ~, Szo_iso, drop_iso] = cross_cov_diagnostics(randn(p, T), Q, Q);
assert(drop_iso < 0.15, sprintf('iso drop large: %.3e', drop_iso));
assert(norm(Szo_iso, 'fro') < 0.5, 'iso Szo large');

% --- Case C: centered_smpc_step flag default OFF, ON uses cross ----------
ny = p; nu = 2; nz = ell; N = 3;
model.A = 0.85 * eye(nz);
model.B = 0.1 * ones(nz, nu);
model.P = P; model.R = R;
model.y_mean = zeros(ny, 1); model.u_mean = zeros(nu, 1);
model.Sigma_eps = 1e-3 * eye(nz);
model.Sigma_obs = 1e-3 * eye(ny);
model.Sigma_zo = Szo;

opt.N = N; opt.Q = eye(ny); opt.Ru = 0.1 * eye(nu);
opt.H = H; opt.h = 5 * [1; 1];
opt.u_min = -10; opt.u_max = 10; opt.alpha_joint = 0.2;
opt.use_terminal_cost = false;

% default / explicit false
opt0 = opt; opt0.use_cross_cov = false;
[~, ~, U0, out0] = centered_smpc_step(zeros(ny, 1), zeros(ny, 1), model, opt0);
assert(out0.use_cross_cov == false, 'default must stay drop-cross');

opt1 = opt; opt1.use_cross_cov = true;
[~, ~, U1, out1] = centered_smpc_step(zeros(ny, 1), zeros(ny, 1), model, opt1);
assert(out1.use_cross_cov == true, 'ON must report use_cross_cov');

% With nonzero Szo, tightening can differ => b_ch may differ
if norm(Szo, 'fro') > 1e-4
    assert(norm(out0.b_ch - out1.b_ch) > 1e-10 || norm(U0 - U1) > 1e-10, ...
        'expected OFF/ON chance rows or U to differ when Szo nonzero');
end

% Missing Sigma_zo with flag true should fall back (use_cross false)
model2 = model; model2 = rmfield(model2, 'Sigma_zo');
[~, ~, ~, out2] = centered_smpc_step(zeros(ny, 1), zeros(ny, 1), model2, opt1);
assert(out2.use_cross_cov == false, 'missing Sigma_zo must disable cross');

fprintf(['PASS opinion06 chance: full_rel=%.3e drop=%.3e ||Szo||=%.3e ' ...
    'vr=[%.3f %.3f] OFF/ON b_diff=%.3e\n'], ...
    xtra.full_recon_rel_err, drop_rel, norm(Szo, 'fro'), ...
    xtra.var_ratio(1), xtra.var_ratio(2), norm(out0.b_ch - out1.b_ch));
end

function [P, R] = local_random_dual_pair(p, ell)
P = randn(p, ell);
M = randn(p); M = M * M' + 0.3 * eye(p);
R = M \ P / (P' * (M \ P));
end
