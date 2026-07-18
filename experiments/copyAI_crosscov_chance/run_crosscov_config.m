function report = run_crosscov_config(cfg, plant, data)
% RUN_CROSSCOV_CONFIG Fair closed-loop: cross OFF vs ON (same seed/plant/ID).
%   cfg.use_cross_cov : false (default baseline) | true (inject Sigma_zo)
% Recovery fixed to 'none' so QP success reflects primary chance rows only.

root = fileparts(mfilename('fullpath'));
addpath(root); addpath(fullfile(root, 'lib'));

n = plant.n; m = plant.m; p = plant.p; ell = plant.ell; tracked = plant.tracked;
A = plant.A; B = plant.B; C = plant.C; Sigma_n = plant.Sigma_n; L_n = plant.L_n;
se = plant.se; u_min = plant.u_min; u_max = plant.u_max;
noise_cycle = plant.noise_cycle;
sw_min = plant.sw_min; sw_max = plant.sw_max;
se_min = plant.se_min; se_max = plant.se_max; noise_phase_e = plant.noise_phase_e;

y_off = data.y_off; u_off = data.u_off;
T_cl = cfg.T_cl; N = cfg.N; y_max = cfg.y_max; alpha_joint = cfg.alpha_joint;
use_cross = isfield(cfg, 'use_cross_cov') && logical(cfg.use_cross_cov);

[Ahat, Bhat, Phat, Rhat, Sigma_eps, stats] = split_control_free_ivr_varx( ...
    y_off, u_off, ell, tracked, Sigma_n, 'input_residualize', false);

Hcc = zeros(2, p); Hcc(1, 1) = 1; Hcc(2, 2) = 1;
[~, ~, Sz_off, So_off, Szo, drop_rel, xtra] = cross_cov_diagnostics(y_off, Phat, Rhat, Hcc);
blocks = sigma_y_chance_blocks(Phat, Sz_off, Szo, So_off, Hcc);

model.A = Ahat; model.B = Bhat; model.P = Phat; model.R = Rhat;
model.y_mean = stats.y_mean; model.u_mean = stats.u_mean;
model.Sigma_eps = (Sigma_eps + Sigma_eps') / 2;
model.Sigma_obs = Sigma_n;
model.Sigma_zo = Szo;   % always stored; used only if opt.use_cross_cov

opt.N = N; opt.Q = zeros(p); opt.Q(1, 1) = 80; opt.Q(2, 2) = 80; opt.Ru = 0.25 * eye(m);
opt.H = Hcc; opt.h = y_max * [1; 1];
opt.u_min = u_min; opt.u_max = u_max; opt.alpha_joint = alpha_joint;
opt.use_cross_cov = use_cross;
opt.use_terminal_cost = false;

seg = T_cl / 5; Rf = zeros(p, T_cl);
rset = [0.6 -0.4; 1.2 0.3; -0.5 1.0; 0.9 -0.8; 0.2 0.5];
for s = 1:5
    idx = (s - 1) * seg + 1:s * seg;
    Rf(1, idx) = rset(s, 1); Rf(2, idx) = rset(s, 2);
end
[swp, sep] = smooth_noise_profile(T_cl, noise_cycle, sw_min, sw_max, se_min, se_max, noise_phase_e);
rng(cfg.cl_seed, 'twister');

x = zeros(n, 1); y = zeros(p, T_cl); u = zeros(m, T_cl);
level = repmat({''}, 1, T_cl);
nw = 40; eb = zeros(ell, nw); ec = 0; ob = zeros(p, nw); oc = 0;
Ipr = eye(p) - Phat * Rhat'; zprev = [];  % I - P R'
n_qp = 0; n_fail = 0;
max_cc_viol = nan(1, T_cl);
active_flags = false(1, T_cl);
cross_flag_hist = false(1, T_cl);

for k = 1:T_cl
    vk = (sep(k) / se) * L_n * randn(p, 1);
    yk = C * x + vk; y(:, k) = yk;
    zk = model.R' * (yk - model.y_mean);
    ores = Ipr * (yk - model.y_mean);
    io = mod(k - 1, nw) + 1; ob(:, io) = ores; oc = min(oc + 1, nw);
    if k >= 2
        er = zk - model.A * zprev - model.B * (u(:, k - 1) - model.u_mean);
        ie = mod(k - 2, nw) + 1; eb(:, ie) = er; ec = min(ec + 1, nw);
    end
    if ec >= 5
        E = eb(:, 1:ec) - mean(eb(:, 1:ec), 2); G = E * E'; Nr = ec;
        model.Sigma_eps = (G / max(Nr - 1, 1)) + 1e-8 * eye(ell);
        model.Sigma_eps = (model.Sigma_eps + model.Sigma_eps') / 2;
    end
    if oc >= 5
        model.Sigma_obs = build_sigma_obs_trial('declared_shape', ob(:, 1:oc), Sigma_n, p, ell, 1e-8);
    end
    rk = Rf(:, min(k + 1, T_cl));
    try
        [~, ~, U, info] = centered_smpc_step(yk, rk, model, opt);
        uk = U(1:m);
        mv = max(info.A_ch * U - info.b_ch);
        max_cc_viol(k) = mv;
        active_flags(k) = (mv > -1e-3);
        cross_flag_hist(k) = isfield(info, 'use_cross_cov') && logical(info.use_cross_cov);
        if ~(info.exitflag > 0 && mv <= 1e-7)
            error('primary row fail');
        end
        level{k} = 'qp_primary';
        n_qp = n_qp + 1;
    catch
        n_fail = n_fail + 1;
        uk = min(max(model.u_mean, u_min), u_max);
        level{k} = 'uncertified_fallback';
        max_cc_viol(k) = NaN;
        active_flags(k) = false;
        cross_flag_hist(k) = use_cross;
    end
    u(:, k) = uk;
    x = A * x + B * uk + swp(k) * randn(n, 1);
    zprev = zk;
end

warm = max(1, floor(0.2 * T_cl)):T_cl;
err = y(tracked, warm) - Rf(tracked, warm);
report = struct();
report.name = cfg.name;
report.seed = cfg.cl_seed;
report.cfg = cfg;
report.use_cross_cov = use_cross;
report.MAE = mean(abs(err), 2);
report.RMSE = sqrt(mean(err .^ 2, 2));
report.joint_viol = mean(any(y(tracked, :) > y_max, 1));
report.joint_cover = 1 - report.joint_viol;
report.qp_success = n_qp / T_cl;
report.qp_fail = n_fail;
report.n_qp = n_qp;
report.uncert_rate = n_fail / T_cl;
report.active_rate = mean(active_flags(warm), 'omitnan');
report.max_cc_viol_warm = max(max_cc_viol(warm), [], 'omitnan');
report.cross_flag_rate = mean(cross_flag_hist);

% Offline diagnostics (identical for OFF/ON pair with same plant/data)
report.Sigma_zo_fro = norm(Szo, 'fro');
report.drop_cross_rel_err = drop_rel;
report.full_recon_rel_err = xtra.full_recon_rel_err;
report.var_full = blocks.var_full(:)';
report.var_drop = blocks.var_drop(:)';
report.var_ratio = blocks.var_ratio(:)';  % full/drop on H rows
report.cross_block_fro = blocks.cross_block_fro;
report.left_err = stats.tracked_left_error;
report.dual_err = stats.dual_error;
end
