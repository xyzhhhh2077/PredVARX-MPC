%% copyAC_open_problems_trial
% Trial solutions for open pure-math opinions 5/6/7/8/9/10.
% Does NOT modify main/, copyX, copyAA, or copyAB.
%
% Hard boundaries (still true after this trial):
% - Soft recovery is NOT the original chance certificate.
% - Cross-term Sigma_y is an experimental inclusion, not a proved Boole law.
% - residual_support / additive Sigma_obs are trial objects, not Cov(o) theorems.
% - Terminal cost default ON here is a soft regularizer only (not stability proof).
% - Input residualize IVR is a candidate, not proved input-conditional PredVARX.

clear; clc; close all;
root_dir = fileparts(mfilename('fullpath'));
addpath(root_dir);
addpath(fullfile(root_dir, 'lib'));
rng(20260710,'twister');

%% ===== OPEN-PROBLEM TRIAL SWITCHES =====
trial = struct();
trial.sigma_eps_mode = 'ols';            % 't2' | 'ml' | 'ols'   (opinion 5)
trial.use_cross_cov = true;             % opinion 6: inject Sigma_zo into Sigma_y
trial.sigma_obs_mode = 'additive';       % 'declared_shape'|'residual_support'|'additive' (op7)
trial.enable_soft_recovery = true;      % opinion 8 L2
trial.use_terminal_cost = true;         % opinion 9 optional (soft only)
trial.input_residualize = true;         % opinion 10 candidate IVR
trial.force_primary_infeasible_once = false; % unit-test style: not used in full run

%% Process-like plant
n = 6; m = 3; p = 30; ell = 5; tracked = [1 2];
T_off = 1500; T_cl = 1200; N = 18;
sw = 0.045; se = 0.055;
noise_cycle = 400;
sw_min = 0.020; sw_max = 0.090;
se_min = 0.025; se_max = 0.100;
noise_phase_e = pi/3;
u_min = -3.0; u_max = 3.0; y_max = 2.00; alpha_joint = 0.10;

A = diag([0.94 0.88 0.78 0.64 0.50 0.35]);
A(1,2) = 0.10; A(2,3) = -0.06; A(3,4) = 0.05; A(4,5) = 0.04;
B = [0.34 -0.10  0.05;
     0.12  0.28 -0.06;
     0.05  0.12  0.24;
    -0.05  0.06  0.18;
     0.02 -0.10  0.14;
     0.08  0.02 -0.08];
C = zeros(p,n);
C(1,1) = 1.00; C(1,3) = 0.16;
C(2,2) = 1.00; C(2,4) = -0.12;
for i = 3:p
    C(i,:) = 0.45*randn(1,n);
    C(i,:) = C(i,:) / max(norm(C(i,:)), 1e-12);
end

sensor_rel = linspace(0.55,1.65,p)';
Corr_n = eye(p);
Corr_n(3,4) = 0.30; Corr_n(4,3) = 0.30;
Corr_n(5,6) = -0.22; Corr_n(6,5) = -0.22;
Corr_n(8,9) = 0.18; Corr_n(9,8) = 0.18;
D_n = diag(se*sensor_rel);
Sigma_n = D_n*Corr_n*D_n;
L_n = chol(Sigma_n,'lower');

%% Offline data
u_off = 1.20*randn(m,T_off);
x_off = zeros(n,T_off+1);
y_off = zeros(p,T_off);
for k = 1:T_off
    y_off(:,k) = C*x_off(:,k) + L_n*randn(p,1);
    x_off(:,k+1) = A*x_off(:,k) + B*u_off(:,k) + sw*randn(n,1);
end

%% Identify (opinion 10 optional residualize)
[Ahat,Bhat,Phat,Rhat,Sigma_eps,stats] = split_control_free_ivr_varx( ...
    y_off,u_off,ell,tracked,Sigma_n, ...
    'input_residualize', trial.input_residualize);
stats.true_A_eigs = eig(A);
stats.id_A_eigs = eig(Ahat);
stats.ell = ell;
stats.tracked = tracked;
stats.trial = trial;

% Opinion 5: choose Sigma_eps denom for controller
switch lower(trial.sigma_eps_mode)
    case 'ml'
        Sigma_eps_use = stats.Sigma_eps_ml;
    case 'ols'
        Sigma_eps_use = stats.Sigma_eps_ols;
    otherwise
        Sigma_eps_use = Sigma_eps; % t2
        trial.sigma_eps_mode = 't2';
end
stats.Sigma_eps_used_mode = trial.sigma_eps_mode;
stats.Sigma_eps_used = Sigma_eps_use;

% Opinion 6 offline diag + fixed offline Sigma_zo for trial injection
[~, o_xcov, Sz_xcov, So_xcov, Szo_xcov, drop_cross_rel_err] = ...
    cross_cov_diagnostics(y_off, Phat, Rhat);
yc_off = y_off - mean(y_off, 2);
den_off = max(size(y_off, 2) - 1, 1);
Sy_off = (yc_off * yc_off') / den_off; Sy_off = (Sy_off+Sy_off')/2;
Sy_full_off = Phat*Sz_xcov*Phat' + Phat*Szo_xcov + Szo_xcov'*Phat' + So_xcov;
Sy_full_off = (Sy_full_off+Sy_full_off')/2;
full_cross_rel_err = norm(Sy_full_off-Sy_off,'fro')/max(norm(Sy_off,'fro'),eps);
stats.cross_cov_full_rel_err = full_cross_rel_err;
stats.cross_cov_drop_rel_err = drop_cross_rel_err;
stats.cross_cov_Sigma_zo_fro = norm(Szo_xcov,'fro');
fprintf(['Opinion6 cross-cov offline: full_rel=%.3e drop=%.3e ||Szo||=%.3e ' ...
    'max|R''o|=%.3e use_cross=%d\n'], full_cross_rel_err, drop_cross_rel_err, ...
    norm(Szo_xcov,'fro'), max(abs(Rhat'*o_xcov),[],'all'), trial.use_cross_cov);

model.A = Ahat; model.B = Bhat; model.P = Phat; model.R = Rhat;
model.y_mean = stats.y_mean; model.u_mean = stats.u_mean;
model.Sigma_eps = (Sigma_eps_use+Sigma_eps_use')/2;
model.Sigma_obs = Sigma_n;
model.Sigma_zo = Szo_xcov;  % ell x p from cross_cov_diagnostics
Sigma_obs_shape = Sigma_n / (trace(Sigma_n)/p);

opt.N = N;
opt.Q = zeros(p); opt.Q(1,1) = 80; opt.Q(2,2) = 80;
opt.Ru = 0.25*eye(m);
opt.H = [1 0; 0 1]; % will be rebuilt as selection of tracked rows via eye
% Use tracked quality bounds: H selects y1,y2
opt.H = zeros(2,p); opt.H(1,1)=1; opt.H(2,2)=1;
opt.h = y_max*[1;1];
opt.u_min = u_min; opt.u_max = u_max;
opt.alpha_joint = alpha_joint;
opt.use_cross_cov = trial.use_cross_cov;
opt.use_terminal_cost = trial.use_terminal_cost;

% References
seg_len = T_cl/5;
Rf = zeros(p,T_cl);
for s = 1:5
    idx = (s-1)*seg_len+1:s*seg_len;
    switch s
        case 1, r1=0.6; r2=-0.4;
        case 2, r1=1.2; r2=0.3;
        case 3, r1=-0.5; r2=1.0;
        case 4, r1=0.9; r2=-0.8;
        case 5, r1=0.2; r2=0.5;
    end
    Rf(1,idx)=r1; Rf(2,idx)=r2;
end
[sigma_w_profile, sigma_e_profile] = smooth_noise_profile( ...
    T_cl, noise_cycle, sw_min, sw_max, se_min, se_max, noise_phase_e);

%% Closed loop
x = zeros(n,1);
y = zeros(p,T_cl); yhat = zeros(p,T_cl); u = zeros(m,T_cl);
exitflag = nan(1,T_cl); max_cc_violation = nan(1,T_cl); costJ = nan(1,T_cl);
estimated_sigma_eps = nan(1,T_cl); estimated_sigma_obs = nan(1,T_cl);
true_sigma_w = nan(1,T_cl); true_sigma_e = nan(1,T_cl);
cc_cert_level = repmat({''},1,T_cl);
soft_mode = repmat({''},1,T_cl);
fallback_det_mean_violation = nan(1,T_cl);
noise_window = 40;
eps_buffer = zeros(ell,noise_window); eps_count = 0;
obs_buffer = zeros(p,noise_window); obs_count = 0;
IminusPR = eye(p)-Phat*Rhat';
infeasible_count = 0; soft_recovery_count = 0; uncertified_count = 0;
last_sigma_obs_diag = struct();
z_prev = [];

for k = 1:T_cl
    vk = (sigma_e_profile(k)/se)*L_n*randn(p,1);
    yk = C*x + vk;
    true_sigma_e(k) = norm(vk)/sqrt(p);
    y(:,k) = yk;
    zk_now = model.R'*(yk-model.y_mean);
    obs_res = IminusPR*(yk-model.y_mean);
    io = mod(k-1,noise_window)+1;
    obs_buffer(:,io) = obs_res;
    obs_count = min(obs_count+1,noise_window);
    if k >= 2
        eps_res = zk_now-model.A*z_prev-model.B*(u(:,k-1)-model.u_mean);
        ie = mod(k-2,noise_window)+1;
        eps_buffer(:,ie) = eps_res;
        eps_count = min(eps_count+1,noise_window);
    end
    if eps_count >= 5
        Ewin = eps_buffer(:,1:eps_count)-mean(eps_buffer(:,1:eps_count),2);
        Geps = Ewin*Ewin';
        Nres = eps_count;
        switch lower(trial.sigma_eps_mode)
            case 'ml', den = max(Nres,1);
            case 'ols', den = max(Nres-(ell+m),1);
            otherwise, den = max(Nres-1,1);
        end
        model.Sigma_eps = (Geps/den + 1e-8*eye(ell));
        model.Sigma_eps = (model.Sigma_eps+model.Sigma_eps')/2;
    end
    if obs_count >= 5
        Owin = obs_buffer(:,1:obs_count);
        [model.Sigma_obs, meta7] = build_sigma_obs_trial( ...
            trial.sigma_obs_mode, Owin, Sigma_n, p, ell, 1e-8);
        last_sigma_obs_diag.k = k;
        last_sigma_obs_diag.meta = meta7;
        try
            [r_o,r_n,tr_o,tr_on,note7] = sigma_obs_support_diag(Owin, Sigma_n, p, ell);
            last_sigma_obs_diag.rank_o=r_o; last_sigma_obs_diag.rank_n=r_n;
            last_sigma_obs_diag.trace_o=tr_o; last_sigma_obs_diag.trace_on=tr_on;
            last_sigma_obs_diag.note=note7;
        catch
        end
    end
    estimated_sigma_eps(k) = sqrt(trace(model.Sigma_eps)/ell);
    estimated_sigma_obs(k) = sqrt(trace(model.Sigma_obs)/p);
    rk = Rf(:,min(k+1,T_cl));

    try
        [~, ypred, U, info] = centered_smpc_step(yk, rk, model, opt);
        uk = U(1:m);
        yhat(:,k) = ypred;
        exitflag(k) = info.exitflag;
        max_cc_violation(k) = max(info.A_ch*U - info.b_ch);
        costJ(k) = info.cost;
        cc_cert_level{k} = 'qp';
        soft_mode{k} = 'primary';
    catch ME
        infeasible_count = infeasible_count + 1;
        recovered = false;
        if trial.enable_soft_recovery
            fb2 = soft_recovery_smpc(yk, rk, model, opt);
            if fb2.success
                recovered = true;
                soft_recovery_count = soft_recovery_count + 1;
                uk = fb2.uk;
                yhat(:,k) = fb2.y_pred;
                exitflag(k) = fb2.exitflag;
                max_cc_violation(k) = fb2.max_cc_violation;
                costJ(k) = fb2.cost;
                cc_cert_level{k} = 'soft_recovery';
                soft_mode{k} = fb2.mode;
                if soft_recovery_count <= 3
                    fprintf('soft_recovery at k=%d: %s | primary: %s\n', ...
                        k, fb2.note, ME.message);
                end
            end
        end
        if ~recovered
            uncertified_count = uncertified_count + 1;
            uk = min(max(model.u_mean, u_min), u_max);
            fb = fallback_certify_step(yk, model, opt, uk);
            exitflag(k) = fb.exitflag;
            max_cc_violation(k) = fb.max_cc_violation;
            cc_cert_level{k} = 'uncertified_fallback';
            soft_mode{k} = 'none';
            fallback_det_mean_violation(k) = fb.det_mean_violation;
            if ~isempty(fb.mean_pred)
                yhat(:,k) = fb.mean_pred;
            elseif k > 1
                yhat(:,k) = yhat(:,k-1);
            else
                yhat(:,k) = model.y_mean;
            end
            if uncertified_count <= 3
                fprintf('uncertified_fallback at k=%d: %s | %s\n', k, ME.message, fb.note);
            end
        end
    end
    u(:,k) = uk;
    wk = sigma_w_profile(k)*randn(n,1);
    true_sigma_w(k) = norm(wk)/sqrt(n);
    x = A*x + B*uk + wk;
    z_prev = zk_now;
end

%% Metrics
warm = 151:T_cl;
err = y(tracked,warm) - Rf(tracked,warm);
MAE = mean(abs(err),2);
RMSE = sqrt(mean(err.^2,2));
Bias = mean(err,2);
upper_violation_rate = sum(y(tracked,:) > y_max, 2) / T_cl;
abs_violation_rate = sum(abs(y(tracked,:)) > y_max, 2) / T_cl;
u_rms = sqrt(mean(u.^2,2));
u_sat_rate = mean(abs(u-u_min)<1e-8 | abs(u-u_max)<1e-8, 2);

qp_certified_count = sum(strcmp(cc_cert_level,'qp'));
soft_recovery_count_chk = sum(strcmp(cc_cert_level,'soft_recovery'));
uncertified_fallback_count = sum(strcmp(cc_cert_level,'uncertified_fallback'));
qp_certified_rate = qp_certified_count / T_cl;
soft_recovery_rate = soft_recovery_count_chk / T_cl;
uncertified_fallback_rate = uncertified_fallback_count / T_cl;
qp_mask = strcmp(cc_cert_level,'qp');
max_recorded_qp_constraint = max(max_cc_violation(qp_mask & ~isnan(max_cc_violation)));
if isempty(max_recorded_qp_constraint), max_recorded_qp_constraint = NaN; end
warm_qp = warm(strcmp(cc_cert_level(warm),'qp'));
if isempty(warm_qp)
    constraint_active_rate = NaN;
else
    constraint_active_rate = mean(max_cc_violation(warm_qp) > -1e-3, 'omitnan');
end

fprintf('\ncopyAC open-problems trial results\n');
fprintf('trial: eps=%s cross=%d obs=%s soft=%d terminal=%d residualize=%d\n', ...
    trial.sigma_eps_mode, trial.use_cross_cov, trial.sigma_obs_mode, ...
    trial.enable_soft_recovery, trial.use_terminal_cost, trial.input_residualize);
fprintf('geometry left/right/col/dual = %.3e / %.3e / %.3e / %.3e\n', ...
    stats.tracked_left_error, stats.tracked_right_error, ...
    stats.tracked_column_error, stats.dual_error);
fprintf('ivr_input_conditional=%d free_oblique=%.3e noise_gain=%.3e\n', ...
    stats.ivr_input_conditional, stats.free_oblique_norm, stats.free_noise_improvement);
fprintf('MAE=[%.4f %.4f] RMSE=[%.4f %.4f] Bias=[%.4f %.4f]\n', ...
    MAE(1),MAE(2),RMSE(1),RMSE(2),Bias(1),Bias(2));
fprintf('upper_viol=[%.4f %.4f] abs_viol=[%.4f %.4f]\n', ...
    upper_violation_rate(1),upper_violation_rate(2), abs_violation_rate(1),abs_violation_rate(2));
fprintf('cert: qp=%.4f (%d) soft=%.4f (%d) uncertified=%.4f (%d)\n', ...
    qp_certified_rate, qp_certified_count, soft_recovery_rate, soft_recovery_count_chk, ...
    uncertified_fallback_rate, uncertified_fallback_count);
fprintf('max QP residual=%.3e active_rate=%.4f avgJ=%.4f\n', ...
    max_recorded_qp_constraint, constraint_active_rate, mean(costJ(warm),'omitnan'));
fprintf('cross offline drop=%.3e ||Szo||=%.3e\n', drop_cross_rel_err, norm(Szo_xcov,'fro'));

%% Save
results_dir = fullfile(root_dir,'results');
if ~exist(results_dir,'dir'), mkdir(results_dir); end
out = struct();
out.schema_version = 'copyAC_open_problems_trial_v1';
out.trial = trial; out.stats = stats;
out.Ahat=Ahat; out.Bhat=Bhat; out.Phat=Phat; out.Rhat=Rhat;
out.Sigma_eps_used=Sigma_eps_use; out.Sigma_n=Sigma_n; out.Szo=Szo_xcov;
out.y=y; out.u=u; out.yhat=yhat; out.Rf=Rf;
out.cc_cert_level=cc_cert_level; out.soft_mode=soft_mode;
out.exitflag=exitflag; out.max_cc_violation=max_cc_violation; out.costJ=costJ;
out.MAE=MAE; out.RMSE=RMSE; out.Bias=Bias;
out.upper_violation_rate=upper_violation_rate; out.abs_violation_rate=abs_violation_rate;
out.qp_certified_rate=qp_certified_rate; out.soft_recovery_rate=soft_recovery_rate;
out.uncertified_fallback_rate=uncertified_fallback_rate;
out.max_recorded_qp_constraint=max_recorded_qp_constraint;
out.constraint_active_rate=constraint_active_rate;
out.sigma_obs_support_diag=last_sigma_obs_diag;
save(fullfile(results_dir,'copyAC_open_problems_trial_data.mat'),'-struct','out','-v7.3');

fid = fopen(fullfile(results_dir,'copyAC_open_problems_trial_metrics.txt'),'w');
fprintf(fid,'copyAC_open_problems_trial metrics\n');
fprintf(fid,'sigma_eps_mode %s\n', trial.sigma_eps_mode);
fprintf(fid,'use_cross_cov %d\n', trial.use_cross_cov);
fprintf(fid,'sigma_obs_mode %s\n', trial.sigma_obs_mode);
fprintf(fid,'enable_soft_recovery %d\n', trial.enable_soft_recovery);
fprintf(fid,'use_terminal_cost %d\n', trial.use_terminal_cost);
fprintf(fid,'input_residualize %d\n', trial.input_residualize);
fprintf(fid,'tracked_left_error %.12e\n', stats.tracked_left_error);
fprintf(fid,'tracked_right_error %.12e\n', stats.tracked_right_error);
fprintf(fid,'tracked_column_error %.12e\n', stats.tracked_column_error);
fprintf(fid,'dual_error %.12e\n', stats.dual_error);
fprintf(fid,'free_oblique_norm %.12e\n', stats.free_oblique_norm);
fprintf(fid,'free_noise_improvement %.12e\n', stats.free_noise_improvement);
fprintf(fid,'cross_drop_rel_err %.12e\n', drop_cross_rel_err);
fprintf(fid,'cross_Sigma_zo_fro %.12e\n', norm(Szo_xcov,'fro'));
fprintf(fid,'MAE %.12f %.12f\n', MAE(1), MAE(2));
fprintf(fid,'RMSE %.12f %.12f\n', RMSE(1), RMSE(2));
fprintf(fid,'upper_violation_rate %.12f %.12f\n', upper_violation_rate(1), upper_violation_rate(2));
fprintf(fid,'qp_certified_rate %.12f\n', qp_certified_rate);
fprintf(fid,'soft_recovery_rate %.12f\n', soft_recovery_rate);
fprintf(fid,'uncertified_fallback_rate %.12f\n', uncertified_fallback_rate);
fprintf(fid,'max_recorded_qp_constraint %.12e\n', max_recorded_qp_constraint);
fprintf(fid,'constraint_active_rate %.12f\n', constraint_active_rate);
fprintf(fid,'avg_costJ %.12f\n', mean(costJ(warm),'omitnan'));
fclose(fid);

%% Figure
fig = figure('Position',[40 40 2200 1600],'Color','w');
tlo = tiledlayout(fig,5,1,'TileSpacing','compact','Padding','compact');
t = 1:T_cl;
ax1=nexttile(tlo,1);
plot(ax1,t,Rf(1,:),'k--',t,y(1,:),'b',t,Rf(2,:),'k:',t,y(2,:),'r'); hold(ax1,'on');
yline(ax1,y_max,'m--'); grid(ax1,'on'); title(ax1,'Tracked outputs');
ax2=nexttile(tlo,2);
plot(ax2,t,estimated_sigma_eps,'r',t,estimated_sigma_obs,'b'); hold(ax2,'on');
plot(ax2,t,sigma_w_profile,'k--',t,sigma_e_profile,'m--');
plot(ax2,t,true_sigma_w,'Color',[0.3 0.3 0.3]);
plot(ax2,t,true_sigma_e,'Color',[0.6 0.1 0.6]);
grid(ax2,'on'); title(ax2,'Noise diagnostics');
ax3=nexttile(tlo,3); plot(ax3,t,u'); yline(ax3,u_min,'k--'); yline(ax3,u_max,'k--'); grid(ax3,'on'); title(ax3,'Inputs');
ax4=nexttile(tlo,4); plot(ax4,t,costJ,'k'); grid(ax4,'on'); title(ax4,'Cost J');
ax5=nexttile(tlo,5); plot(ax5,t,max_cc_violation,'k'); yline(ax5,0,'r--'); grid(ax5,'on');
title(ax5,sprintf('QP residual | qp=%.3f soft=%.3f uncert=%.3f',qp_certified_rate,soft_recovery_rate,uncertified_fallback_rate));
linkaxes([ax1 ax2 ax3 ax4 ax5],'x');
title(tlo,sprintf('copyAC open problems trial (eps=%s cross=%d obs=%s)',trial.sigma_eps_mode,trial.use_cross_cov,trial.sigma_obs_mode));
print(fig,fullfile(results_dir,'copyAC_open_problems_trial_fig'),'-dpng','-r140');
