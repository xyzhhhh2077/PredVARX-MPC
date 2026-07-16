%% copyY_no_whitening_direct_update — copyX without whitening normalization
% Synthetic industrial-process validation for high-dimensional process data.
% copyX and main/ are untouched. This version uses direct eigenspace
% replacement in raw centered coordinates, without Y*=Y*U*D^(-1/2).
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
rng(20260710,'twister');

%% Process-like plant and controller parameters
n = 6;              % true latent process states
m = 3;              % manipulated variables
p = 30;             % high-dimensional sensors/quality measurements
ell = 5;            % reduced latent coordinates for controller
tracked = [1 2];    % controlled quality variables
T_off = 1500;
T_cl = 1200;
N = 18;

sw = 0.045;         % offline process disturbance scale (stationary identification data)
se = 0.055;         % offline sensor noise scale (stationary identification data)
% Closed-loop extension: smooth nonstationary standard-deviation envelopes.
noise_cycle = 400;
sw_min = 0.020; sw_max = 0.090;
se_min = 0.025; se_max = 0.100;
noise_phase_e = pi/3;
u_min = -3.0;
u_max =  3.0;
y_max = 2.00;        % stress-test quality limit requested by user
alpha_joint = 0.10;

% Stable, slow process dynamics with cross-coupling, typical of process plants.
A = diag([0.94 0.88 0.78 0.64 0.50 0.35]);
A(1,2) = 0.10; A(2,3) = -0.06; A(3,4) = 0.05; A(4,5) = 0.04;
B = [0.34 -0.10  0.05;
     0.12  0.28 -0.06;
     0.05  0.12  0.24;
    -0.05  0.06  0.18;
     0.02 -0.10  0.14;
     0.08  0.02 -0.08];

% Observation matrix: first two outputs are quality variables with direct coverage.
C = zeros(p,n);
C(1,1) = 1.00; C(1,3) = 0.16;    % quality y1
C(2,2) = 1.00; C(2,4) = -0.12;   % quality y2
for i = 3:p
    C(i,:) = 0.45*randn(1,n);
end
% Normalize non-quality sensor rows to avoid ill-scaled measurements.
for i = 3:p
    C(i,:) = C(i,:) / max(norm(C(i,:)), 1e-12);
end

%% Offline excitation data
u_off = 1.20*randn(m,T_off);
x_off = zeros(n,T_off+1);
y_off = zeros(p,T_off);
for k = 1:T_off
    y_off(:,k) = C*x_off(:,k) + se*randn(p,1);
    x_off(:,k+1) = A*x_off(:,k) + B*u_off(:,k) + sw*randn(n,1);
end

%% Identify reduced control-aware PredVARX model
% Held-out screening selected the smallest genuinely oblique setting.  The
% full covariance-weighted dual (alpha=1) was numerically too aggressive.
oblique_alpha = 0.02;
[Ahat,Bhat,Phat,Rhat,Sigma_eps,stats] = control_aware_direct_update_varx(y_off,u_off,ell,tracked,oblique_alpha);
stats.true_A_eigs = eig(A);
stats.id_A_eigs = eig(Ahat);
stats.ell = ell;
stats.tracked = tracked;

model.A = Ahat;
model.B = Bhat;
model.P = Phat;
model.R = Rhat;
model.y_mean = stats.y_mean;
model.u_mean = stats.u_mean;
model.Sigma_eps = (Sigma_eps + Sigma_eps')/2;
model.Sigma_obs = se^2 * eye(p);

opt.N = N;
opt.Q = zeros(p);
opt.Q(1,1) = 80;
opt.Q(2,2) = 80;
opt.Ru = 0.18 * eye(m);
opt.u_min = u_min;
opt.u_max = u_max;
opt.H = zeros(numel(tracked),p);
opt.H(:,tracked) = eye(numel(tracked));
opt.h = y_max * ones(numel(tracked),1);
opt.alpha_joint = alpha_joint;

%% Closed-loop reference: quality grade transitions / batch phases
Rf = zeros(p,T_cl);
levels = [0.25 1.50 0.65 1.85 0.45; 0.35 1.25 1.75 0.80 1.55];
seg_len = floor(T_cl/5);
for s = 1:5
    ix = (s-1)*seg_len+1:min(s*seg_len,T_cl);
    Rf(1,ix) = levels(1,s);
    Rf(2,ix) = levels(2,s);
end
if seg_len*5 < T_cl
    Rf(1,seg_len*5+1:end) = levels(1,end);
    Rf(2,seg_len*5+1:end) = levels(2,end);
end

%% Closed-loop SMPC simulation
[sigma_w_profile, sigma_e_profile] = smooth_noise_profile( ...
    T_cl, noise_cycle, sw_min, sw_max, se_min, se_max, noise_phase_e);
x = zeros(n,1);
y = zeros(p,T_cl);
yhat = zeros(p,T_cl);
u = zeros(m,T_cl);
exitflag = zeros(1,T_cl);
max_cc_violation = nan(1,T_cl);
costJ = nan(1,T_cl);
estimated_sigma_eps = nan(1,T_cl);
estimated_sigma_obs = nan(1,T_cl);
true_sigma_w = nan(1,T_cl);   % realized process-noise RMS ||w_k||/sqrt(n)
true_sigma_e = nan(1,T_cl);   % realized measurement/input-to-estimator noise RMS ||v_k||/sqrt(p)
noise_window = 40;
eps_buffer = zeros(ell,noise_window); eps_count = 0;
obs_buffer = zeros(p,noise_window); obs_count = 0;
IminusPR = eye(p)-Phat*Rhat';
infeasible_count = 0;

for k = 1:T_cl
    vk = sigma_e_profile(k)*randn(p,1);
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
        model.Sigma_eps = (Ewin*Ewin')/max(eps_count-1,1)+1e-8*eye(ell);
        model.Sigma_eps = (model.Sigma_eps+model.Sigma_eps')/2;
    end
    if obs_count >= 5
        Owin = obs_buffer(:,1:obs_count)-mean(obs_buffer(:,1:obs_count),2);
        sigma_obs_k = norm(Owin,'fro')/sqrt(max((p-ell)*(obs_count-1),1));
        model.Sigma_obs = max(sigma_obs_k^2,1e-8)*eye(p);
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

    catch ME
        infeasible_count = infeasible_count + 1;
        exitflag(k) = -1;
        % Safe fallback: hold nearest mean input inside bounds.
        uk = min(max(model.u_mean, u_min), u_max);
        if k > 1
            yhat(:,k) = yhat(:,k-1);
        else
            yhat(:,k) = model.y_mean;
        end
        if infeasible_count <= 3
            fprintf('QP fallback at k=%d: %s\n', k, ME.message);
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
upper_violation_count = sum(y(tracked,:) > y_max, 2);
upper_violation_rate = upper_violation_count / T_cl;
abs_violation_count = sum(abs(y(tracked,:)) > y_max, 2);
abs_violation_rate = abs_violation_count / T_cl;
u_rms = sqrt(mean(u.^2,2));
u_sat_rate = mean(abs(u - u_min) < 1e-8 | abs(u - u_max) < 1e-8, 2);
qp_success_rate = mean(exitflag > 0);
max_recorded_qp_constraint = max(max_cc_violation(~isnan(max_cc_violation)));
if isempty(max_recorded_qp_constraint), max_recorded_qp_constraint = NaN; end
constraint_active_rate = mean(max_cc_violation(warm) > -1e-3, 'omitnan');

fprintf('\ncopyY no-whitening direct-update PredVARX-SMPC results\n');
fprintf('normalization applied = %d, update rule = %s\n', stats.normalization_applied, stats.update_rule);
fprintf('tracked projection error = %.3e\n', stats.tracked_projection_error);
fprintf('tracked oblique error = %.3e, dual error = %.3e, PR asymmetry = %.3e\n', stats.tracked_oblique_error, stats.dual_error, stats.pr_asymmetry);
fprintf('reconstruction residual  = %.3f\n', stats.reconstruction_residual);
fprintf('MAE  = [%.4f %.4f]\n', MAE(1), MAE(2));
fprintf('RMSE = [%.4f %.4f]\n', RMSE(1), RMSE(2));
fprintf('Bias = [%.4f %.4f]\n', Bias(1), Bias(2));
fprintf('upper violation rate = [%.4f %.4f]\n', upper_violation_rate(1), upper_violation_rate(2));
fprintf('abs violation rate   = [%.4f %.4f]\n', abs_violation_rate(1), abs_violation_rate(2));
fprintf('QP success rate = %.4f, fallback count = %d, max logged QP constraint = %.3e\n', qp_success_rate, infeasible_count, max_recorded_qp_constraint);
fprintf('constraint active rate = %.4f (residual > -1e-3)\n',constraint_active_rate);
fprintf('avg J = %.4f, estimated sigma_eps = %.4f, estimated sigma_obs = %.4f, realized sigma_w = %.4f, realized sigma_e = %.4f\n', mean(costJ(warm),'omitnan'), mean(estimated_sigma_eps(warm),'omitnan'), mean(estimated_sigma_obs(warm),'omitnan'), mean(true_sigma_w(warm)), mean(true_sigma_e(warm)));
fprintf('smooth sigma ranges: sw=[%.3f %.3f], se=[%.3f %.3f], cycle=%d\n', min(sigma_w_profile), max(sigma_w_profile), min(sigma_e_profile), max(sigma_e_profile), noise_cycle);
fprintf('u RMS = [%.4f %.4f %.4f], saturation rate = [%.4f %.4f %.4f]\n', u_rms(1), u_rms(2), u_rms(3), u_sat_rate(1), u_sat_rate(2), u_sat_rate(3));

%% Save artifacts
results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(results_dir,'dir'), mkdir(results_dir); end

out.schema_version = 'copyY_no_whitening_direct_update_v1';
out.description = 'copyX-derived PredVARX with no whitening normalization and direct eigenspace replacement';
out.A_true = A; out.B_true = B; out.C_true = C;
out.tracked = tracked; out.n = n; out.m = m; out.p = p; out.ell = ell;
out.T_off = T_off; out.T_cl = T_cl; out.N = N;
out.sw = sw; out.se = se; out.oblique_alpha = oblique_alpha; out.sigma_w_profile = sigma_w_profile; out.sigma_e_profile = sigma_e_profile;
out.noise_cycle = noise_cycle; out.sw_min = sw_min; out.sw_max = sw_max; out.se_min = se_min; out.se_max = se_max; out.noise_phase_e = noise_phase_e;
out.u_min = u_min; out.u_max = u_max; out.y_max = y_max; out.alpha_joint = alpha_joint;
out.x_off = x_off; out.y_off = y_off; out.u_off = u_off;
out.Ahat = Ahat; out.Bhat = Bhat; out.Phat = Phat; out.Rhat = Rhat; out.Sigma_eps = Sigma_eps; out.stats = stats;
out.model = model; out.opt = opt;
out.Rf = Rf; out.y = y; out.yhat = yhat; out.u = u;
out.exitflag = exitflag; out.max_cc_violation = max_cc_violation;
out.costJ = costJ; out.estimated_sigma_eps = estimated_sigma_eps; out.estimated_sigma_obs = estimated_sigma_obs;
out.true_sigma_w = true_sigma_w; out.true_sigma_e = true_sigma_e;
out.noise_window = noise_window; out.constraint_active_rate = constraint_active_rate;
out.MAE = MAE; out.RMSE = RMSE; out.Bias = Bias;
out.upper_violation_count = upper_violation_count; out.upper_violation_rate = upper_violation_rate;
out.abs_violation_count = abs_violation_count; out.abs_violation_rate = abs_violation_rate;
out.u_rms = u_rms; out.u_sat_rate = u_sat_rate; out.qp_success_rate = qp_success_rate; out.infeasible_count = infeasible_count;
save(fullfile(results_dir,'copyY_no_whitening_direct_update_data.mat'), '-struct', 'out', '-v7.3');

fid = fopen(fullfile(results_dir,'copyY_no_whitening_direct_update_metrics.txt'), 'w');
fprintf(fid, 'copyY_no_whitening_direct_update metrics\n');
fprintf(fid, 'normalization_applied %d\n', stats.normalization_applied);
fprintf(fid, 'update_rule %s\n', stats.update_rule);
fprintf(fid, 'tracked_projection_error %.12e\n', stats.tracked_projection_error);
fprintf(fid, 'tracked_oblique_error %.12e\n', stats.tracked_oblique_error);
fprintf(fid, 'dual_error %.12e\n', stats.dual_error);
fprintf(fid, 'pr_asymmetry %.12e\n', stats.pr_asymmetry);
fprintf(fid, 'cond_dual_gram %.12e\n', stats.cond_dual_gram);
fprintf(fid, 'reconstruction_residual %.12f\n', stats.reconstruction_residual);
fprintf(fid, 'MAE %.12f %.12f\n', MAE(1), MAE(2));
fprintf(fid, 'RMSE %.12f %.12f\n', RMSE(1), RMSE(2));
fprintf(fid, 'Bias %.12f %.12f\n', Bias(1), Bias(2));
fprintf(fid, 'upper_violation_rate %.12f %.12f\n', upper_violation_rate(1), upper_violation_rate(2));
fprintf(fid, 'abs_violation_rate %.12f %.12f\n', abs_violation_rate(1), abs_violation_rate(2));
fprintf(fid, 'qp_success_rate %.12f\n', qp_success_rate);
fprintf(fid, 'infeasible_count %d\n', infeasible_count);
fprintf(fid, 'max_recorded_qp_constraint %.12e\n', max_recorded_qp_constraint);
fprintf(fid, 'constraint_active_rate %.12f\n', constraint_active_rate);
fprintf(fid, 'avg_costJ %.12f\n', mean(costJ(warm),'omitnan'));
fprintf(fid, 'estimated_sigma_eps_mean %.12f\n', mean(estimated_sigma_eps(warm),'omitnan'));
fprintf(fid, 'estimated_sigma_obs_mean %.12f\n', mean(estimated_sigma_obs(warm),'omitnan'));
fprintf(fid, 'realized_sigma_w_mean %.12f\n', mean(true_sigma_w(warm)));
fprintf(fid, 'realized_sigma_e_mean %.12f\n', mean(true_sigma_e(warm)));
fprintf(fid, 'sigma_w_profile_range %.12f %.12f\n', min(sigma_w_profile), max(sigma_w_profile));
fprintf(fid, 'sigma_e_profile_range %.12f %.12f\n', min(sigma_e_profile), max(sigma_e_profile));
fprintf(fid, 'u_rms %.12f %.12f %.12f\n', u_rms(1), u_rms(2), u_rms(3));
fprintf(fid, 'u_sat_rate %.12f %.12f %.12f\n', u_sat_rate(1), u_sat_rate(2), u_sat_rate(3));
fclose(fid);

%% Figure: paper-style 5 rows x 1 column, aligned with baseline diagnostics
fig = figure('Position',[50 50 2400 1800], 'Color', 'w');
t = 1:T_cl;
tlo = tiledlayout(fig, 5, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tlo, 1);
plot(ax1, t, Rf(1,:), 'k--', 'LineWidth', 1.2); hold(ax1, 'on');
plot(ax1, t, y(1,:), 'b', 'LineWidth', 0.8);
plot(ax1, t, Rf(2,:), 'Color', [0.20 0.20 0.20], 'LineStyle', ':', 'LineWidth', 1.4);
plot(ax1, t, y(2,:), 'Color', [0.85 0.20 0.10], 'LineWidth', 0.8);
yline(ax1, y_max, 'm--', 'y_{max}', 'LabelHorizontalAlignment', 'left');
for s = 1:5, xline(ax1, (s-1)*seg_len+1, 'Color', [.75 .75 .75], 'HandleVisibility', 'off'); end
grid(ax1, 'on'); ylabel(ax1, 'quality outputs');
title(ax1, 'Tracked quality variables y_1,y_2');
legend(ax1, {'r_1','y_1','r_2','y_2','upper chance limit'}, 'Location', 'eastoutside');

ax2 = nexttile(tlo, 2);
plot(ax2, t, estimated_sigma_eps, 'Color', [0.85 0.20 0.10], 'LineWidth', 0.8, 'DisplayName', sprintf('estimated sigma_{eps} (avg=%.3f)', mean(estimated_sigma_eps(warm),'omitnan'))); hold(ax2, 'on');
plot(ax2, t, estimated_sigma_obs, 'Color', [0.15 0.45 0.85], 'LineWidth', 0.8, 'DisplayName', sprintf('estimated sigma_{obs} (avg=%.3f)', mean(estimated_sigma_obs(warm),'omitnan')));
plot(ax2, t, sigma_w_profile, 'k--', 'LineWidth', 1.4, 'DisplayName', 'specified sigma_w(k) envelope');
plot(ax2, t, sigma_e_profile, 'm--', 'LineWidth', 1.4, 'DisplayName', 'specified sigma_e(k) envelope');
plot(ax2, t, true_sigma_w, 'Color', [0.25 0.25 0.25], 'LineWidth', 0.55, 'DisplayName', sprintf('realized process noise RMS (avg=%.3f)', mean(true_sigma_w(warm))));
plot(ax2, t, true_sigma_e, 'Color', [0.65 0.10 0.65], 'LineWidth', 0.55, 'DisplayName', sprintf('realized measurement/input noise RMS (avg=%.3f)', mean(true_sigma_e(warm))));
for s = 1:5, xline(ax2, (s-1)*seg_len+1, 'Color', [.75 .75 .75], 'HandleVisibility', 'off'); end
grid(ax2, 'on'); ylabel(ax2, 'sigma');
title(ax2, sprintf('Smooth nonstationary noise (cycle=%d) and %d-step rolling estimates',noise_cycle,noise_window));
legend(ax2, 'Location', 'eastoutside');

ax3 = nexttile(tlo, 3);
plot(ax3, t, u', 'LineWidth', 0.8); hold(ax3, 'on');
yline(ax3, u_min, 'k--', 'u bounds', 'LabelHorizontalAlignment', 'left');
yline(ax3, u_max, 'k--', 'HandleVisibility', 'off');
for s = 1:5, xline(ax3, (s-1)*seg_len+1, 'Color', [.75 .75 .75], 'HandleVisibility', 'off'); end
grid(ax3, 'on'); ylabel(ax3, 'u');
title(ax3, sprintf('Manipulated variables: RMS=[%.2f %.2f %.2f], saturation=[%.2f %.2f %.2f]', u_rms(1), u_rms(2), u_rms(3), u_sat_rate(1), u_sat_rate(2), u_sat_rate(3)));
legend(ax3, {'u_1','u_2','u_3','u_{min/max}'}, 'Location', 'eastoutside');

ax4 = nexttile(tlo, 4);
plot(ax4, t, costJ, 'Color', [0.30 0.30 0.30], 'LineWidth', 0.8, 'DisplayName', sprintf('J (avg=%.2f)', mean(costJ(warm),'omitnan'))); hold(ax4, 'on');
for s = 1:5, xline(ax4, (s-1)*seg_len+1, 'Color', [.75 .75 .75], 'HandleVisibility', 'off'); end
grid(ax4, 'on'); ylabel(ax4, 'J');
title(ax4, 'Centered MPC raw cost J: Q-tracking + R on (U-U_0)');
legend(ax4, 'Location', 'eastoutside');

ax5 = nexttile(tlo, 5);
plot(ax5, t, max_cc_violation, 'k', 'LineWidth', 0.8); hold(ax5, 'on');
yline(ax5, 0, 'r--', 'constraint boundary', 'LabelHorizontalAlignment', 'left');
for s = 1:5, xline(ax5, (s-1)*seg_len+1, 'Color', [.75 .75 .75], 'HandleVisibility', 'off'); end
grid(ax5, 'on'); xlabel(ax5, 'time step'); ylabel(ax5, 'max(AU-b)');
title(ax5, sprintf('QP chance-constraint residual: success=%.3f, active=%.3f, fallbacks=%d, max=%.2e', qp_success_rate, constraint_active_rate, infeasible_count, max_recorded_qp_constraint));

linkaxes([ax1 ax2 ax3 ax4 ax5], 'x');
title(tlo, sprintf('copyY no-whitening direct-update PredVARX-SMPC (p=%d, ell=%d, coverage %.1e)', p, ell, stats.tracked_oblique_error));
print(fig, fullfile(results_dir,'copyY_no_whitening_direct_update_fig'), '-dpng', '-r160');
