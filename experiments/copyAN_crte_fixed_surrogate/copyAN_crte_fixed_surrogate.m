%% copyAN_crte_fixed_surrogate — paper CRTE spectral surrogate in the free complement
% Implements the executable fixed spectral surrogate from CRTE draft Sec 3.3.
% It is NOT the profiled teacher objective: surrogate uses a weighted Frobenius-
% normalized linear combination of S_yu, A_T, and a free-output residual noise
% proxy, then solves a symmetric generalized eigenproblem in the metric C_mu.
%
% Boundaries:
% - True sensor-noise Sigma_n is unknown and never passed in.
% - Each candidate rebuilds the complete metric-dual basis and refits VARX.
% - Selection is from a frozen (mu, alpha, beta) grid using validation-set
%   one-step NRMSE on the controlled outputs; closed-loop simulation uses
%   the chosen candidate rebuilt on all offline data.
% - Coverage is reported as a warm-horizon QP-success time-fraction proxy
%   only. It is NOT a chance-constraint probability certificate.
% - No claim of recursive feasibility, closed-loop stability, or global
%   surrogate-vs-teacher optimality is made.
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
rng(20260710,'twister');  % same plant seed as copyAA / copyAM for fair compare

%% Process-like plant and controller parameters
n = 6;
m = 3;
p = 30;
ell = 5;
tracked = [1 2];
T_off = 1500;
T_cl = 1200;
N = 18;

sw = 0.045;
se = 0.055;
sensor_noise_mode = 'heteroscedastic_correlated_plant_only';
noise_cycle = 400;
sw_min = 0.020; sw_max = 0.090;
se_min = 0.025; se_max = 0.100;
noise_phase_e = pi/3;
u_min = -3.0;
u_max =  3.0;
y_max = 2.00;
alpha_joint = 0.10;

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
end
for i = 3:p
    C(i,:) = C(i,:) / max(norm(C(i,:)), 1e-12);
end

%% Plant sensor-noise (for data generation only; algorithm never receives it)
sensor_rel = linspace(0.55,1.65,p)';
Corr_n = eye(p);
Corr_n(3,4) = 0.30; Corr_n(4,3) = 0.30;
Corr_n(5,6) = -0.22; Corr_n(6,5) = -0.22;
Corr_n(8,9) = 0.18; Corr_n(9,8) = 0.18;
D_n = diag(se*sensor_rel);
Sigma_n_plant = D_n*Corr_n*D_n;
L_n = chol(Sigma_n_plant,'lower');

%% Offline excitation
u_off = 1.20*randn(m,T_off);
x_off = zeros(n,T_off+1);
y_off = zeros(p,T_off);
for k = 1:T_off
    y_off(:,k) = C*x_off(:,k) + L_n*randn(p,1);
    x_off(:,k+1) = A*x_off(:,k) + B*u_off(:,k) + sw*randn(n,1);
end

%% CRTE candidate grid search + final refit
mu_grid = [0 0.25 0.5 0.75 1];
alpha_grid = [0 0.5 1];
beta_grid = [0 0.25 0.5];
[Ahat,Bhat,Phat,Rhat,Sigma_eps,stats] = crte_fixed_surrogate_varx( ...
    y_off,u_off,ell,tracked, ...
    struct('mu_grid',mu_grid,'alpha_grid',alpha_grid,'beta_grid',beta_grid, ...
          'val_fraction',0.25,'ridge',1e-8, ...
          'task_gate_fraction',0.10,'noise_gate_factor',2.0, ...
          'reach_gate_fraction',0.05,'reach_horizon',N));
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
model.Sigma_eps = (Sigma_eps+Sigma_eps')/2;
IminusPR = eye(p)-Phat*Rhat';
obs_off = IminusPR*(y_off-stats.y_mean);
obs_off = obs_off-mean(obs_off,2);
Sigma_obs_proxy = (obs_off*obs_off')/max(size(obs_off,2)-1,1);
Sigma_obs_proxy = (Sigma_obs_proxy+Sigma_obs_proxy')/2;
ridge_obs = 1e-8*max(trace(Sigma_obs_proxy)/max(p,1),1);
model.Sigma_obs = Sigma_obs_proxy+ridge_obs*eye(p);

opt.N = N;
opt.Q = zeros(p);
opt.Q(1,1) = 80;
opt.Q(2,2) = 80;
opt.Ru = 0.18*eye(m);
opt.u_min = u_min;
opt.u_max = u_max;
opt.H = zeros(numel(tracked),p);
opt.H(:,tracked) = eye(numel(tracked));
opt.h = y_max*ones(numel(tracked),1);
opt.alpha_joint = alpha_joint;

%% Closed-loop reference
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

%% Closed-loop SMPC
[sigma_w_profile,sigma_e_profile] = smooth_noise_profile( ...
    T_cl,noise_cycle,sw_min,sw_max,se_min,se_max,noise_phase_e);
x = zeros(n,1);
y = zeros(p,T_cl);
yhat = zeros(p,T_cl);
u = zeros(m,T_cl);
exitflag = zeros(1,T_cl);
max_cc_violation = nan(1,T_cl);
costJ = nan(1,T_cl);
estimated_sigma_eps = nan(1,T_cl);
estimated_sigma_obs = nan(1,T_cl);
true_sigma_w = nan(1,T_cl);
true_sigma_e = nan(1,T_cl);
noise_window = 40;
eps_buffer = zeros(ell,noise_window); eps_count = 0;
obs_buffer = zeros(p,noise_window); obs_count = 0;
infeasible_count = 0;

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
        Ewin = eps_buffer(:,1:eps_count);
        model.Sigma_eps = (Ewin*Ewin')/max(eps_count,1)+1e-8*eye(ell);
        model.Sigma_eps = (model.Sigma_eps+model.Sigma_eps')/2;
    end
    if obs_count >= 5
        Owin = obs_buffer(:,1:obs_count)-mean(obs_buffer(:,1:obs_count),2);
        sigma_obs_k = norm(Owin,'fro')/sqrt(max((p-ell)*(obs_count-1),1));
        model.Sigma_obs = max(sigma_obs_k^2,1e-8)*mean(diag(Sigma_obs_proxy))*eye(p)/p+1e-8*eye(p);
    end
    estimated_sigma_eps(k) = sqrt(trace(model.Sigma_eps)/ell);
    estimated_sigma_obs(k) = sqrt(trace(model.Sigma_obs)/p);
    rk = Rf(:,min(k+1,T_cl));
    try
        [~,ypred,U,info] = centered_smpc_step(yk,rk,model,opt);
        uk = U(1:m);
        yhat(:,k) = ypred;
        exitflag(k) = info.exitflag;
        max_cc_violation(k) = max(info.A_ch*U-info.b_ch);
        costJ(k) = info.cost;
    catch ME
        infeasible_count = infeasible_count+1;
        exitflag(k) = -1;
        uk = min(max(model.u_mean,u_min),u_max);
        if k>1, yhat(:,k)=yhat(:,k-1); else, yhat(:,k)=model.y_mean; end
        if infeasible_count<=3
            fprintf('QP fallback at k=%d: %s\n',k,ME.message);
        end
    end
    u(:,k) = uk;
    wk = sigma_w_profile(k)*randn(n,1);
    true_sigma_w(k) = norm(wk)/sqrt(n);
    x = A*x+B*uk+wk;
    z_prev = zk_now;
end

%% Metrics
warm = 151:T_cl;
err = y(tracked,warm)-Rf(tracked,warm);
MAE = mean(abs(err),2);
RMSE = sqrt(mean(err.^2,2));
Bias = mean(err,2);
upper_violation_count = sum(y(tracked,:)>y_max,2);
upper_violation_rate = upper_violation_count/T_cl;
abs_violation_count = sum(abs(y(tracked,:))>y_max,2);
abs_violation_rate = abs_violation_count/T_cl;
u_rms = sqrt(mean(u.^2,2));
u_sat_rate = mean(abs(u-u_min)<1e-8 | abs(u-u_max)<1e-8,2);
qp_success_rate = mean(exitflag>0);
max_recorded_qp_constraint = max(max_cc_violation(~isnan(max_cc_violation)));
if isempty(max_recorded_qp_constraint), max_recorded_qp_constraint = NaN; end
constraint_active_rate = mean(max_cc_violation(warm)>-1e-3,'omitnan');
cover = mean(exitflag(warm)>0);

fprintf('\ncopyAN CRTE fixed spectral surrogate (free complement only) results\n');
fprintf('selected mu=%.3f alpha=%.3f beta=%.3f validation NRMSE=%.4f\n', ...
    stats.selected_mu,stats.selected_alpha,stats.selected_beta,stats.selected_validation_nrmse);
fprintf('num candidates=%d valid after gates=%d\n', numel(stats.candidates), sum(stats.valid_candidates));
Ec = zeros(p,numel(tracked)); Ec(tracked,:) = eye(numel(tracked));
stats.tracked_projection_error = norm(Phat*Phat'*Ec-Ec,'fro');
stats.tracked_right_error = norm(Phat*Rhat'*Ec-Ec,'fro');
stats.tracked_left_error = norm(Ec'*Phat*Rhat'-Ec','fro');
stats.tracked_column_error = norm(Rhat(:,1:numel(tracked))-Ec,'fro');
fprintf('tracked projection error = %.3e\n', stats.tracked_projection_error);
fprintf('tracked right = %.3e, left = %.3e, columns = %.3e, dual = %.3e\n', ...
    stats.tracked_right_error,stats.tracked_left_error, ...
    stats.tracked_column_error,stats.dual_error);
fprintf('4-piece dual completion errors (R Pbar Rbar P) max = %.3e\n', ...
    max(stats.dual_errors_4piece));
fprintf('spectral radius max|eig(Ahat)|=%.4f, Pi idempotency err=%.3e\n', ...
    stats.spectral_radius,stats.Pi_idempotency_err);
fprintf('reconstruction residual = %.3f\n', stats.reconstruction_residual);
fprintf('MAE  = [%.4f %.4f]\n',MAE(1),MAE(2));
fprintf('RMSE = [%.4f %.4f]\n',RMSE(1),RMSE(2));
fprintf('Bias = [%.4f %.4f]\n',Bias(1),Bias(2));
fprintf('upper violation rate = [%.4f %.4f]\n',upper_violation_rate(1),upper_violation_rate(2));
fprintf('abs violation rate   = [%.4f %.4f]\n',abs_violation_rate(1),abs_violation_rate(2));
fprintf('QP success rate = %.4f, fallback count = %d, max logged QP constraint = %.3e\n', ...
    qp_success_rate,infeasible_count,max_recorded_qp_constraint);
fprintf('constraint active rate = %.4f\n',constraint_active_rate);
fprintf('cover (warm QP-success time fraction proxy) = %.4f\n',cover);
fprintf('avg J = %.4f, estimated sigma_eps = %.4f, estimated sigma_obs = %.4f, realized sigma_w = %.4f, realized sigma_e = %.4f\n', ...
    mean(costJ(warm),'omitnan'), mean(estimated_sigma_eps(warm),'omitnan'), ...
    mean(estimated_sigma_obs(warm),'omitnan'), mean(true_sigma_w(warm)), mean(true_sigma_e(warm)));
fprintf('u RMS = [%.4f %.4f %.4f], saturation rate = [%.4f %.4f %.4f]\n', ...
    u_rms(1),u_rms(2),u_rms(3), u_sat_rate(1),u_sat_rate(2),u_sat_rate(3));

%% Save artifacts
results_dir = fullfile(fileparts(mfilename('fullpath')),'results');
if ~exist(results_dir,'dir'), mkdir(results_dir); end

out.schema_version = 'copyAN_crte_fixed_surrogate_v1';
out.description = ['fixed CRTE spectral surrogate in free complement only; ' ...
    'no oracle Sigma_n; selection by validation NRMSE'];
out.A_true = A; out.B_true = B; out.C_true = C;
out.tracked = tracked; out.n = n; out.m = m; out.p = p; out.ell = ell;
out.T_off = T_off; out.T_cl = T_cl; out.N = N;
out.sw = sw; out.se = se; out.sensor_noise_mode = sensor_noise_mode;
out.Sigma_n_plant = Sigma_n_plant;
out.sensor_rel = sensor_rel; out.Corr_n = Corr_n;
out.sigma_w_profile = sigma_w_profile; out.sigma_e_profile = sigma_e_profile;
out.noise_cycle = noise_cycle;
out.sw_min = sw_min; out.sw_max = sw_max;
out.se_min = se_min; out.se_max = se_max; out.noise_phase_e = noise_phase_e;
out.u_min = u_min; out.u_max = u_max; out.y_max = y_max; out.alpha_joint = alpha_joint;
out.mu_grid = mu_grid; out.alpha_grid = alpha_grid; out.beta_grid = beta_grid;
out.x_off = x_off; out.y_off = y_off; out.u_off = u_off;
out.Ahat = Ahat; out.Bhat = Bhat; out.Phat = Phat; out.Rhat = Rhat;
out.Sigma_eps = Sigma_eps; out.stats = stats;
out.model = model; out.opt = opt;
out.Rf = Rf; out.y = y; out.yhat = yhat; out.u = u;
out.exitflag = exitflag; out.max_cc_violation = max_cc_violation;
out.costJ = costJ;
out.estimated_sigma_eps = estimated_sigma_eps;
out.estimated_sigma_obs = estimated_sigma_obs;
out.true_sigma_w = true_sigma_w; out.true_sigma_e = true_sigma_e;
out.noise_window = noise_window;
out.constraint_active_rate = constraint_active_rate; out.cover = cover;
out.MAE = MAE; out.RMSE = RMSE; out.Bias = Bias;
out.upper_violation_count = upper_violation_count;
out.upper_violation_rate = upper_violation_rate;
out.abs_violation_count = abs_violation_count;
out.abs_violation_rate = abs_violation_rate;
out.u_rms = u_rms; out.u_sat_rate = u_sat_rate;
out.qp_success_rate = qp_success_rate; out.infeasible_count = infeasible_count;
save(fullfile(results_dir,'copyAN_crte_fixed_surrogate_data.mat'),'-struct','out','-v7.3');

fid = fopen(fullfile(results_dir,'copyAN_crte_fixed_surrogate_metrics.txt'),'w');
fprintf(fid,'copyAN_crte_fixed_surrogate metrics\n');
fprintf(fid,'metric_name %s\n','CRTE fixed spectral surrogate (free complement only)');
fprintf(fid,'uses_true_Sigma_n %d\n',0);
fprintf(fid,'selected_mu %.12f\n',stats.selected_mu);
fprintf(fid,'selected_alpha %.12f\n',stats.selected_alpha);
fprintf(fid,'selected_beta %.12f\n',stats.selected_beta);
fprintf(fid,'selected_validation_nrmse %.12f\n',stats.selected_validation_nrmse);
fprintf(fid,'num_candidates %d\n',numel(stats.candidates));
fprintf(fid,'valid_candidates %d\n',sum(stats.valid_candidates));
fprintf(fid,'tracked_right_error %.12e\n',stats.tracked_right_error);
fprintf(fid,'tracked_left_error %.12e\n',stats.tracked_left_error);
fprintf(fid,'tracked_column_error %.12e\n',stats.tracked_column_error);
fprintf(fid,'dual_error %.12e\n',stats.dual_error);
fprintf(fid,'dual_error_RPbar %.12e\n',stats.dual_errors_4piece(2));
fprintf(fid,'dual_error_RbarP %.12e\n',stats.dual_errors_4piece(3));
fprintf(fid,'dual_error_RbarPbar %.12e\n',stats.dual_errors_4piece(4));
fprintf(fid,'dual_basis_completion %d\n',stats.dual_basis_completion);
fprintf(fid,'Pi_idempotency_err %.12e\n',stats.Pi_idempotency_err);
fprintf(fid,'pr_asymmetry %.12e\n',stats.pr_asymmetry);
fprintf(fid,'spectral_radius %.12f\n',stats.spectral_radius);
fprintf(fid,'reconstruction_residual %.12f\n',stats.reconstruction_residual);
fprintf(fid,'MAE %.12f %.12f\n',MAE(1),MAE(2));
fprintf(fid,'RMSE %.12f %.12f\n',RMSE(1),RMSE(2));
fprintf(fid,'Bias %.12f %.12f\n',Bias(1),Bias(2));
fprintf(fid,'upper_violation_rate %.12f %.12f\n',upper_violation_rate(1),upper_violation_rate(2));
fprintf(fid,'abs_violation_rate %.12f %.12f\n',abs_violation_rate(1),abs_violation_rate(2));
fprintf(fid,'qp_success_rate %.12f\n',qp_success_rate);
fprintf(fid,'infeasible_count %d\n',infeasible_count);
fprintf(fid,'max_recorded_qp_constraint %.12e\n',max_recorded_qp_constraint);
fprintf(fid,'constraint_active_rate %.12f\n',constraint_active_rate);
fprintf(fid,'cover_warm_qp_success_time_fraction %.12f\n',cover);
fprintf(fid,'avg_costJ %.12f\n',mean(costJ(warm),'omitnan'));
fprintf(fid,'estimated_sigma_eps_mean %.12f\n',mean(estimated_sigma_eps(warm),'omitnan'));
fprintf(fid,'estimated_sigma_obs_mean %.12f\n',mean(estimated_sigma_obs(warm),'omitnan'));
fprintf(fid,'realized_sigma_w_mean %.12f\n',mean(true_sigma_w(warm)));
fprintf(fid,'realized_sigma_e_mean %.12f\n',mean(true_sigma_e(warm)));
fprintf(fid,'u_rms %.12f %.12f %.12f\n',u_rms(1),u_rms(2),u_rms(3));
fprintf(fid,'u_sat_rate %.12f %.12f %.12f\n',u_sat_rate(1),u_sat_rate(2),u_sat_rate(3));
fclose(fid);

%% Figure
fig = figure('Position',[50 50 2400 1800], 'Color', 'w');
t = 1:T_cl;
tlo = tiledlayout(fig, 5, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tlo,1);
plot(ax1,t,Rf(1,:),'k--','LineWidth',1.2); hold(ax1,'on');
plot(ax1,t,y(1,:),'b','LineWidth',0.8);
plot(ax1,t,Rf(2,:),'Color',[0.20 0.20 0.20],'LineStyle',':','LineWidth',1.4);
plot(ax1,t,y(2,:),'Color',[0.85 0.20 0.10],'LineWidth',0.8);
yline(ax1,y_max,'m--','y_{max}','LabelHorizontalAlignment','left');
for s=1:5, xline(ax1,(s-1)*seg_len+1,'Color',[.75 .75 .75],'HandleVisibility','off'); end
grid(ax1,'on'); ylabel(ax1,'quality outputs');
title(ax1,'Tracked quality variables y_1,y_2');
legend(ax1,{'r_1','y_1','r_2','y_2','upper chance limit'},'Location','eastoutside');

ax2 = nexttile(tlo,2);
plot(ax2,t,estimated_sigma_eps,'Color',[0.85 0.20 0.10],'LineWidth',0.8); hold(ax2,'on');
plot(ax2,t,estimated_sigma_obs,'Color',[0.15 0.45 0.85],'LineWidth',0.8);
plot(ax2,t,sigma_w_profile,'k--','LineWidth',1.4);
plot(ax2,t,sigma_e_profile,'m--','LineWidth',1.4);
plot(ax2,t,true_sigma_w,'Color',[0.25 0.25 0.25],'LineWidth',0.55);
plot(ax2,t,true_sigma_e,'Color',[0.65 0.10 0.65],'LineWidth',0.55);
for s=1:5, xline(ax2,(s-1)*seg_len+1,'Color',[.75 .75 .75],'HandleVisibility','off'); end
grid(ax2,'on'); ylabel(ax2,'sigma');
title(ax2,sprintf('Noise diagnostics (cycle=%d, window=%d)',noise_cycle,noise_window));
legend(ax2,{'est eps','est obs','sw env','se env','real w','real e'},'Location','eastoutside');

ax3 = nexttile(tlo,3);
plot(ax3,t,u','LineWidth',0.8); hold(ax3,'on');
yline(ax3,u_min,'k--'); yline(ax3,u_max,'k--');
for s=1:5, xline(ax3,(s-1)*seg_len+1,'Color',[.75 .75 .75],'HandleVisibility','off'); end
grid(ax3,'on'); ylabel(ax3,'u');
title(ax3,sprintf('Manipulated vars RMS=[%.2f %.2f %.2f]',u_rms(1),u_rms(2),u_rms(3)));

ax4 = nexttile(tlo,4);
plot(ax4,t,costJ,'Color',[0.30 0.30 0.30],'LineWidth',0.8);
for s=1:5, xline(ax4,(s-1)*seg_len+1,'Color',[.75 .75 .75],'HandleVisibility','off'); end
grid(ax4,'on'); ylabel(ax4,'J');
title(ax4,sprintf('MPC cost J (avg=%.2f)',mean(costJ(warm),'omitnan')));

ax5 = nexttile(tlo,5);
plot(ax5,t,max_cc_violation,'k','LineWidth',0.8); hold(ax5,'on');
yline(ax5,0,'r--');
for s=1:5, xline(ax5,(s-1)*seg_len+1,'Color',[.75 .75 .75],'HandleVisibility','off'); end
grid(ax5,'on'); xlabel(ax5,'time step'); ylabel(ax5,'max(AU-b)');
title(ax5,sprintf('QP residual: success=%.3f cover=%.3f active=%.3f fallbacks=%d', ...
    qp_success_rate,cover,constraint_active_rate,infeasible_count));

linkaxes([ax1 ax2 ax3 ax4 ax5],'x');
title(tlo,sprintf('copyAN CRTE fixed surrogate: mu=%.2f alpha=%.2f beta=%.2f val NRMSE=%.3f (valid %d/%d)', ...
    stats.selected_mu,stats.selected_alpha,stats.selected_beta, ...
    stats.selected_validation_nrmse,sum(stats.valid_candidates),numel(stats.candidates)));
print(fig,fullfile(results_dir,'copyAN_crte_fixed_surrogate_fig'),'-dpng','-r160');
close(fig);