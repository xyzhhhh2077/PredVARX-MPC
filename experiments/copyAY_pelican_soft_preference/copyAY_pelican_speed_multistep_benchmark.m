%% copyAY_pelican_speed_multistep_benchmark
% FAIR multi-step benchmark vs Mohajerin & Waslander ICRA 2018.
% The paper predicts TRANSLATIONAL VELOCITY and BODY RATES (not position):
% "hybrid model produces predictions, over 1.9 second, which remain within
%  9 cm/s and 0.12 rad/s of the measured translational and rotational
%  velocities, with 99% confidence on the test dataset" (speed max ~4 m/s).
% Protocol: initialize from a measurement, feed the RECORDED input sequence,
% roll the identified model forward H=190 steps (1.9 s at 100 Hz), report the
% per-horizon error distribution (mean + p99) exactly as the paper does.
%
% NOTE (2026-08-02): the position-state model (copyAY main) scores R2 0.98
% one-step but collapses in multi-step (velocity error ~100 cm/s ~= persistence)
% because position is an INTEGRATED quantity and a 5-dim latent cannot carry
% fast dynamics.  This script uses the SPEED STATE y=[Vel(3);pqr(3);Motors(4)]
% (the paper's own prediction targets; C_y conditioning 6.2 non-singular since
% Vel/pqr never appear together with their integrals) and ell=10 (full-rank
% latent, matching the paper's TDL=10 tapped delay line).  One-step R2 at
% ell=10 is 0.996; the ell sweep was 5:0.41, 7:0.70, 9:0.92, 10:0.996.
% Protocol: 10-step init window (0.1 s) -> recursive prediction over H=190
% steps (1.9 s) with the recorded input; report per-horizon error
% distributions (mean and 99th percentile over all start points) and compare
% against the paper's hybrid-model numbers 9 cm/s and 0.12 rad/s.  The paper
% also reports MAE per axis (eq. 30): e = 1/3 sum|e_i|.
%
% Units: Vel m/s, pqr rad/s, Motors dimensionless.
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
here = fileparts(mfilename('fullpath'));
D = load(fullfile(here,'data','copyAY_pelican_dataset.mat'));

y_spd = double(D.y_speed_state); u_spd = double(D.u_speed_state);
seg_spd = double(D.segment_id_speed_state(:)');
train_mask = ismember(seg_spd,1:36);  val_mask = ismember(seg_spd,37:54);

% ---- scale with training statistics only (same contract as copyAY) ----
y_offset = mean(y_spd(:,train_mask),2); y_scale = std(y_spd(:,train_mask),0,2);
u_offset = mean(u_spd(:,train_mask),2); u_scale = std(u_spd(:,train_mask),0,2);
assert(all(y_scale>1e-10) && all(u_scale>1e-10),'zero training variance');
y_spd_tr = (y_spd(:,train_mask)-y_offset)./y_scale;
u_spd_tr = (u_spd(:,train_mask)-u_offset)./u_scale;

% ---- identify on speed state, same pipeline, soft preference on speed ----
% Preference: translational velocity and body rates are the task channels;
% motor speeds are internal (low weight). q=2 task directions.
opt = struct('weights',[1 1 1 1 1 1 0.05 0.05 0.05 0.05], ...
    'preference_strength',0.7,'reach_horizon',18,'Ru',eye(4));
ELL = 10;   % full-rank latent for speed state (ell sweep: 5->0.41 ... 10->0.996)
[E,learn_stats] = learn_segmented_preferred_output_directions( ...
    y_spd_tr,u_spd_tr,seg_spd(train_mask),2,opt);

[Ahat,Bhat,Phat,Rhat,Sigma_eps,fit_stats] = fit_segmented_anchored_varx( ...
    y_spd_tr,u_spd_tr,seg_spd(train_mask),E,ELL, ...
    struct('ridge',1e-8,'mu',0.10,'ntr_epsilon',10e-6));
model = struct('A',Ahat,'B',Bhat,'P',Phat,'R',Rhat, ...
    'Sigma_eps',Sigma_eps,'y_mean',fit_stats.y_mean,'u_mean',fit_stats.u_mean);

scales = struct('y_offset',y_offset,'y_scale',y_scale, ...
    'u_offset',u_offset,'u_scale',u_scale);

fprintf('speed-state fit: spectral radius %.6f, dual error %.2e\n', ...
    fit_stats.spectral_radius, fit_stats.dual_error);

% ---- recursive multi-step prediction on validation flights ----
H = 190; % 1.9 s at 100 Hz (paper long horizon; 0.4 s = 40 also reported)
out = evaluate_multistep_speed_state(y_spd(:,val_mask),u_spd(:,val_mask), ...
    seg_spd(val_mask),model,scales,H);
N = out.N;
fprintf('multistep: %d starts, H=%d (%.2f s)\n',N,out.H,out.horizon_s);

m_vel = mean(out.err_vel,1);        p99_vel = prctile(out.err_vel,99,1);
m_rb  = mean(out.err_rate,1);       p99_rb  = prctile(out.err_rate,99,1);
m_pvel = mean(out.per_vel,1);       m_prb  = mean(out.per_rate,1);
% paper norm (30): MAE per axis, deg/s for rates
ma_vel = mean(out.mae_vel,1);       ma_rb  = mean(out.mae_rate,1);
ma_pvel = mean(out.mae_pvel,1);     ma_prb = mean(out.mae_prate,1);

% one-step R2 on the validation set (standardized, full output)
y_val  = (y_spd(:,val_mask)-y_offset)./y_scale;
u_val  = (u_spd(:,val_mask)-u_offset)./u_scale;
zval = model.R'*(y_val - model.y_mean);
yp1  = model.y_mean + model.P*(model.A*zval + model.B*u_val);
one_r2 = 1 - sum((y_val-yp1).^2, 'all') / sum((y_val-model.y_mean).^2, 'all');
fprintf('speed-state one-step validation R2 (ell=%d): %.4f\n', ELL, one_r2);

h190 = H-1; h40 = 40; h18 = 18;  % 18 steps = 0.18 s ~ copyAY MPC horizon
fprintf('--- headline @1.9 s (paper hybrid: 9 cm/s, 0.12 rad/s, 99% conf) ---\n');
fprintf('speed  mean %.2f cm/s  p99 %.2f cm/s | persistence mean %.2f cm/s\n', ...
    100*m_vel(h190),100*p99_vel(h190),100*m_pvel(h190));
fprintf('bodyrate mean %.4f rad/s p99 %.4f rad/s | persistence mean %.4f rad/s\n', ...
    m_rb(h190),p99_rb(h190),m_prb(h190));
fprintf('--- @0.4 s ---\n');
fprintf('speed mean %.2f cm/s p99 %.2f cm/s | bodyrate mean %.4f rad/s p99 %.4f rad/s\n', ...
    100*m_vel(h40),100*p99_vel(h40),m_rb(h40),p99_rb(h40));
fprintf('--- @0.18 s (MPC horizon) ---\n');
fprintf('speed mean %.2f cm/s p99 %.2f cm/s | bodyrate mean %.4f rad/s p99 %.4f rad/s\n', ...
    100*m_vel(h18),100*p99_vel(h18),m_rb(h18),p99_rb(h18));
fprintf('--- MAE (paper eq.30) @1.9 s: speed %.2f cm/s (pers %.2f) | rate %.4f rad/s (pers %.4f) ---\n', ...
    100*ma_vel(h190),100*ma_pvel(h190),ma_rb(h190),ma_prb(h190));

% ---- figures ----
results_dir = fullfile(here,'results');
fig = figure('Position',[80 60 1200 800],'Color','w','Visible','off');
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
h_ax = 1:H;
nexttile;
plot(h_ax(1:H-1),100*m_vel,'g','LineWidth',1.5); hold on;
plot(h_ax(1:H-1),100*m_pvel,'m--','LineWidth',1.2);
yline(9,'r:','paper hybrid 9 cm/s'); grid on;
legend('speed err (cm/s)','persistence','Location','northwest');
xlabel('horizon step (10 ms)'); ylabel('cm/s');
title(sprintf('Speed-state multi-step prediction (1.9 s, %d starts)',N));
nexttile;
plot(h_ax(1:H-1),m_rb,'g','LineWidth',1.5); hold on;
plot(h_ax(1:H-1),m_prb,'m--','LineWidth',1.2);
yline(0.12,'r:','paper hybrid 0.12 rad/s'); grid on;
legend('body rate err (rad/s)','persistence','Location','northwest');
xlabel('horizon step (10 ms)'); ylabel('rad/s');
title('Body-rate multi-step prediction error vs horizon');
print(fig,fullfile(results_dir,'copyAY_pelican_speed_multistep_fig.png'),'-dpng','-r160'); close(fig);

save(fullfile(results_dir,'copyAY_pelican_speed_multistep_data.mat'), ...
    'out','m_vel','p99_vel','m_rb','p99_rb','m_pvel','m_prb', ...
    'learn_stats','fit_stats','-v7.3');
fid=fopen(fullfile(results_dir,'copyAY_pelican_speed_multistep_metrics.txt'),'w');
c=onCleanup(@() fclose(fid));
fprintf(fid,'copyAY Pelican SPEED-STATE multi-step benchmark (ICRA 2018 protocol)\n');
fprintf(fid,'state [Vel(3);pqr(3);Motors(4)]  latent_dim %d\n',ELL);
fprintf(fid,'horizon_steps %d\nhorizon_seconds %.2f\nstarts %d\n',out.H,out.horizon_s,N);
fprintf(fid,'one_step_validation_R2 %.6f\n',one_r2);
fprintf(fid,'speed_error_mean_1p9s_cmps %.4f\nspeed_error_p99_1p9s_cmps %.4f\n', ...
    100*m_vel(h190),100*p99_vel(h190));
fprintf(fid,'body_rate_error_mean_1p9s_radps %.6f\nbody_rate_error_p99_1p9s_radps %.6f\n', ...
    m_rb(h190),p99_rb(h190));
fprintf(fid,'speed_error_mean_0p4s_cmps %.4f\nspeed_error_p99_0p4s_cmps %.4f\n', ...
    100*m_vel(h40),100*p99_vel(h40));
fprintf(fid,'body_rate_error_mean_0p4s_radps %.6f\nbody_rate_error_p99_0p4s_radps %.6f\n', ...
    m_rb(h40),p99_rb(h40));
fprintf(fid,'speed_error_mean_0p18s_cmps %.4f\nspeed_error_p99_0p18s_cmps %.4f\n', ...
    100*m_vel(h18),100*p99_vel(h18));
fprintf(fid,'body_rate_error_mean_0p18s_radps %.6f\nbody_rate_error_p99_0p18s_radps %.6f\n', ...
    m_rb(h18),p99_rb(h18));
fprintf(fid,'mae_speed_1p9s_cmps %.4f\nmae_bodyrate_1p9s_radps %.6f\n', ...
    100*ma_vel(h190),ma_rb(h190));
fprintf(fid,'paper_hybrid_speed_99pct_cmps 9.0\npaper_hybrid_bodyrate_99pct_radps 0.12\n');
fprintf(fid,'paper_blackbox_lstm_bodyrate_mean_1p9s_degps 3.5\n');
fprintf(fid,'speed_persistence_mean_1p9s_cmps %.4f\nbodyrate_persistence_mean_1p9s_radps %.6f\n', ...
    100*m_pvel(h190),m_prb(h190));
fprintf(fid,'fit_spectral_radius %.10f\n',fit_stats.spectral_radius);
fprintf(fid,'model copyAY pipeline on speed state (ell=%d latent, first-order, standardized)\n',ELL);
