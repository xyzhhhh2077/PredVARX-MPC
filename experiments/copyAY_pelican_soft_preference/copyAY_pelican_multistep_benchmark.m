%% copyAY_pelican_multistep_benchmark
% Multi-step prediction benchmark on the Waterloo Pelican quadrotor dataset,
% following the protocol of Mohajerin & Waslander, ICRA 2018 (arXiv
% 1806.00526): recursively iterate the learned model over H steps with the
% recorded input sequence; report per-horizon error distributions (mean and
% 99th percentile over all start points) and the paper's headline numbers
% (hybrid model: speed within 9 cm/s and body rate within 0.12 rad/s over
% 1.9 s with 99% confidence on the test dataset).
%
% Units: Pos m, Euler rad, Motors dimensionless; speed = 100*diff(Pos) m/s,
% body rate = T(attitude)*euler_rate rad/s (dataset stores pqr this way).
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
here = fileparts(mfilename('fullpath'));
D = load(fullfile(here,'data','copyAY_pelican_dataset.mat'));
R = load(fullfile(here,'results','copyAY_pelican_soft_preference_data.mat'));

segment_id = double(D.segment_id(:)');
validation_mask = ismember(segment_id,37:54);
y_val_raw = double(D.y(:,validation_mask));
u_val_raw = double(D.u(:,validation_mask));
seg_val = segment_id(validation_mask);

scales = struct('y_offset',R.y_offset,'y_scale',R.y_scale, ...
    'u_offset',R.u_offset,'u_scale',R.u_scale);

H = 190; % 1.9 s at 100 Hz (paper's long horizon; 40 = 0.4 s also reported)
out = evaluate_multistep_prediction(y_val_raw,u_val_raw,seg_val,R.model,scales,H);
N = out.N;
fprintf('multistep: %d starts, H=%d (%.2f s)\n',N,out.H,out.horizon_s);

% ---- per-horizon mean and 99th percentile over starts ----
m_vel = mean(out.err_vel,1);        p99_vel = prctile(out.err_vel,99,1);
m_rb  = mean(out.err_rate_body,1);  p99_rb  = prctile(out.err_rate_body,99,1);
m_pos = mean(out.err_pos,1);        m_att  = mean(out.err_att,1);
m_pvel = mean(out.per_vel,1);       m_ppos = mean(out.per_pos,1);
m_patt = mean(out.per_att,1);       m_prb  = mean(out.per_rate_body,1);

h190 = H-1; h40 = 40;
fprintf('--- headline @1.9 s (paper hybrid: 9 cm/s speed, 0.12 rad/s rate, 99%% conf) ---\n');
fprintf('speed  mean %.2f cm/s  p99 %.2f cm/s | persistence mean %.2f cm/s\n', ...
    100*m_vel(h190),100*p99_vel(h190),100*m_pvel(h190));
fprintf('bodyrate mean %.4f rad/s p99 %.4f rad/s | persistence mean %.4f rad/s\n', ...
    m_rb(h190),p99_rb(h190),m_prb(h190));
fprintf('pos    mean %.2f cm | persistence %.2f cm\n',100*m_pos(H),100*m_ppos(H));
fprintf('att    mean %.2f deg | persistence %.2f deg\n',180/pi*m_att(H),180/pi*m_patt(H));
fprintf('--- @0.4 s ---\n');
fprintf('speed mean %.2f cm/s p99 %.2f cm/s | bodyrate mean %.4f rad/s\n', ...
    100*m_vel(h40),100*p99_vel(h40),m_rb(h40));
fprintf('pos mean %.2f cm | att mean %.2f deg\n',100*m_pos(40),180/pi*m_att(40));

% ---- figures ----
results_dir = fullfile(here,'results');
fig = figure('Position',[80 60 1900 1150],'Color','w','Visible','off');
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
h_ax = 1:H;
nexttile;
plot(h_ax,100*m_pos,'b','LineWidth',1.4); hold on;
plot(h_ax,100*m_ppos,'k--','LineWidth',1.2);
plot(h_ax(1:H-1),100*m_vel,'g','LineWidth',1.4);
plot(h_ax(1:H-1),100*m_pvel,'m--','LineWidth',1.2);
yline(9,'r:','paper hybrid 9 cm/s'); grid on;
legend('pos err (cm)','persist pos','speed err (cm/s)','persist speed','Location','northwest');
xlabel('horizon step (10 ms)'); ylabel('cm');
title(sprintf('Multi-step prediction error vs horizon (1.9 s, %d starts)',N));
nexttile;
plot(h_ax,180/pi*m_att,'b','LineWidth',1.4); hold on;
plot(h_ax,180/pi*m_patt,'k--','LineWidth',1.2);
plot(h_ax(1:H-1),m_rb,'g','LineWidth',1.4);
plot(h_ax(1:H-1),m_prb,'m--','LineWidth',1.2);
yline(0.12,'r:','paper hybrid 0.12 rad/s'); grid on;
legend('att err (deg)','persist att','body rate err (rad/s)','persist rate','Location','northwest');
xlabel('horizon step (10 ms)'); ylabel('deg / rad/s');
title('Attitude and body-rate error vs horizon');
nexttile;
g = m_pos./max(m_pos(1),eps);
plot(h_ax,g,'b','LineWidth',1.4); grid on; yline(1,'k--');
xlabel('horizon step (10 ms)'); ylabel('pos err / 1-step pos err');
title(sprintf('Position error growth (spectral radius %.6f)',R.fitstats.spectral_radius));
print(fig,fullfile(results_dir,'copyAY_pelican_multistep_fig.png'),'-dpng','-r160'); close(fig);

save(fullfile(results_dir,'copyAY_pelican_multistep_data.mat'),'out','m_vel','p99_vel', ...
    'm_rb','p99_rb','m_pos','m_att','m_pvel','m_ppos','m_patt','m_prb','-v7.3');
fid=fopen(fullfile(results_dir,'copyAY_pelican_multistep_metrics.txt'),'w');
c=onCleanup(@() fclose(fid));
fprintf(fid,'copyAY Pelican multi-step benchmark (ICRA 2018 protocol)\n');
fprintf(fid,'horizon_steps %d\nhorizon_seconds %.2f\nstarts %d\n',out.H,out.horizon_s,N);
fprintf(fid,'speed_error_mean_1p9s_cmps %.4f\nspeed_error_p99_1p9s_cmps %.4f\n', ...
    100*m_vel(h190),100*p99_vel(h190));
fprintf(fid,'body_rate_error_mean_1p9s_radps %.6f\nbody_rate_error_p99_1p9s_radps %.6f\n', ...
    m_rb(h190),p99_rb(h190));
fprintf(fid,'pos_error_mean_1p9s_cm %.4f\npos_persistence_mean_1p9s_cm %.4f\n', ...
    100*m_pos(H),100*m_ppos(H));
fprintf(fid,'att_error_mean_1p9s_rad %.6f\natt_persistence_mean_1p9s_rad %.6f\n', ...
    m_att(H),m_patt(H));
fprintf(fid,'speed_error_mean_0p4s_cmps %.4f\nspeed_error_p99_0p4s_cmps %.4f\n', ...
    100*m_vel(h40),100*p99_vel(h40));
fprintf(fid,'paper_hybrid_speed_99pct_cmps 9.0\npaper_hybrid_bodyrate_99pct_radps 0.12\n');
fprintf(fid,'spectral_radius %.10f\n',R.fitstats.spectral_radius);
fprintf(fid,'model copyAY one-step latent VARX (5 latent, first-order, standardized)\n');
