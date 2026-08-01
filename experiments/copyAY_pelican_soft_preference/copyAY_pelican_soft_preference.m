%% copyAY_pelican_soft_preference
% Apply the copyAU soft-preference CRTE construction to the Waterloo Pelican
% quadrotor flight dataset (real indoor flight data, Vicon motion capture).
%
%   y = [Pos(3); Euler(3); Motors(4)]      (10 raw measurements, 100 Hz)
%   u = Motors_CMD(4)                       (commanded motor speeds)
%
% Alignment (standard causality): u(:,k) drives y(:,k+1). Vel and pqr are
% excluded (derived linearly from Pos/Euler -> singular C_y for C_y^{-1/2}).
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
here = fileparts(mfilename('fullpath'));
D = load(fullfile(here,'data','copyAY_pelican_dataset.mat'));

assert(isequal(size(D.y,1),10) && isequal(size(D.u,1),4), ...
    'Expected 10 Pelican outputs (Pos+Euler+Motors) and 4 commanded inputs.');
assert(size(D.y,2)==size(D.u,2) && size(D.y,2)==numel(D.segment_id), ...
    'Pelican sample counts do not match.');
assert(all(isfinite(D.y),'all') && all(isfinite(D.u),'all'), ...
    'Pelican dataset contains non-finite values.');

segment_id = double(D.segment_id(:)');
train_mask = ismember(segment_id,1:36);
validation_mask = ismember(segment_id,37:54);
y_train_raw = double(D.y(:,train_mask));
u_train_raw = double(D.u(:,train_mask));
y_validation_raw = double(D.y(:,validation_mask));
u_validation_raw = double(D.u(:,validation_mask));
train_segment_id = segment_id(train_mask);
validation_segment_id = segment_id(validation_mask);

% Scale with training statistics only. No validation sample sets model parameters.
y_offset = mean(y_train_raw,2); y_scale = std(y_train_raw,0,2);
u_offset = mean(u_train_raw,2); u_scale = std(u_train_raw,0,2);
assert(all(y_scale>1e-10) && all(u_scale>1e-10), ...
    'A Pelican channel has zero training variance.');
y_train = (y_train_raw-y_offset)./y_scale;
u_train = (u_train_raw-u_offset)./u_scale;
y_validation = (y_validation_raw-y_offset)./y_scale;
u_validation = (u_validation_raw-u_offset)./u_scale;

p=10; m=4; q=2; ell=5; reach_horizon=18;
preference_strength=0.78; metric_mu=0.10; ntr_epsilon=10e-6;
% Soft preference over the 10 raw channels: flight state (position, Euler
% attitude) is the control task; motor speeds are internal actuator state.
%   Pos   (3): 1.00  -- position tracking is the outer-loop task
%   Euler (3): 0.90  -- attitude serves position (inner loop)
%   Motors(4): 0.70  -- actuator internals, still represented (no hard lock)
weights=[1.0*ones(3,1);0.9*ones(3,1);0.7*ones(4,1)];
assert(numel(weights)==p,'weights must cover all 10 outputs.');
[Etask,dirstats]=learn_segmented_preferred_output_directions( ...
    y_train,u_train,train_segment_id,q,struct('weights',weights, ...
    'preference_strength',preference_strength,'reach_horizon',reach_horizon, ...
    'Ru',eye(m),'ridge',1e-8));
[Ahat,Bhat,Phat,Rhat,Sigma_eps,fitstats]=fit_segmented_anchored_varx( ...
    y_train,u_train,train_segment_id,Etask,ell,struct('ridge',1e-8, ...
    'mu',metric_mu,'ntr_epsilon',ntr_epsilon));

model=struct('A',Ahat,'B',Bhat,'P',Phat,'R',Rhat, ...
    'Sigma_eps',Sigma_eps,'y_mean',fitstats.y_mean,'u_mean',fitstats.u_mean);
validation=evaluate_segmented_prediction(y_validation,u_validation, ...
    validation_segment_id,Etask,model);

results_dir=fullfile(here,'results');
if ~exist(results_dir,'dir'), mkdir(results_dir); end
metrics_path=fullfile(results_dir,'copyAY_pelican_soft_preference_metrics.txt');
fid=fopen(metrics_path,'w'); cleaner=onCleanup(@() fclose(fid));
fprintf(fid,'copyAY Waterloo Pelican quadrotor soft preference offline validation\n');
fprintf(fid,'environment real_flight_data_vicon_100hz\n');
fprintf(fid,'native_measurements %d\nmanipulated_inputs %d\nlatent_dimension %d\n',p,m,ell);
fprintf(fid,'training_segments %s\nvalidation_segments %s\n', ...
    mat2str(unique(train_segment_id)),mat2str(unique(validation_segment_id)));
fprintf(fid,'training_samples %d\nvalidation_samples %d\n', ...
    numel(train_segment_id),numel(validation_segment_id));
fprintf(fid,'training_transitions %d\nvalidation_transitions %d\n', ...
    fitstats.transition_count,validation.transition_count);
fprintf(fid,'preference_weights pos=1.0 euler=0.9 motors=0.7\n');
fprintf(fid,'preference_strength %.8f\nselected_mu %.12f\n', ...
    preference_strength,fitstats.selected_mu);
fprintf(fid,'ntr_epsilon %.12e\ndual_error %.12e\n', ...
    fitstats.ntr_epsilon,fitstats.dual_error);
fprintf(fid,'spectral_radius %.12f\npreference_capture %.12f\n', ...
    fitstats.spectral_radius,dirstats.preference_capture);
fprintf(fid,'latent_RMSE %s\nlatent_persistence_RMSE %s\n', ...
    mat2str(validation.latent_rmse',10),mat2str(validation.latent_persistence_rmse',10));
fprintf(fid,'task_RMSE %s\ntask_persistence_RMSE %s\n', ...
    mat2str(validation.task_rmse',10),mat2str(validation.task_persistence_rmse',10));
fprintf(fid,'task_R2 %s\nfull_output_standardized_RMSE %.12f\n', ...
    mat2str(validation.task_r2',10),validation.full_output_rmse);
fprintf(fid,'scope offline_one_step_prediction_only\n');
fprintf(fid,'closed_loop_mpc_run 0\nexogenous_disturbance_modeled 0\n');
clear cleaner;

save(fullfile(results_dir,'copyAY_pelican_soft_preference_data.mat'), ...
    'Etask','weights','preference_strength','metric_mu','ntr_epsilon', ...
    'dirstats','fitstats','model','validation','y_offset','y_scale', ...
    'u_offset','u_scale','-v7.3');

fig=figure('Position',[80 60 1900 1250],'Color','w','Visible','off');
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
tt=1:validation.transition_count;
ax1=nexttile; plot(tt,validation.task_actual(1,:),'k', ...
    tt,validation.task_predicted(1,:),'b','LineWidth',1.0); grid on;
ylabel('task 1'); title('Held-out Pelican flights (37-54): one-step task prediction');
legend('actual','copyAY prediction','Location','best');
ax2=nexttile; plot(tt,validation.task_actual(2,:),'k', ...
    tt,validation.task_predicted(2,:),'r','LineWidth',1.0); grid on;
ylabel('task 2'); title('Held-out Pelican flights (37-54): one-step task prediction');
legend('actual','copyAY prediction','Location','best');
ax3=nexttile; bar(1:p,dirstats.contribution,'FaceColor',[0.20 0.45 0.72]); hold on;
plot(1:p,weights/max(weights)*max(dirstats.contribution),'r--','LineWidth',1.2);
grid on; xlim([0 p+1]); xlabel('channel index (1-3 pos, 4-6 euler, 7-10 motors)');
ylabel('direction contribution');
title('10 channels: pos 1.00 | euler 0.90 | motors 0.70');
legend('learned contribution','scaled soft preference','Location','best');
linkaxes([ax1 ax2],'x'); xlabel(ax2,'validation transition index');
title(tl,sprintf('copyAY | Pelican 10 channels -> 5 latent variables | task RMSE=%s', ...
    mat2str(validation.task_rmse',3)));
print(fig,fullfile(results_dir,'copyAY_pelican_soft_preference_fig.png'),'-dpng','-r160'); close(fig);
fprintf('copyAY Pelican complete: task RMSE=%s, persistence=%s, R2=%s\n', ...
    mat2str(validation.task_rmse',4), ...
    mat2str(validation.task_persistence_rmse',4),mat2str(validation.task_r2',4));
