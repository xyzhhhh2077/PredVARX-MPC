%% copyAX_cdr_soft_preference
% Apply the copyAU soft-preference CRTE construction to ControlGym CDR data.
% CDR: convection-diffusion-reaction PDE, linear analytical dynamics
%   x_{k+1} = A x_k + B2 u_k + w_k,  y_k = C x_k + v_k,
% 200-state field observed at 30 sensors, 8 localized actuators.
% Alignment (standard causality): u(:,k) drives y(:,k+1).
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
here = fileparts(mfilename('fullpath'));
D = load(fullfile(here,'data','copyAX_cdr_dataset.mat'));

assert(isequal(size(D.y,1),30) && isequal(size(D.u,1),8), ...
    'Expected 30 CDR sensor outputs and 8 manipulated inputs.');
assert(size(D.y,2)==size(D.u,2) && size(D.y,2)==numel(D.segment_id), ...
    'CDR sample counts do not match.');
assert(all(isfinite(D.y),'all') && all(isfinite(D.u),'all'), ...
    'CDR dataset contains non-finite values.');

segment_id = double(D.segment_id(:)');
train_mask = ismember(segment_id,1:6);
validation_mask = ismember(segment_id,7:8);
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
    'A CDR channel has zero training variance.');
y_train = (y_train_raw-y_offset)./y_scale;
u_train = (u_train_raw-u_offset)./u_scale;
y_validation = (y_validation_raw-y_offset)./y_scale;
u_validation = (u_validation_raw-u_offset)./u_scale;

p=30; m=8; q=2; ell=5; reach_horizon=18;
preference_strength=0.78; metric_mu=0.10; ntr_epsilon=10e-6;
% Soft spatial preference: sensors in the central bulk of the domain (grid
% indices 40..160 of 200, i.e. sensor idx 7..26) weight 1.0; the 10 edge
% sensors weight 0.85. All outputs remain represented; no hard task lock.
weights=[0.85*ones(7,1);ones(20,1);0.85*ones(3,1)];
assert(numel(weights)==p,'weights must cover all 30 sensors.');
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
metrics_path=fullfile(results_dir,'copyAX_cdr_soft_preference_metrics.txt');
fid=fopen(metrics_path,'w'); cleaner=onCleanup(@() fclose(fid));
fprintf(fid,'copyAX ControlGym CDR soft preference offline validation\n');
fprintf(fid,'environment convection_diffusion_reaction (linear analytical)\n');
fprintf(fid,'native_measurements %d\nmanipulated_inputs %d\nlatent_dimension %d\n',p,m,ell);
fprintf(fid,'training_segments %s\nvalidation_segments %s\n', ...
    mat2str(unique(train_segment_id)),mat2str(unique(validation_segment_id)));
fprintf(fid,'training_samples %d\nvalidation_samples %d\n', ...
    numel(train_segment_id),numel(validation_segment_id));
fprintf(fid,'training_transitions %d\nvalidation_transitions %d\n', ...
    fitstats.transition_count,validation.transition_count);
fprintf(fid,'preference_weights central_sensors=1.0 edge_sensors=0.85\n');
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
fprintf(fid,'scope active_cdr_simulator_offline_one_step_prediction_only\n');
fprintf(fid,'closed_loop_mpc_run 0\nexogenous_disturbance_modeled 0\n');
clear cleaner;

save(fullfile(results_dir,'copyAX_cdr_soft_preference_data.mat'), ...
    'Etask','weights','preference_strength','metric_mu','ntr_epsilon', ...
    'dirstats','fitstats','model','validation','y_offset','y_scale', ...
    'u_offset','u_scale','-v7.3');

fig=figure('Position',[80 60 1900 1250],'Color','w','Visible','off');
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
tt=1:validation.transition_count;
ax1=nexttile; plot(tt,validation.task_actual(1,:),'k', ...
    tt,validation.task_predicted(1,:),'b','LineWidth',1.0); grid on;
ylabel('task 1'); title('Held-out CDR validation segments: one-step task prediction');
legend('actual','copyAX prediction','Location','best');
ax2=nexttile; plot(tt,validation.task_actual(2,:),'k', ...
    tt,validation.task_predicted(2,:),'r','LineWidth',1.0); grid on;
ylabel('task 2'); title('Held-out CDR validation segments: one-step task prediction');
legend('actual','copyAX prediction','Location','best');
ax3=nexttile; bar(1:p,dirstats.contribution,'FaceColor',[0.20 0.45 0.72]); hold on;
plot(1:p,weights/max(weights)*max(dirstats.contribution),'r--','LineWidth',1.2);
grid on; xlim([0 p+1]); xlabel('CDR sensor index (grid 6*idx, 0..174)'); ylabel('direction contribution');
title('30 sensors: 20 central weight 1.0 | 10 edge weight 0.85');
legend('learned contribution','scaled soft preference','Location','best');
linkaxes([ax1 ax2],'x'); xlabel(ax2,'validation transition index');
title(tl,sprintf('copyAX | ControlGym CDR 30 sensors -> 5 latent variables | task RMSE=%s', ...
    mat2str(validation.task_rmse',3)));
print(fig,fullfile(results_dir,'copyAX_cdr_soft_preference_fig.png'),'-dpng','-r160'); close(fig);
fprintf('copyAX CDR complete: task RMSE=%s, persistence=%s, R2=%s\n', ...
    mat2str(validation.task_rmse',4), ...
    mat2str(validation.task_persistence_rmse',4),mat2str(validation.task_r2',4));
