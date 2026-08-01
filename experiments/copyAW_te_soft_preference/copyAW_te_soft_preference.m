%% copyAW_te_soft_preference
% Apply the copyAU soft-preference CRTE construction to native 41-output TEP data.
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
here=fileparts(mfilename('fullpath'));
D=load(fullfile(here,'data','tep_fault_free_runs_001_005.mat'));
assert(strcmp(strtrim(D.source_doi),'10.7910/DVN/6C3JR1'),'TEP DOI provenance mismatch.');
assert(strcmp(strtrim(D.source_file_id),'3031241'),'TEP Dataverse file ID mismatch.');
assert(strcmp(strtrim(D.source_md5),'ec126484534331f85001d8c4ebce6d17'), ...
    'TEP source MD5 provenance mismatch.');
assert(isequal(size(D.y),[41 2500]) && isequal(size(D.u),[11 2500]), ...
    'Unexpected TEP snapshot dimensions.');
assert(all(D.fault_number==0),'copyAW accepts only fault-free runs.');

train_mask=ismember(D.run_id,1:3); validation_mask=ismember(D.run_id,4:5);
y_train_raw=D.y(:,train_mask); u_train_raw=D.u(:,train_mask);
y_validation_raw=D.y(:,validation_mask); u_validation_raw=D.u(:,validation_mask);
train_run_id=D.run_id(train_mask); validation_run_id=D.run_id(validation_mask);

% Scale with training statistics only. Validation data never sets model parameters.
y_offset=mean(y_train_raw,2); y_scale=std(y_train_raw,0,2);
u_offset=mean(u_train_raw,2); u_scale=std(u_train_raw,0,2);
assert(all(y_scale>0) && all(u_scale>0),'A TEP channel has zero training variance.');
y_train=(y_train_raw-y_offset)./y_scale;
u_train=(u_train_raw-u_offset)./u_scale;
y_validation=(y_validation_raw-y_offset)./y_scale;
u_validation=(u_validation_raw-u_offset)./u_scale;

p=41; m=11; q=2; ell=5; reach_horizon=18;
preference_strength=0.78; metric_mu=0.10; ntr_epsilon=10e-6;
weights=0.10*ones(p,1);
preferred_measurements=[7 8 9 11 12 13 15 16 18];
weights(preferred_measurements)=1.0;
[Etask,dirstats]=learn_segmented_preferred_output_directions( ...
    y_train,u_train,train_run_id,q,struct('weights',weights, ...
    'preference_strength',preference_strength,'reach_horizon',reach_horizon, ...
    'Ru',eye(m),'ridge',1e-8));
[Ahat,Bhat,Phat,Rhat,Sigma_eps,fitstats]=fit_segmented_anchored_varx( ...
    y_train,u_train,train_run_id,Etask,ell,struct('ridge',1e-8, ...
    'mu',metric_mu,'ntr_epsilon',ntr_epsilon));

model=struct('A',Ahat,'B',Bhat,'P',Phat,'R',Rhat, ...
    'Sigma_eps',Sigma_eps,'y_mean',fitstats.y_mean,'u_mean',fitstats.u_mean);
validation=evaluate_segmented_prediction(y_validation,u_validation, ...
    validation_run_id,Etask,model);

results_dir=fullfile(here,'results');
if ~exist(results_dir,'dir'), mkdir(results_dir); end
metrics_path=fullfile(results_dir,'copyAW_te_soft_preference_metrics.txt');
fid=fopen(metrics_path,'w'); cleaner=onCleanup(@() fclose(fid));
fprintf(fid,'copyAW TEP soft preference offline validation\n');
fprintf(fid,'source_doi %s\nsource_file_id %s\nsource_md5 %s\n', ...
    strtrim(D.source_doi),strtrim(D.source_file_id),strtrim(D.source_md5));
fprintf(fid,'native_measurements %d\nmanipulated_inputs %d\nlatent_dimension %d\n',p,m,ell);
fprintf(fid,'training_runs %s\nvalidation_runs %s\ntraining_samples %d\nvalidation_samples %d\n', ...
    mat2str(unique(train_run_id)),mat2str(unique(validation_run_id)), ...
    numel(train_run_id),numel(validation_run_id));
fprintf(fid,'training_transitions %d\nvalidation_transitions %d\n', ...
    fitstats.transition_count,validation.transition_count);
fprintf(fid,'preferred_measurements %s\npreference_strength %.8f\n', ...
    mat2str(preferred_measurements),preference_strength);
fprintf(fid,'selected_mu %.12f\nntr_epsilon %.12e\ndual_error %.12e\n', ...
    fitstats.selected_mu,fitstats.ntr_epsilon,fitstats.dual_error);
fprintf(fid,'spectral_radius %.12f\npreference_capture %.12f\n', ...
    fitstats.spectral_radius,dirstats.preference_capture);
fprintf(fid,'latent_RMSE %s\nlatent_persistence_RMSE %s\n', ...
    mat2str(validation.latent_rmse',10),mat2str(validation.latent_persistence_rmse',10));
fprintf(fid,'task_RMSE %s\ntask_persistence_RMSE %s\n', ...
    mat2str(validation.task_rmse',10),mat2str(validation.task_persistence_rmse',10));
fprintf(fid,'task_R2 %s\nfull_output_standardized_RMSE %.12f\n', ...
    mat2str(validation.task_r2',10),validation.full_output_rmse);
fprintf(fid,'scope offline_one_step_prediction_only\nclosed_loop_costep_run 0\n');
clear cleaner;

save(fullfile(results_dir,'copyAW_te_soft_preference_data.mat'), ...
    'Etask','weights','preferred_measurements','preference_strength', ...
    'metric_mu','ntr_epsilon','dirstats','fitstats','model','validation', ...
    'y_offset','y_scale','u_offset','u_scale','-v7.3');

fig=figure('Position',[80 60 1900 1250],'Color','w','Visible','off');
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
idx=validation.run_id==4; tt=1:sum(idx);
ax1=nexttile; plot(tt,validation.task_actual(1,idx),'k', ...
    tt,validation.task_predicted(1,idx),'b'); grid on;
ylabel('task 1'); title('Held-out TEP run 4: one-step task prediction');
legend('actual','copyAW prediction','Location','eastoutside');
ax2=nexttile; plot(tt,validation.task_actual(2,idx),'k', ...
    tt,validation.task_predicted(2,idx),'r'); grid on;
ylabel('task 2'); title('Held-out TEP run 4: one-step task prediction');
legend('actual','copyAW prediction','Location','eastoutside');
ax3=nexttile; bar(1:p,dirstats.contribution,'FaceColor',[0.20 0.45 0.72]); hold on;
plot(1:p,weights/max(weights)*max(dirstats.contribution),'r--','LineWidth',1.2);
grid on; xlim([0 p+1]); xlabel('TEP measurement index'); ylabel('direction contribution');
title('Learned output-direction contribution and physical soft preference');
legend('learned contribution','scaled preference','Location','eastoutside');
linkaxes([ax1 ax2],'x');
title(tl,sprintf('copyAW | TEP 41 outputs -> 5 latent variables | held-out runs 4-5 | task RMSE=%s', ...
    mat2str(validation.task_rmse',3)));
print(fig,fullfile(results_dir,'copyAW_te_soft_preference_fig.png'),'-dpng','-r160'); close(fig);
fprintf('copyAW complete: task RMSE=%s, persistence=%s, R2=%s\n', ...
    mat2str(validation.task_rmse',4), ...
    mat2str(validation.task_persistence_rmse',4),mat2str(validation.task_r2',4));