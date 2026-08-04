%% copyBD_pelican_3d_fixed_anchor
% Three prescribed standardized position axes on real Waterloo Pelican data.
% Flights 1:36 identify the model; flights 37:54 remain held out.
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','copyBA_pelican_position_task'));
data_dir = fullfile(here,'..','copyAY_pelican_soft_preference','data');
D = load(fullfile(data_dir,'copyAY_pelican_dataset.mat'));
assert(isequal(size(D.y,1),10) && isequal(size(D.u,1),4), ...
    'Expected 10 outputs and 4 commanded inputs.');

segment_id = double(D.segment_id(:)');
train_mask = ismember(segment_id,1:36);
validation_mask = ismember(segment_id,37:54);
y_train_raw = double(D.y(:,train_mask));
u_train_raw = double(D.u(:,train_mask));
y_validation_raw = double(D.y(:,validation_mask));
u_validation_raw = double(D.u(:,validation_mask));
train_segment_id = segment_id(train_mask);
validation_segment_id = segment_id(validation_mask);

% Training statistics only.
y_offset = mean(y_train_raw,2); y_scale = std(y_train_raw,0,2);
u_offset = mean(u_train_raw,2); u_scale = std(u_train_raw,0,2);
y_train = (y_train_raw-y_offset)./y_scale;
u_train = (u_train_raw-u_offset)./u_scale;
y_validation = (y_validation_raw-y_offset)./y_scale;
u_validation = (u_validation_raw-u_offset)./u_scale;

p=10; m=4; q=3; ell=5; metric_mu=0.10; ntr_epsilon=10e-6;
Eanchor=zeros(p,q); Eanchor(1:3,1:3)=eye(3);
assert(norm(Eanchor'*Eanchor-eye(q),'fro')<eps, ...
    'The prescribed x/y/z anchor must be orthonormal.');
[Ahat,Bhat,Phat,Rhat,Sigma_eps,fitstats]=fit_segmented_anchored_varx( ...
    y_train,u_train,train_segment_id,Eanchor,ell,struct('ridge',1e-8, ...
    'mu',metric_mu,'ntr_epsilon',ntr_epsilon));
model=struct('A',Ahat,'B',Bhat,'P',Phat,'R',Rhat, ...
    'Sigma_eps',Sigma_eps,'y_mean',fitstats.y_mean,'u_mean',fitstats.u_mean);
validation=evaluate_segmented_prediction(y_validation,u_validation, ...
    validation_segment_id,Eanchor,model);

train_task=Eanchor'*y_train;
task_bound=max(3.0,max(abs(train_task),[],'all'));
results_dir=fullfile(here,'results');
if ~exist(results_dir,'dir'), mkdir(results_dir); end
metrics_path=fullfile(results_dir,'copyBD_pelican_3d_fixed_anchor_metrics.txt');
fid=fopen(metrics_path,'w'); cleaner=onCleanup(@() fclose(fid));
fprintf(fid,'copyBD three-dimensional fixed-anchor Pelican identification\n');
fprintf(fid,'training_source real_waterloo_pelican_flights_1_36\n');
fprintf(fid,'validation_source held_out_flights_37_54\n');
fprintf(fid,'controlled_directions prescribed_standardized_x_y_z_position\n');
fprintf(fid,'latent_dimension %d\nfree_dimension %d\n',ell,ell-q);
fprintf(fid,'task_bound_training_only %.12f\n',task_bound);
fprintf(fid,'dual_error %.12e\nanchor_preservation_error %.12e\n', ...
    fitstats.dual_error,fitstats.anchor_preservation_error);
fprintf(fid,'spectral_radius %.12f\n',fitstats.spectral_radius);
fprintf(fid,'task_RMSE %s\ntask_persistence_RMSE %s\n', ...
    mat2str(validation.task_rmse',10),mat2str(validation.task_persistence_rmse',10));
fprintf(fid,'task_R2 %s\n',mat2str(validation.task_r2',10));
fprintf(fid,'scope offline_identification_before_3d_model_in_loop_control\n');
clear cleaner;

save(fullfile(results_dir,'copyBD_pelican_3d_fixed_anchor_data.mat'), ...
    'Eanchor','metric_mu','ntr_epsilon','task_bound','fitstats','model', ...
    'validation','y_offset','y_scale','u_offset','u_scale','-v7.3');

fig=figure('Position',[80 80 1800 1050],'Color','w','Visible','off');
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
tt=1:validation.transition_count; labels={'x','y','z'};
for j=1:3
    ax=nexttile;
    plot(tt,validation.task_actual(j,:),'k','LineWidth',1.0); hold on;
    plot(tt,validation.task_predicted(j,:),'Color',[0.03 0.49 0.55],'LineWidth',1.0);
    grid on; ylabel(sprintf('standardized %s',labels{j}));
    title(sprintf('Held-out one-step prediction: fixed %s-position direction',labels{j}));
    legend('actual','prediction','Location','best');
end
xlabel(ax,'validation transition index');
title(tl,sprintf('copyBD 3D fixed-anchor identification | RMSE=%s | R^2=%s', ...
    mat2str(validation.task_rmse',3),mat2str(validation.task_r2',3)));
print(fig,fullfile(results_dir,'copyBD_pelican_3d_fixed_anchor_prediction.png'),'-dpng','-r160');
close(fig);
fprintf('copyBD identification complete: RMSE=%s, persistence=%s, R2=%s\n', ...
    mat2str(validation.task_rmse',4),mat2str(validation.task_persistence_rmse',4), ...
    mat2str(validation.task_r2',4));
