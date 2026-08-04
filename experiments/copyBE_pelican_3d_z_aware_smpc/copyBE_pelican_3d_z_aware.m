%% copyBE_pelican_3d_z_aware
% Fixed x/y/z anchors with second-order latent VARX dynamics.
% The extra latent lag represents position-rate memory needed by vertical motion.
clear; clc; close all;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','copyBA_pelican_position_task'));
data_dir = fullfile(here,'..','copyAY_pelican_soft_preference','data');
D = load(fullfile(data_dir,'copyAY_pelican_dataset.mat'));
assert(isequal(size(D.y,1),10) && isequal(size(D.u,1),4), ...
    'Expected 10 outputs and 4 commanded motor inputs.');

segment_id = double(D.segment_id(:)');
train_mask = ismember(segment_id,1:36);
validation_mask = ismember(segment_id,37:54);
y_train_raw = double(D.y(:,train_mask));
u_train_raw = double(D.u(:,train_mask));
y_validation_raw = double(D.y(:,validation_mask));
u_validation_raw = double(D.u(:,validation_mask));
train_segment_id = segment_id(train_mask);
validation_segment_id = segment_id(validation_mask);

y_offset = mean(y_train_raw,2); y_scale = std(y_train_raw,0,2);
u_offset = mean(u_train_raw,2); u_scale = std(u_train_raw,0,2);
y_train = (y_train_raw-y_offset)./y_scale;
u_train = (u_train_raw-u_offset)./u_scale;
y_validation = (y_validation_raw-y_offset)./y_scale;
u_validation = (u_validation_raw-u_offset)./u_scale;

p=10; q=3; ell=5; ridge=1e-8; metric_mu=0.10; ntr_epsilon=10e-6;
Eanchor=zeros(p,q); Eanchor(1:3,1:3)=eye(3);
[~,~,Phat,Rhat,~,anchor_stats]=fit_segmented_anchored_varx( ...
    y_train,u_train,train_segment_id,Eanchor,ell,struct('ridge',ridge, ...
    'mu',metric_mu,'ntr_epsilon',ntr_epsilon));

% Refit dynamics with one additional latent-state lag. For transition t -> t+1,
% regress z(t+1) on z(t), z(t-1), and u(t), never crossing flight boundaries.
y_mean=mean(y_train,2); u_mean=mean(u_train,2);
zc=Rhat'*(y_train-y_mean); uc=u_train-u_mean;
valid=find(train_segment_id(1:end-2)==train_segment_id(2:end-1) & ...
    train_segment_id(2:end-1)==train_segment_id(3:end))+1;
Phi=[zc(:,valid);zc(:,valid-1);uc(:,valid)]; target=zc(:,valid+1);
scale=max(trace(Phi*Phi'/size(Phi,2))/size(Phi,1),1e-12);
lambda=ridge*max(scale,1)+1e-12;
Theta=(Phi*Phi'+lambda*eye(size(Phi,1)))\(Phi*target');
A1=Theta(1:ell,:)'; A2=Theta(ell+1:2*ell,:)'; Bhat=Theta(2*ell+1:end,:)';
Eps=target-A1*zc(:,valid)-A2*zc(:,valid-1)-Bhat*uc(:,valid);
Sigma_eps=Eps*Eps'/size(Eps,2); Sigma_eps=(Sigma_eps+Sigma_eps')/2;
A_aug=[A1,A2;eye(ell),zeros(ell)]; B_aug=[Bhat;zeros(ell,4)];
P_aug=[Phat,zeros(p,ell)];
model=struct('A1',A1,'A2',A2,'A',A_aug,'B_base',Bhat,'B',B_aug, ...
    'P_base',Phat,'P',P_aug,'R',Rhat,'Sigma_eps_base',Sigma_eps, ...
    'y_mean',y_mean,'u_mean',u_mean);

validation=evaluate_second_order(y_validation,u_validation,validation_segment_id, ...
    Eanchor,model);
collective=ones(4,1)/2;
collective_one_step=Eanchor'*Phat*Bhat*collective;
reach_horizon=18; d=10;
Ad=A_aug^d; Bd=zeros(size(B_aug)); Ai=eye(size(A_aug));
for k=1:d, Bd=Bd+Ai*B_aug; Ai=Ai*A_aug; end
collective_horizon=zeros(3,reach_horizon);
for j=1:reach_horizon
    collective_horizon(:,j)=Eanchor'*P_aug*(Ad^(j-1))*Bd*collective;
end
collective_horizon_gain=sqrt(sum(collective_horizon.^2,2));
all_input_horizon_gain=zeros(3,1);
for axis=1:3
    G=[];
    for j=1:reach_horizon
        Gj=zeros(1,reach_horizon*4);
        for i=1:j
            Gj(1,(i-1)*4+(1:4))=Eanchor(:,axis)'*P_aug*(Ad^(j-i))*Bd;
        end
        G=[G;Gj]; %#ok<AGROW>
    end
    all_input_horizon_gain(axis)=norm(G,'fro');
end

train_task=Eanchor'*y_train; task_bound=max(3.0,max(abs(train_task),[],'all'));
fitstats=anchor_stats; fitstats.order=2; fitstats.spectral_radius=max(abs(eig(A_aug)));
fitstats.collective_one_step=collective_one_step;
fitstats.collective_horizon_gain=collective_horizon_gain;
fitstats.all_input_horizon_gain=all_input_horizon_gain;
fitstats.transition_count_second_order=numel(valid);
results_dir=fullfile(here,'results'); if ~exist(results_dir,'dir'), mkdir(results_dir); end
fid=fopen(fullfile(results_dir,'copyBE_pelican_3d_z_aware_metrics.txt'),'w'); cleaner=onCleanup(@() fclose(fid));
fprintf(fid,'copyBE xyz-anchor second-order latent VARX\n');
fprintf(fid,'model_order 2\ncontrolled_directions standardized_x_y_z\n');
fprintf(fid,'dual_error %.12e\nanchor_error %.12e\n',fitstats.dual_error,fitstats.anchor_preservation_error);
fprintf(fid,'spectral_radius_augmented %.12f\n',fitstats.spectral_radius);
fprintf(fid,'task_RMSE %s\ntask_persistence_RMSE %s\n',mat2str(validation.task_rmse',10),mat2str(validation.task_persistence_rmse',10));
fprintf(fid,'collective_one_step_xyz %s\n',mat2str(collective_one_step',10));
fprintf(fid,'collective_horizon_gain_xyz %s\n',mat2str(collective_horizon_gain',10));
fprintf(fid,'all_input_horizon_gain_xyz %s\n',mat2str(all_input_horizon_gain',10));
clear cleaner;
save(fullfile(results_dir,'copyBE_pelican_3d_z_aware_data.mat'), ...
    'Eanchor','task_bound','fitstats','model','validation', ...
    'y_offset','y_scale','u_offset','u_scale','-v7.3');
fprintf('copyBE identification complete: RMSE=%s, persistence=%s, reach=%s\n', ...
    mat2str(validation.task_rmse',4),mat2str(validation.task_persistence_rmse',4), ...
    mat2str(all_input_horizon_gain',4));

function out=evaluate_second_order(y,u,run_id,E,model)
valid=find(run_id(1:end-2)==run_id(2:end-1) & run_id(2:end-1)==run_id(3:end))+1;
yc=y-model.y_mean; uc=u-model.u_mean; z=model.R'*yc;
zpred=model.A1*z(:,valid)+model.A2*z(:,valid-1)+model.B_base*uc(:,valid);
ypred=model.y_mean+model.P_base*zpred;
actual=E'*y(:,valid+1); predicted=E'*ypred; persistence=E'*y(:,valid);
out=struct('transition_count',numel(valid),'task_actual',actual,'task_predicted',predicted, ...
    'task_rmse',sqrt(mean((actual-predicted).^2,2)), ...
    'task_persistence_rmse',sqrt(mean((actual-persistence).^2,2)));
end
