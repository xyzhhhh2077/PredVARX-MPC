%% copyAZ_pybullet_identify
% Identify a copyAY-style soft-preference latent VARX on gym-pybullet-drones
% CF2X data (CtrlAviary, direct RPM). y = [pos(3); rpy(3); rpm(4)] mirrors the
% Pelican copyAY structure Pos+Euler+Motors; u = 4 commanded RPMs.
%
% This is phase-2 (closed-loop) validation: the simulator differs from the
% Pelican plant (CF2X 33 g vs Pelican 1.3 kg), so this validates the METHOD
% end-to-end on a dynamic simulator, not the Pelican plant itself.
%
% Output: results/copyAZ_pybullet_model.mat + metrics txt (one-step val R2).
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
here = fileparts(mfilename('fullpath'));
D = load(fullfile(here,'data','copyAZ_pybullet_dataset.mat'));
y = double(D.y); u = double(D.u); seg = double(D.segment_id(:)');

p = size(y,1); m = size(u,1); T = size(y,2);
q = 2; ELL = 5;   % latent dim like copyAY main (ell sweep optional)
train_mask = ismember(seg,1:4);   % 4 of 6 flights train
val_mask   = ismember(seg,5:6);   % 2 flights validation
fprintf('copyAZ: y %dx%d, u %dx%d, train %d samples, val %d samples\n', ...
    p,T,m,T,sum(train_mask),sum(val_mask));

% ---- scale with training statistics only ----
y_offset = mean(y(:,train_mask),2); y_scale = std(y(:,train_mask),0,2);
u_offset = mean(u(:,train_mask),2); u_scale = std(u(:,train_mask),0,2);
assert(all(y_scale>1e-6) && all(u_scale>1e-6),'zero training variance');
y_tr = (y(:,train_mask)-y_offset)./y_scale;
u_tr = (u(:,train_mask)-u_offset)./u_scale;

% ---- soft preference: position/attitude are the task, rpm internal ----
opt = struct('weights',[1 1 1 1 1 1 0.05 0.05 0.05 0.05], ...
    'preference_strength',0.7,'reach_horizon',18,'Ru',eye(m));
[E,learn_stats] = learn_segmented_preferred_output_directions( ...
    y_tr,u_tr,seg(train_mask),q,opt);
[Ahat,Bhat,Phat,Rhat,Sigma_eps,fit_stats] = fit_segmented_anchored_varx( ...
    y_tr,u_tr,seg(train_mask),E,ELL, ...
    struct('ridge',1e-8,'mu',0.10,'ntr_epsilon',10e-6));
model = struct('A',Ahat,'B',Bhat,'P',Phat,'R',Rhat, ...
    'Sigma_eps',Sigma_eps,'y_mean',fit_stats.y_mean,'u_mean',fit_stats.u_mean);
scales = struct('y_offset',y_offset,'y_scale',y_scale, ...
    'u_offset',u_offset,'u_scale',u_scale);
fprintf('fit: spectral radius %.6f, dual error %.2e\n', ...
    fit_stats.spectral_radius, fit_stats.dual_error);

% ---- one-step validation R2 (standardized, full output) ----
y_val = (y(:,val_mask)-y_offset)./y_scale;
u_val = (u(:,val_mask)-u_offset)./u_scale;
zval = model.R'*(y_val - model.y_mean);
yp1  = model.y_mean + model.P*(model.A*zval + model.B*u_val);
one_r2 = 1 - sum((y_val-yp1).^2,'all') / sum((y_val-model.y_mean).^2,'all');
fprintf('one-step validation R2 (ell=%d): %.4f\n',ELL,one_r2);

% ---- task-axis R2 (copyAY contract) ----
G = E'*E;
s_val = E'*y_val; sp1 = E'*yp1;
task_r2 = 1 - sum((s_val-sp1).^2,2)./sum((s_val-mean(s_val,2)).^2,2);
fprintf('task-axis R2: %s\n',mat2str(task_r2',4));

results_dir = fullfile(here,'results');
if ~exist(results_dir,'dir'), mkdir(results_dir); end
save(fullfile(results_dir,'copyAZ_pybullet_model.mat'), ...
    'model','scales','E','ELL','q','one_r2','task_r2', ...
    'fit_stats','learn_stats','-v7.3');

fid=fopen(fullfile(results_dir,'copyAZ_pybullet_metrics.txt'),'w');
c=onCleanup(@() fclose(fid));
fprintf(fid,'copyAZ pybullet CF2X identification (copyAY pipeline)\n');
fprintf(fid,'y [pos(3);rpy(3);rpm(4)] u [rpm_cmd(4)] ell %d\n',ELL);
fprintf(fid,'T %d train_flights 4 val_flights 2\n',T);
fprintf(fid,'one_step_validation_R2 %.6f\n',one_r2);
fprintf(fid,'task_axis_R2 %s\n',mat2str(task_r2',10));
fprintf(fid,'spectral_radius %.10f\ndual_error %.3e\n', ...
    fit_stats.spectral_radius, fit_stats.dual_error);
fprintf(fid,'note method-validation-only: CF2X dynamics != Pelican plant\n');
