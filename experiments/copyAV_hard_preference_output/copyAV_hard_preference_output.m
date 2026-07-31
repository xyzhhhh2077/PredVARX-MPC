%% copyAV_hard_preference_output
% Lock the highest-weight original outputs as the two final outputs.
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
parts=strsplit(path,pathsep); path(strjoin(parts(~contains(lower(parts),'mosek')),pathsep));
addpath(fullfile(matlabroot,'toolbox','optim','optim'),'-begin');
here=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(here));
D=load(fullfile(repo,'experiments','copyAR_crte_paper_spectral_validation_unknown_noise', ...
    'results','copyAR_crte_paper_spectral_validation_unknown_noise_data.mat'));
y_off=D.y_off; u_off=D.u_off; A=D.A_true; B=D.B_true; C=D.C_true;
[p,T_off]=size(y_off); m=size(u_off,1); n=size(A,1); q=2; T=D.T_cl; N=D.N;
weights=exp(-0.18*(0:p-1))'; [Etask,prefstats]=select_hard_preference_outputs(weights,q);
% Hard output preference reuses the verified CRTE model because it locks y1,y2.
metric_mu=0.10;
[model,source_provenance]=validate_crte_source(D,metric_mu);
Ru=0.18*eye(m); Q=80*(Etask*Etask'); limit=D.y_max;
opt=struct('N',N,'Q',Q,'Ru',Ru,'u_min',D.u_min,'u_max',D.u_max, ...
    'H',Etask','h',limit*ones(q,1),'alpha_joint',D.alpha_joint);
Sref=D.Rf(prefstats.selected_indices,:);
rng(20260730,'twister');
[swp,sep]=smooth_noise_profile(T,D.noise_cycle,D.sw_min,D.sw_max,D.se_min,D.se_max,D.noise_phase_e);
Wnoise=randn(n,T).*swp; L=chol(D.Sigma_n_plant,'lower'); V0=randn(p,T); Vnoise=zeros(p,T);
for k=1:T, Vnoise(:,k)=(sep(k)/D.se)*L*V0(:,k); end
out=run_preference_branch(A,B,C,model,opt,Sref,Wnoise,Vnoise,Etask,limit);
assert(out.qp_success==1 && out.fallback==0,'Hard-preference 1200-step run failed.');
results_dir=fullfile(here,'results'); if ~exist(results_dir,'dir'),mkdir(results_dir);end
metrics_path=fullfile(results_dir,'copyAV_hard_preference_output_metrics.txt'); fid=fopen(metrics_path,'w');
fprintf(fid,'copyAV hard preference final outputs\nold_training_samples %d\nnew_training_samples 0\nclosed_loop_steps %d\n',T_off,T);
fprintf(fid,'weights %s\nselected_indices %s\nhard_locked 1\n',mat2str(weights',6),mat2str(prefstats.selected_indices));
fprintf(fid,'selected_mu %.12f\nntr_mode %s\nntr_epsilon %.12e\nntr_formula %s\n', ...
    source_provenance.selected_mu,source_provenance.ntr_mode, ...
    source_provenance.ntr_epsilon,source_provenance.ntr_formula);
fprintf(fid,'MAE %s\nRMSE %s\nBias %s\nQP_success %.12f\nfallback %d\nmax_qp %.12e\n', ...
    mat2str(out.MAE',10),mat2str(out.RMSE',10),mat2str(out.Bias',10),out.qp_success,out.fallback,out.max_qp); fclose(fid);
save(fullfile(results_dir,'copyAV_hard_preference_output_data.mat'),'Etask','weights','prefstats', ...
    'metric_mu','source_provenance', ...
    'out','model','opt','Wnoise','Vnoise','-v7.3');
fig=figure('Position',[60 40 2100 1600],'Color','w','Visible','off'); tl=tiledlayout(fig,5,1,'TileSpacing','compact'); tt=1:T;
ax=nexttile; plot(tt,Sref(1,:),'k--',tt,out.s(1,:),'b',tt,Sref(2,:),'k:',tt,out.s(2,:),'r'); yline(limit,'m--'); grid on; ylabel('locked outputs'); title(sprintf('Hard preference y_%d,y_%d | MAE=%s',prefstats.selected_indices(1),prefstats.selected_indices(2),mat2str(out.MAE',3))); legend('r_1','y_1','r_2','y_2','limit','Location','eastoutside');
ax2=nexttile; plot(tt,out.estimated_sigma_eps,'b',tt,out.estimated_sigma_obs,'r',tt,out.true_sigma_w,'k--',tt,out.true_sigma_e,'Color',[.5 .2 .7]); grid on; ylabel('noise RMS'); title('Rolling estimates and realized noise'); legend('est innovation','est residual','realized w','realized v','Location','eastoutside');
ax3=nexttile; plot(tt,out.u'); grid on; ylabel('u'); title('Manipulated inputs');
ax4=nexttile; plot(tt,out.cost,'k'); grid on; ylabel('J'); title('Full MPC cost');
ax5=nexttile; plot(tt,out.maxcc,'k'); yline(0,'r--'); grid on; ylabel('max(AU-b)'); xlabel('time step'); title(sprintf('QP success=%.3f fallback=%d',out.qp_success,out.fallback));
linkaxes([ax ax2 ax3 ax4 ax5],'x'); title(tl,sprintf( ...
    'copyAV hard preference | selected outputs=%s | CRTE mu=%.2f | corrected Ntr', ...
    mat2str(prefstats.selected_indices),metric_mu));
print(fig,fullfile(results_dir,'copyAV_hard_preference_output_fig.png'),'-dpng','-r160'); close(fig);
fprintf('copyAV complete: selected=%s MAE=%s QP=%.3f fallback=%d\n',mat2str(prefstats.selected_indices),mat2str(out.MAE',4),out.qp_success,out.fallback);
