%% copyAU_soft_preference_output
% Learn two final output directions from old data using decaying preferences.
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
parts=strsplit(path,pathsep); path(strjoin(parts(~contains(lower(parts),'mosek')),pathsep));
addpath(fullfile(matlabroot,'toolbox','optim','optim'),'-begin');
here=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(here));
D=load(fullfile(repo,'experiments','copyAR_crte_paper_spectral_validation_unknown_noise', ...
    'results','copyAR_crte_paper_spectral_validation_unknown_noise_data.mat'));
y_off=D.y_off; u_off=D.u_off; A=D.A_true; B=D.B_true; C=D.C_true;
[p,T_off]=size(y_off); m=size(u_off,1); n=size(A,1); ell=D.ell; q=2;
T=D.T_cl; N=D.N; Ru=0.18*eye(m); limit=D.y_max;
weights=exp(-0.18*(0:p-1))'; preference_strength=0.78;
[Etask,dirstats]=learn_preferred_output_directions(y_off,u_off,q,struct( ...
    'weights',weights,'preference_strength',preference_strength,'reach_horizon',N,'Ru',Ru));
[Ahat,Bhat,Phat,Rhat,Sigma_eps,fitstats]=fit_anchored_varx(y_off,u_off,Etask,ell,struct('ridge',1e-8));
model=struct('A',Ahat,'B',Bhat,'P',Phat,'R',Rhat,'y_mean',fitstats.y_mean, ...
    'u_mean',fitstats.u_mean,'Sigma_eps',(Sigma_eps+Sigma_eps')/2);
O=(eye(p)-Phat*Rhat')*(y_off-fitstats.y_mean); O=O-mean(O,2);
model.Sigma_obs=O*O'/max(T_off-1,1)+1e-8*eye(p);
Q=80*(Etask*Etask'); opt=struct('N',N,'Q',Q,'Ru',Ru,'u_min',D.u_min, ...
    'u_max',D.u_max,'H',Etask','h',limit*ones(q,1),'alpha_joint',D.alpha_joint);
% Same five-segment scalar references, now attached to learned outputs s1,s2.
Sref=D.Rf(1:2,:);
rng(20260730,'twister');
[swp,sep]=smooth_noise_profile(T,D.noise_cycle,D.sw_min,D.sw_max,D.se_min,D.se_max,D.noise_phase_e);
Wnoise=randn(n,T).*swp; L=chol(D.Sigma_n_plant,'lower'); V0=randn(p,T); Vnoise=zeros(p,T);
for k=1:T, Vnoise(:,k)=(sep(k)/D.se)*L*V0(:,k); end
out=run_preference_branch(A,B,C,model,opt,Sref,Wnoise,Vnoise,Etask,limit);
assert(out.qp_success==1 && out.fallback==0,'Soft-preference 1200-step run failed.');
results_dir=fullfile(here,'results'); if ~exist(results_dir,'dir'),mkdir(results_dir);end
metrics_path=fullfile(results_dir,'copyAU_soft_preference_output_metrics.txt'); fid=fopen(metrics_path,'w');
fprintf(fid,'copyAU soft preference learned final outputs\nold_training_samples %d\nnew_training_samples 0\nclosed_loop_steps %d\n',T_off,T);
fprintf(fid,'weights %s\npreference_strength %.8f\n',mat2str(weights',6),preference_strength);
fprintf(fid,'direction_contribution %s\npreference_capture %.12f\n',mat2str(dirstats.contribution',8),dirstats.preference_capture);
fprintf(fid,'MAE %s\nRMSE %s\nBias %s\nQP_success %.12f\nfallback %d\nmax_qp %.12e\n', ...
    mat2str(out.MAE',10),mat2str(out.RMSE',10),mat2str(out.Bias',10),out.qp_success,out.fallback,out.max_qp); fclose(fid);
save(fullfile(results_dir,'copyAU_soft_preference_output_data.mat'),'Etask','weights','preference_strength', ...
    'dirstats','fitstats','out','model','opt','Wnoise','Vnoise','-v7.3');
fig=figure('Position',[60 40 2100 1600],'Color','w','Visible','off'); tl=tiledlayout(fig,5,1,'TileSpacing','compact'); tt=1:T;
ax=nexttile; plot(tt,Sref(1,:),'k--',tt,out.s(1,:),'b',tt,Sref(2,:),'k:',tt,out.s(2,:),'r'); yline(limit,'m--'); grid on; ylabel('learned outputs'); title(sprintf('Soft preference final outputs | MAE=%s',mat2str(out.MAE',3))); legend('r_1','s_1','r_2','s_2','limit','Location','eastoutside');
ax2=nexttile; plot(tt,out.estimated_sigma_eps,'b',tt,out.estimated_sigma_obs,'r',tt,out.true_sigma_w,'k--',tt,out.true_sigma_e,'Color',[.5 .2 .7]); grid on; ylabel('noise RMS'); title('Rolling estimates and realized noise'); legend('est innovation','est residual','realized w','realized v','Location','eastoutside');
ax3=nexttile; plot(tt,out.u'); grid on; ylabel('u'); title('Manipulated inputs');
ax4=nexttile; plot(tt,out.cost,'k'); grid on; ylabel('J'); title('Full MPC cost');
ax5=nexttile; plot(tt,out.maxcc,'k'); yline(0,'r--'); grid on; ylabel('max(AU-b)'); xlabel('time step'); title(sprintf('QP success=%.3f fallback=%d',out.qp_success,out.fallback));
linkaxes([ax ax2 ax3 ax4 ax5],'x'); title(tl,sprintf('copyAU soft preference | strength=%.2f | no new training set',preference_strength));
print(fig,fullfile(results_dir,'copyAU_soft_preference_output_fig.png'),'-dpng','-r160'); close(fig);
fprintf('copyAU complete: MAE=%s QP=%.3f fallback=%d\n',mat2str(out.MAE',4),out.qp_success,out.fallback);
