%% copyAT_learned_output_directions
% Same-data comparison of three output anchors:
% 1) fixed physical outputs e1,e2;
% 2) supervised learner that recovers the declared task-output span;
% 3) finite-horizon input-authority learner from the same old (u,y) data.
% No new training trajectories are generated. This is an experimental extension.
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
parts=strsplit(path,pathsep); path(strjoin(parts(~contains(lower(parts),'mosek')),pathsep));
addpath(fullfile(matlabroot,'toolbox','optim','optim'),'-begin');

here=fileparts(mfilename('fullpath'));
repo=fileparts(fileparts(here));
basefile=fullfile(repo,'experiments','copyAR_crte_paper_spectral_validation_unknown_noise', ...
    'results','copyAR_crte_paper_spectral_validation_unknown_noise_data.mat');
assert(exist(basefile,'file')==2,'Verified copyAR data artifact is required.');
D=load(basefile);

% Reuse the exact copyAR plant, old offline dataset, reference, and settings.
y_off=D.y_off; u_off=D.u_off; A=D.A_true; B=D.B_true; C=D.C_true;
[p,T_off]=size(y_off); m=size(u_off,1); n=size(A,1); ell=D.ell; q=2;
T_cl_full=D.T_cl; T_cl=min(T_cl_full,300); N=D.N; Rf=D.Rf(:,1:T_cl); tracked=D.tracked;
u_min=D.u_min; u_max=D.u_max; y_max=D.y_max; alpha_joint=D.alpha_joint;
Sigma_n_plant=D.Sigma_n_plant; L_n=chol(Sigma_n_plant,'lower');
Ru=0.18*eye(m); Q=zeros(p); Q(1,1)=80; Q(2,2)=80;
Ec=zeros(p,q); Ec(tracked,:)=eye(q);

% Learn both candidates from the same fixed old dataset.
[E_sup,dir_sup]=learn_output_directions(y_off,u_off,q, ...
    struct('mode','supervised','task_outputs',Ec,'reach_horizon',N,'Ru',Ru));
[E_auth,dir_auth]=learn_output_directions(y_off,u_off,q, ...
    struct('mode','authority','reach_horizon',N,'Ru',Ru));

anchors={Ec,E_sup,E_auth};
names={'fixed y1,y2','supervised output-span','input-authority'};
colors={[0.15 0.45 0.85],[0.10 0.65 0.35],[0.85 0.25 0.10]};
nb=numel(names); branches=cell(1,nb);

% Generate one common closed-loop disturbance realization, not training data.
rng(20260729,'twister');
[sigma_w_profile,sigma_e_profile]=smooth_noise_profile( ...
    T_cl,D.noise_cycle,D.sw_min,D.sw_max,D.se_min,D.se_max,D.noise_phase_e);
Wnoise=randn(n,T_cl).*sigma_w_profile;
Vstandard=randn(p,T_cl);
Vnoise=zeros(p,T_cl);
for k=1:T_cl
    Vnoise(:,k)=(sigma_e_profile(k)/D.se)*L_n*Vstandard(:,k);
end

for ib=1:nb
    if ib<=2
        model=D.model;
        fitstats=struct('y_mean',model.y_mean,'u_mean',model.u_mean, ...
            'dual_error',norm(model.R'*model.P-eye(ell),'fro'), ...
            'anchor_preservation_error',norm(model.P*model.R'*anchors{ib}-anchors{ib},'fro'), ...
            'spectral_radius',max(abs(eig(model.A))),'uses_new_training_data',false);
    else
        [Ahat,Bhat,Phat,Rhat,Sigma_eps,fitstats]=fit_anchored_varx( ...
            y_off,u_off,anchors{ib},ell,struct('ridge',1e-8));
        model=struct('A',Ahat,'B',Bhat,'P',Phat,'R',Rhat, ...
            'y_mean',fitstats.y_mean,'u_mean',fitstats.u_mean, ...
            'Sigma_eps',(Sigma_eps+Sigma_eps')/2);
        IminusPR=eye(p)-Phat*Rhat';
        O=IminusPR*(y_off-fitstats.y_mean); O=O-mean(O,2);
        Sobs=O*O'/max(T_off-1,1); Sobs=(Sobs+Sobs')/2;
        model.Sigma_obs=Sobs+1e-8*max(trace(Sobs)/p,1)*eye(p);
    end
    opt=struct('N',N,'Q',Q,'Ru',Ru,'u_min',u_min,'u_max',u_max, ...
        'H',Ec','h',y_max*ones(q,1),'alpha_joint',alpha_joint);
    branches{ib}=run_branch(A,B,C,model,opt,Rf,Wnoise,Vnoise,tracked,y_max);
    branches{ib}.name=names{ib}; branches{ib}.fitstats=fitstats;
    branches{ib}.anchor=anchors{ib};
    branches{ib}.angles_to_fixed=subspace_angles(anchors{ib},Ec);
end

% Direction-only diagnostics.
dir_fixed_authority=anchor_authority_fraction(Ec,dir_auth.W_y,dir_auth.C_y);
dir_sup_authority=anchor_authority_fraction(E_sup,dir_auth.W_y,dir_auth.C_y);
dir_auth_authority=dir_auth.authority_fraction;

fprintf('\ncopyAT same-old-data learned output directions\n');
fprintf('old training samples=%d; new training samples=0\n',T_off);
for ib=1:nb
    b=branches{ib};
    fprintf('%s: angles=%s MAE=%s RMSE=%s Bias=%s QP=%.4f fallback=%d maxQP=%.3e\n', ...
        names{ib},mat2str(b.angles_to_fixed',4),mat2str(b.MAE',4), ...
        mat2str(b.RMSE',4),mat2str(b.Bias',4),b.qp_success,b.fallback,b.max_qp);
end
fprintf('authority fractions fixed=%.4f supervised=%.4f authority=%.4f\n', ...
    dir_fixed_authority,dir_sup_authority,dir_auth_authority);

results_dir=fullfile(here,'results'); if ~exist(results_dir,'dir'), mkdir(results_dir); end
metrics_path=fullfile(results_dir,'copyAT_learned_output_directions_metrics.txt');
fid=fopen(metrics_path,'w');
fprintf(fid,'copyAT same-old-data learned output directions\n');
fprintf(fid,'not_original_CRTE 1\nuses_true_Sigma_n 0\n');
fprintf(fid,'old_training_samples %d\nnew_training_samples 0\n',T_off);
fprintf(fid,'common_closed_loop_disturbance 1\n');
fprintf(fid,'authority_fraction_fixed %.12f\n',dir_fixed_authority);
fprintf(fid,'authority_fraction_supervised %.12f\n',dir_sup_authority);
fprintf(fid,'authority_fraction_authority %.12f\n',dir_auth_authority);
for ib=1:nb
    b=branches{ib}; key={'fixed','supervised','authority'}; k=key{ib};
    fprintf(fid,'%s_angles_to_fixed_deg %s\n',k,mat2str(b.angles_to_fixed',8));
    fprintf(fid,'%s_MAE %s\n',k,mat2str(b.MAE',10));
    fprintf(fid,'%s_RMSE %s\n',k,mat2str(b.RMSE',10));
    fprintf(fid,'%s_Bias %s\n',k,mat2str(b.Bias',10));
    fprintf(fid,'%s_qp_success %.12f\n',k,b.qp_success);
    fprintf(fid,'%s_fallback %d\n',k,b.fallback);
    fprintf(fid,'%s_max_qp %.12e\n',k,b.max_qp);
    fprintf(fid,'%s_dual_error %.12e\n',k,b.fitstats.dual_error);
    fprintf(fid,'%s_anchor_preservation %.12e\n',k,b.fitstats.anchor_preservation_error);
end
fclose(fid);

save(fullfile(results_dir,'copyAT_learned_output_directions_data.mat'), ...
    'branches','dir_sup','dir_auth','Ec','E_sup','E_auth','names', ...
    'dir_fixed_authority','dir_sup_authority','dir_auth_authority','-v7.3');

fig=figure('Position',[50 50 2200 1500],'Color','w','Visible','off');
tlo=tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact'); t=1:T_cl;
for ib=1:3
    ax=nexttile(tlo,ib); b=branches{ib};
    plot(ax,t,Rf(1,:),'k--','LineWidth',1.0); hold(ax,'on');
    plot(ax,t,b.y(1,:),'Color',colors{ib},'LineWidth',0.75);
    plot(ax,t,Rf(2,:),'Color',[.25 .25 .25],'LineStyle',':','LineWidth',1.0);
    plot(ax,t,b.y(2,:),'Color',0.65*colors{ib},'LineWidth',0.75);
    yline(ax,y_max,'m--'); grid(ax,'on'); ylabel(ax,'output');
    title(ax,sprintf('%s | angle=%s deg | MAE=%s | QP=%.3f', ...
        names{ib},mat2str(b.angles_to_fixed',3),mat2str(b.MAE',3),b.qp_success));
    legend(ax,{'r_1','y_1','r_2','y_2','limit'},'Location','eastoutside');
end
ax=nexttile(tlo,4);
M=[branches{1}.MAE branches{2}.MAE branches{3}.MAE]';
bar(ax,M); grid(ax,'on'); ylabel(ax,'MAE'); xlabel(ax,'anchor method');
set(ax,'XTick',1:3,'XTickLabel',names); legend(ax,{'y_1','y_2'},'Location','eastoutside');
title(ax,sprintf('Same old training set; authority fractions=[%.3f %.3f %.3f]', ...
    dir_fixed_authority,dir_sup_authority,dir_auth_authority));
title(tlo,'copyAT: fixed versus learned output directions (no new training set)');
fig_path=fullfile(results_dir,'copyAT_learned_output_directions_fig.png');
print(fig,fig_path,'-dpng','-r160'); close(fig);
fprintf('metrics: %s\nfigure: %s\n',metrics_path,fig_path);

function out=run_branch(A,B,C,model,opt,Rf,Wnoise,Vnoise,tracked,y_max)
T=size(Rf,2); m=size(B,2); n=size(A,1); p=size(C,1); ell=size(model.A,1);
x=zeros(n,1); y=zeros(p,T); u=zeros(m,T); exitflag=zeros(1,T);
maxcc=nan(1,T); cost=nan(1,T); fallback=0; noise_window=40;
epsbuf=zeros(ell,noise_window); obsbuf=zeros(p,noise_window); ne=0; no=0;
IminusPR=eye(p)-model.P*model.R'; base_obs=max(mean(diag(model.Sigma_obs)),1e-8);
for k=1:T
    yk=C*x+Vnoise(:,k); y(:,k)=yk;
    zk=model.R'*(yk-model.y_mean); ores=IminusPR*(yk-model.y_mean);
    io=mod(k-1,noise_window)+1; obsbuf(:,io)=ores; no=min(no+1,noise_window);
    if k>=2
        er=zk-model.A*zprev-model.B*(u(:,k-1)-model.u_mean);
        ie=mod(k-2,noise_window)+1; epsbuf(:,ie)=er; ne=min(ne+1,noise_window);
    end
    if ne>=5
        E=epsbuf(:,1:ne); model.Sigma_eps=E*E'/ne+1e-8*eye(ell);
        model.Sigma_eps=(model.Sigma_eps+model.Sigma_eps')/2;
    end
    if no>=5
        O=obsbuf(:,1:no)-mean(obsbuf(:,1:no),2);
        so=norm(O,'fro')/sqrt(max((p-ell)*(no-1),1));
        model.Sigma_obs=max(so^2,1e-8)*base_obs*eye(p)+1e-8*eye(p);
    end
    rk=Rf(:,min(k+1,T));
    try
        [~,~,U,info]=centered_smpc_step(yk,rk,model,opt);
        uk=U(1:m); exitflag(k)=info.exitflag;
        maxcc(k)=max(info.A_ch*U-info.b_ch); cost(k)=info.cost;
    catch
        fallback=fallback+1; exitflag(k)=-1;
        uk=min(max(model.u_mean,opt.u_min),opt.u_max);
    end
    u(:,k)=uk; x=A*x+B*uk+Wnoise(:,k); zprev=zk;
end
warm=151:T; err=y(tracked,warm)-Rf(tracked,warm);
out.y=y; out.u=u; out.exitflag=exitflag; out.maxcc=maxcc; out.cost=cost;
out.MAE=mean(abs(err),2); out.RMSE=sqrt(mean(err.^2,2)); out.Bias=mean(err,2);
out.qp_success=mean(exitflag>0); out.fallback=fallback;
out.max_qp=max(maxcc,[],'omitnan');
out.upper_violation_rate=sum(y(tracked,:)>y_max,2)/T;
end

function a=subspace_angles(E,F)
[Qe,~]=qr(E,0); [Qf,~]=qr(F,0); s=svd(Qe'*Qf); s=min(max(s,-1),1); a=acosd(s);
end

function f=anchor_authority_fraction(E,Wy,Cy)
[Qe,~]=qr(E,0); En=Qe/real_sqrtm_local(Qe'*Cy*Qe);
num=sum(max(real(eig((En'*Wy*En+En'*Wy'*En)/2)),0));
[~,D]=eig((Wy+Wy')/2,Cy); den=sum(max(real(diag(D)),0));
f=num/max(den,1e-12);
end

function S=real_sqrtm_local(A)
A=(A+A')/2; [U,D]=eig(A); d=max(real(diag(D)),1e-12); S=U*diag(sqrt(d))*U';
end
