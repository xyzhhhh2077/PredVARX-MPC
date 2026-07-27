%% copyAS_learned_task_anchor_smpc
% Experimental extension: learn task-output anchor from a full-output task
% reference, then learn free CRTE-style directions and run the same SMPC.
% This is NOT the CRTE draft algorithm and has no stability/RF theorem.
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
parts=strsplit(path,pathsep); path(strjoin(parts(~contains(lower(parts),'mosek')),pathsep));
addpath(fullfile(matlabroot,'toolbox','optim','optim'),'-begin');
rng(20260728,'twister');

n=6; m=3; p=30; ell=5; q=2; T_off=1200; T_cl=700; N=18;
sw=0.045; se=0.055; u_min=-3; u_max=3; y_max=2.0;
A=diag([0.94 0.88 0.78 0.64 0.50 0.35]);
A(1,2)=0.10; A(2,3)=-0.06; A(3,4)=0.05; A(4,5)=0.04;
B=[0.34 -0.10 0.05;0.12 0.28 -0.06;0.05 0.12 0.24; ...
   -0.05 0.06 0.18;0.02 -0.10 0.14;0.08 0.02 -0.08];
C=randn(p,n); C(1,:)=[1 0 0.16 0 0 0]; C(2,:)=[0 1 0 -0.12 0 0];
for i=3:p, C(i,:)=C(i,:)/max(norm(C(i,:)),1e-12); end
sensor_rel=linspace(0.55,1.65,p)'; Dn=diag(se*sensor_rel);
Corr=eye(p); Corr(3,4)=0.30; Corr(4,3)=0.30; Corr(5,6)=-0.22; Corr(6,5)=-0.22;
Sigma_n_plant=Dn*Corr*Dn; Ln=chol(Sigma_n_plant,'lower');

u_off=1.2*randn(m,T_off); x_off=zeros(n,T_off+1); y_off=zeros(p,T_off);
for k=1:T_off
    y_off(:,k)=C*x_off(:,k)+Ln*randn(p,1);
    x_off(:,k+1)=A*x_off(:,k)+B*u_off(:,k)+sw*randn(n,1);
end

% Full-output task reference: deliberately mixed directions, no tracked indices.
t=linspace(0,10*pi,T_off);
Rtask_off=C*[0.9*sin(t);0.7*cos(0.7*t);0.5*sin(0.4*t+0.3); ...
    0.35*cos(0.25*t);zeros(2,T_off)];
[Ahat,Bhat,Phat,Rhat,Sigma_eps,stats]=learned_task_anchor_varx( ...
    y_off,u_off,ell,q,struct('task_reference',Rtask_off, ...
    'Ru',0.18*eye(m),'anchor_weights',[1 0.6 0.5 0.4], ...
    'mu_grid',[0 0.25 0.5 0.75 1],'reach_horizon',N));

model=struct('A',Ahat,'B',Bhat,'P',Phat,'R',Rhat, ...
    'y_mean',stats.y_mean,'u_mean',stats.u_mean, ...
    'Sigma_eps',(Sigma_eps+Sigma_eps')/2,'Sigma_obs',1e-6*eye(p));
opt.N=N; opt.Q=eye(p); opt.Ru=0.18*eye(m); opt.u_min=u_min; opt.u_max=u_max;
% Physical constraints remain on original outputs 1 and 2.
opt.H=zeros(2,p); opt.H(1,1)=1; opt.H(2,2)=1;
opt.h=y_max*ones(2,1); opt.alpha_joint=0.10;

Rf=zeros(p,T_cl); seg=floor(T_cl/5);
levels=[0.25 1.50 0.65 1.85 0.45;0.35 1.25 1.75 0.80 1.55];
for s=1:5
    ix=(s-1)*seg+1:min(s*seg,T_cl); Rf(1,ix)=levels(1,s); Rf(2,ix)=levels(2,s);
end
if seg*5<T_cl, Rf(1,seg*5+1:end)=levels(1,end); Rf(2,seg*5+1:end)=levels(2,end); end
opt.Q=zeros(p); opt.Q(1,1)=80; opt.Q(2,2)=80;

x=zeros(n,1); y=zeros(p,T_cl); u=zeros(m,T_cl); exitflag=zeros(1,T_cl);
maxcc=nan(1,T_cl); cost=nan(1,T_cl); fallback=0;
for k=1:T_cl
    yk=C*x+Ln*randn(p,1); y(:,k)=yk; rk=Rf(:,min(k+1,T_cl));
    try
        [~,~,U,info]=centered_smpc_step(yk,rk,model,opt);
        uk=U(1:m); exitflag(k)=info.exitflag; maxcc(k)=max(info.A_ch*U-info.b_ch); cost(k)=info.cost;
    catch ME
        fallback=fallback+1; exitflag(k)=-1; uk=zeros(m,1);
        if fallback<=3, fprintf('fallback k=%d: %s\n',k,ME.message); end
    end
    u(:,k)=uk; x=A*x+B*uk+sw*randn(n,1);
end
warm=101:T_cl; err=y(1:2,warm)-Rf(1:2,warm);
MAE=mean(abs(err),2); RMSE=sqrt(mean(err.^2,2)); Bias=mean(err,2);
qp_success=mean(exitflag>0); max_qp=max(maxcc,[],'omitnan'); cover=mean(exitflag(warm)>0);
origEc=zeros(p,2); origEc(1,1)=1; origEc(2,2)=1;
angles=acosd(min(max(svd(stats.E_task_anchor'*origEc),-1),1));

fprintf('\ncopyAS learned task anchor SMPC\n');
fprintf('anchor gap=%.3e selected mu=%.2f val=%.4f\n',stats.anchor_eigengap,stats.selected_mu,stats.selected_validation_nrmse);
fprintf('dual=%.3e preserve=%.3e angles-to-[e1,e2]=%s deg\n',stats.dual_error,stats.task_preservation_error,mat2str(angles',4));
fprintf('MAE=%s RMSE=%s Bias=%s\n',mat2str(MAE',4),mat2str(RMSE',4),mat2str(Bias',4));
fprintf('QP success=%.4f fallback=%d cover=%.4f maxQP=%.3e\n',qp_success,fallback,cover,max_qp);

results_dir=fullfile(fileparts(mfilename('fullpath')),'results'); if ~exist(results_dir,'dir'), mkdir(results_dir); end
metrics_path=fullfile(results_dir,'copyAS_learned_task_anchor_smpc_metrics.txt');
fid=fopen(metrics_path,'w');
fprintf(fid,'copyAS learned-task-anchor experimental extension\n');
fprintf(fid,'not_original_CRTE 1\nuses_true_Sigma_n 0\n');
fprintf(fid,'anchor_eigengap %.12g\nselected_mu %.12g\nvalidation_nrmse %.12g\n',stats.anchor_eigengap,stats.selected_mu,stats.selected_validation_nrmse);
fprintf(fid,'dual_error %.12g\ntask_preservation_error %.12g\n',stats.dual_error,stats.task_preservation_error);
fprintf(fid,'principal_angles_deg %s\n',mat2str(angles',8));
fprintf(fid,'MAE %s\nRMSE %s\nBias %s\n',mat2str(MAE',8),mat2str(RMSE',8),mat2str(Bias',8));
fprintf(fid,'qp_success_rate %.12g\nfallback_count %d\ncover %.12g\nmax_qp_constraint %.12g\n',qp_success,fallback,cover,max_qp);
fclose(fid);
save(fullfile(results_dir,'copyAS_learned_task_anchor_smpc_smoke.mat'), ...
    'Ahat','Bhat','Phat','Rhat','Sigma_eps','stats','MAE','RMSE','Bias', ...
    'qp_success','fallback','cover','max_qp','angles','-v7.3');
