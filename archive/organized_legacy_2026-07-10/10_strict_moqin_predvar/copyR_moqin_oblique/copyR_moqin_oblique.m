%% copyR_moqin_oblique -- strict Mo--Qin Algorithm-1 realization + centered SMPC
% This is a theory-faithful PredVAR baseline, not a control-constrained method.
% P,R,Pbar,Rbar follow direct normalized-space IVR de-normalization; main,
% copyO, copyP and copyQ are intentionally unchanged.
clear; clc; here=fileparts(mfilename('fullpath')); addpath(here); addpath(fullfile(fileparts(here),'copyP_centered_smpc'));
rng(42,'twister');
n=6; m=3; p=15; ell=4; T_off=1000; T_cl=1200; N=15; tracked=[1 2];
sw=0.05; se=0.10; u_min=-4; u_max=4; y_max=2.5;
A=diag([.85,.70,.55,.40,.30,.20]); A(1,2)=.10; A(3,4)=.01;
B=randn(n,m); C=randn(p,n); [C,~]=qr(C,0);
Rf=zeros(p,T_cl); rv=[0.4,1.8,0.3,1.8,0.6];
for s=1:5, ix=(s-1)*240+1:min(s*240,T_cl); Rf(tracked,ix)=rv(s); end

% Offline PE data, then Algorithm-1 normalized IVR and the explicit VARX extension.
u_off=4*randn(m,T_off); x_off=zeros(n,T_off+1); y_off=zeros(p,T_off);
for k=1:T_off
    y_off(:,k)=C*x_off(:,k)+se*randn(p,1);
    x_off(:,k+1)=A*x_off(:,k)+B*u_off(:,k)+sw*randn(n,1);
end
[Ahat,Bhat,Phat,Pbarhat,Rhat,Rbarhat,~,Sigma_eps,Sigma_ebar,~,~,~,info] = ...
    predvarx_identify_moqin(y_off,u_off,ell);
model.A=Ahat; model.B=Bhat; model.P=Phat; model.R=Rhat;
model.y_mean=info.y_mean; model.u_mean=info.u_mean;
model.Sigma_eps=(Sigma_eps+Sigma_eps')/2;
model.Sigma_obs=Pbarhat*Sigma_ebar*Pbarhat'; model.Sigma_obs=(model.Sigma_obs+model.Sigma_obs')/2;
opt.N=N; opt.Q=zeros(p); opt.Q(tracked,tracked)=40*eye(numel(tracked)); opt.Ru=.20*eye(m);
opt.u_min=u_min; opt.u_max=u_max; opt.H=zeros(numel(tracked),p); opt.H(:,tracked)=eye(numel(tracked));
opt.h=y_max*ones(numel(tracked),1); opt.alpha_joint=.20;

x=zeros(n,1); y=nan(p,T_cl); u=nan(m,T_cl); ypred=nan(p,T_cl);
exitflag=zeros(1,T_cl); max_cc_violation=nan(1,T_cl); projection_error=zeros(numel(tracked),1); failure_message=''; completed_steps=0;
for q=1:numel(tracked), ei=zeros(p,1); ei(tracked(q))=1; projection_error(q)=norm(Phat*(Rhat'*ei)-ei); end
for k=1:T_cl
    yk=C*x+se*randn(p,1); y(:,k)=yk; rk=Rf(:,min(k+1,T_cl));
    try
        [~,yp,U,step]=centered_smpc_step(yk,rk,model,opt);
    catch ME
        failure_message=ME.message; fprintf('copyR stops at k=%d: %s\n',k,ME.message); break;
    end
    u(:,k)=U(1:m); ypred(:,k)=yp; exitflag(k)=step.exitflag;
    max_cc_violation(k)=max(step.A_ch*U-step.b_ch);
    x=A*x+B*u(:,k)+sw*randn(n,1); completed_steps=k;
end
warm=100:completed_steps;
if isempty(warm)
    MAE=nan(numel(tracked),1); RMSE=nan(numel(tracked),1); upper_violation=nan(numel(tracked),1); abs_violation=nan(numel(tracked),1);
else
    err=y(tracked,warm)-Rf(tracked,warm); MAE=mean(abs(err),2); RMSE=sqrt(mean(err.^2,2));
    upper_violation=sum(y(tracked,1:completed_steps)>y_max,2); abs_violation=sum(abs(y(tracked,1:completed_steps))>y_max,2);
end
fprintf('copyR strict Mo-Qin ell=%d: completed=%d/%d, MAE=[%.4f %.4f], QP=%d/%d, cover=[%.2e %.2e], dual=[%.2e %.2e %.2e %.2e]\n', ...
 ell,completed_steps,T_cl,MAE(1),MAE(2),sum(exitflag>0),T_cl,projection_error(1),projection_error(2),info.dual_errors);

out.schema_version='copyR_moqin_oblique_v1'; out.A_true=A; out.B_true=B; out.C_true=C; out.tracked=tracked;
out.Ahat=Ahat; out.Bhat=Bhat; out.Phat=Phat; out.Pbarhat=Pbarhat; out.Rhat=Rhat; out.Rbarhat=Rbarhat; out.info=info; out.model=model; out.opt=opt;
out.x_off=x_off; out.y_off=y_off; out.u_off=u_off; out.y=y; out.u=u; out.ypred=ypred; out.Rf=Rf;
out.MAE=MAE; out.RMSE=RMSE; out.upper_violation=upper_violation; out.abs_violation=abs_violation;
out.exitflag=exitflag; out.max_cc_violation=max_cc_violation; out.projection_error=projection_error; out.completed_steps=completed_steps; out.failure_message=failure_message;
results_dir=fullfile(fileparts(mfilename('fullpath')),'results'); if ~exist(results_dir,'dir'), mkdir(results_dir); end
save(fullfile(results_dir,'copyR_moqin_oblique_data.mat'),'-struct','out','-v7.3');

figure('Position',[50 50 2200 1200],'Color','w'); t=1:T_cl;
subplot(4,1,1); plot(t,Rf(1,:),'k--',t,y(1,:),'b',t,Rf(2,:),'k:',t,y(2,:),'r'); yline(y_max,'m--'); grid on; title('copyR strict Mo--Qin PredVARX-SMPC: tracked outputs'); legend('r_1','y_1','r_2','y_2','y_{max}','Location','best');
subplot(4,1,2); plot(t,abs(y(1,:)-Rf(1,:)),'b',t,abs(y(2,:)-Rf(2,:)),'r'); grid on; title('absolute tracking error'); legend('|e_1|','|e_2|');
subplot(4,1,3); plot(t,u'); yline(u_min,'k--'); yline(u_max,'k--'); grid on; title('control input'); legend('u_1','u_2','u_3');
subplot(4,1,4); plot(t,max_cc_violation,'k'); yline(0,'r--'); grid on; title('max(A_{cc}U-b_{cc}): must be <= 0'); xlabel('time');
sgtitle(sprintf('copyR strict Mo--Qin: ell=%d; MAE=[%.3f %.3f]; coverage=[%.2e %.2e]',ell,MAE(1),MAE(2),projection_error(1),projection_error(2)));
print(fullfile(results_dir,'copyR_moqin_oblique_fig'),'-dpng','-r150');
