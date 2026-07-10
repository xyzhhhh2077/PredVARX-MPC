%% copyQ_control_aware — reduced control-aware subspace SMPC
clear; clc; addpath(fileparts(mfilename('fullpath')));
rng(42,'twister');
n=6; m=3; p=15; ell=4; tracked=[1 2]; T_off=1000; T_cl=1200; N=15;
sw=0.05; se=0.10; u_min=-4; u_max=4; y_max=2.5;
A=diag([.85,.70,.55,.40,.30,.20]); A(1,2)=.10; A(3,4)=.01;
B=randn(n,m); C=randn(p,n); [C,~]=qr(C,0);

Rf=zeros(p,T_cl); rv=[0.4,1.8,0.3,1.8,0.6];
for s=1:5
    ix=(s-1)*240+1:min(s*240,T_cl); Rf(tracked,ix)=rv(s);
end
u_off=4*randn(m,T_off); x_off=zeros(n,T_off+1); y_off=zeros(p,T_off);
for k=1:T_off
    y_off(:,k)=C*x_off(:,k)+se*randn(p,1);
    x_off(:,k+1)=A*x_off(:,k)+B*u_off(:,k)+sw*randn(n,1);
end
[Ahat,Bhat,Phat,Rhat,Sigma_eps,pw]=control_aware_subspace_varx(y_off,u_off,ell,tracked);

model.A=Ahat; model.B=Bhat; model.P=Phat; model.R=Rhat;
model.y_mean=pw.y_mean; model.u_mean=pw.u_mean;
model.Sigma_eps=(Sigma_eps+Sigma_eps')/2; model.Sigma_obs=se^2*eye(p);
opt.N=N; opt.Q=zeros(p); opt.Q(tracked,tracked)=40*eye(numel(tracked)); opt.Ru=.20*eye(m);
opt.u_min=u_min; opt.u_max=u_max; opt.H=zeros(numel(tracked),p); opt.H(:,tracked)=eye(numel(tracked)); opt.h=y_max*ones(numel(tracked),1); opt.alpha_joint=.20;

x=zeros(n,1); y=zeros(p,T_cl); u=zeros(m,T_cl); yhat=zeros(p,T_cl); exitflag=zeros(1,T_cl); max_cc_violation=zeros(1,T_cl);
for k=1:T_cl
    yk=C*x+se*randn(p,1); y(:,k)=yk; rk=Rf(:,min(k+1,T_cl));
    [~,ypred,U,info]=centered_smpc_step(yk,rk,model,opt);
    u(:,k)=U(1:m); yhat(:,k)=ypred; exitflag(k)=info.exitflag; max_cc_violation(k)=max(info.A_ch*U-info.b_ch);
    x=A*x+B*u(:,k)+sw*randn(n,1);
end
warm=100:T_cl; err=y(tracked,warm)-Rf(tracked,warm); MAE=mean(abs(err),2); RMSE=sqrt(mean(err.^2,2));
upper_violation=sum(y(tracked,:)>y_max,2); abs_violation=sum(abs(y(tracked,:))>y_max,2);
fprintf('copyQ ell=%d: MAE=[%.3f %.3f], RMSE=[%.3f %.3f], upper=[%d %d], abs=[%d %d], maxQPviol=%.2e, projErr=%.2e\n',ell,MAE(1),MAE(2),RMSE(1),RMSE(2),upper_violation(1),upper_violation(2),abs_violation(1),abs_violation(2),max(max_cc_violation),pw.tracked_projection_error);

out.schema_version='copyQ_control_aware_v1'; out.A_true=A; out.B_true=B; out.C_true=C; out.tracked=tracked;
out.Ahat=Ahat; out.Bhat=Bhat; out.Phat=Phat; out.Rhat=Rhat; out.model=model; out.opt=opt; out.stats=pw;
out.x_off=x_off; out.y_off=y_off; out.u_off=u_off; out.y=y; out.u=u; out.yhat=yhat; out.Rf=Rf;
out.MAE=MAE; out.RMSE=RMSE; out.upper_violation=upper_violation; out.abs_violation=abs_violation; out.exitflag=exitflag; out.max_cc_violation=max_cc_violation;
save('copyQ_control_aware_data.mat','-struct','out','-v7.3');

figure('Position',[50 50 2000 1100],'Color','w'); t=1:T_cl;
subplot(4,1,1); plot(t,Rf(1,:),'k--',t,y(1,:),'b',t,Rf(2,:),'k:',t,y(2,:),'r'); yline(y_max,'m--'); grid on; title('copyQ reduced control-aware SMPC: outputs'); legend('r_1','y_1','r_2','y_2','y_{max}','Location','best');
subplot(4,1,2); plot(t,abs(y(1,:)-Rf(1,:)),'b',t,abs(y(2,:)-Rf(2,:)),'r'); grid on; title('absolute tracking error'); legend('|e_1|','|e_2|');
subplot(4,1,3); plot(t,u'); yline(u_min,'k--'); yline(u_max,'k--'); grid on; title('control input'); legend('u_1','u_2','u_3');
subplot(4,1,4); plot(t,max_cc_violation,'k'); yline(0,'r--'); grid on; title('max(A_{cc}U-b_{cc}): must be <= 0'); xlabel('time');
sgtitle(sprintf('copyQ ell=%d: MAE y_1=%.3f, y_2=%.3f; joint alpha=%.2f',ell,MAE(1),MAE(2),opt.alpha_joint));
print('copyQ_control_aware_fig','-dpng','-r150');
