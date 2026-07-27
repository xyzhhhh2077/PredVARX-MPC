function out=run_preference_branch(A,B,C,model,opt,Sref,Wnoise,Vnoise,Etask,limit)
%RUN_PREFERENCE_BRANCH Run one fixed-noise 1200-step learned-output SMPC.
T=size(Sref,2); m=size(B,2); n=size(A,1); p=size(C,1); ell=size(model.A,1);
x=zeros(n,1); y=zeros(p,T); s=zeros(size(Etask,2),T); u=zeros(m,T);
exitflag=zeros(1,T); maxcc=nan(1,T); cost=nan(1,T); fallback=0;
noise_window=40; epsbuf=zeros(ell,noise_window); obsbuf=zeros(p,noise_window);
ne=0; no=0; IminusPR=eye(p)-model.P*model.R';
base_obs=max(mean(diag(model.Sigma_obs)),1e-8);
est_eps=nan(1,T); est_obs=nan(1,T);
for k=1:T
    yk=C*x+Vnoise(:,k); y(:,k)=yk; s(:,k)=Etask'*yk;
    zk=model.R'*(yk-model.y_mean); ores=IminusPR*(yk-model.y_mean);
    io=mod(k-1,noise_window)+1; obsbuf(:,io)=ores; no=min(no+1,noise_window);
    if k>=2
        er=zk-model.A*zprev-model.B*(u(:,k-1)-model.u_mean);
        ie=mod(k-2,noise_window)+1; epsbuf(:,ie)=er; ne=min(ne+1,noise_window);
    end
    if ne>=5
        Er=epsbuf(:,1:ne)-mean(epsbuf(:,1:ne),2);
        model.Sigma_eps=Er*Er'/max(ne-1,1)+1e-8*eye(ell);
        model.Sigma_eps=(model.Sigma_eps+model.Sigma_eps')/2;
    end
    if no>=5
        Or=obsbuf(:,1:no)-mean(obsbuf(:,1:no),2);
        so=norm(Or,'fro')/sqrt(max((p-ell)*(no-1),1));
        model.Sigma_obs=max(so^2,1e-8)*base_obs*eye(p)+1e-8*eye(p);
    end
    est_eps(k)=sqrt(trace(model.Sigma_eps)/ell);
    est_obs(k)=sqrt(trace(model.Sigma_obs)/p);
    r_full=Etask*Sref(:,min(k+1,T));
    try
        [~,~,U,info]=centered_smpc_step(yk,r_full,model,opt);
        uk=U(1:m); exitflag(k)=info.exitflag;
        maxcc(k)=max(info.A_ch*U-info.b_ch); cost(k)=info.cost;
    catch
        fallback=fallback+1; exitflag(k)=-1;
        uk=min(max(model.u_mean,opt.u_min),opt.u_max);
    end
    u(:,k)=uk; x=A*x+B*uk+Wnoise(:,k); zprev=zk;
end
warm=151:T; err=s(:,warm)-Sref(:,warm);
out=struct('y',y,'s',s,'u',u,'Sref',Sref,'exitflag',exitflag, ...
    'maxcc',maxcc,'cost',cost,'estimated_sigma_eps',est_eps, ...
    'estimated_sigma_obs',est_obs,'true_sigma_w',vecnorm(Wnoise,2,1)/sqrt(n), ...
    'true_sigma_e',vecnorm(Vnoise,2,1)/sqrt(p));
out.MAE=mean(abs(err),2); out.RMSE=sqrt(mean(err.^2,2)); out.Bias=mean(err,2);
out.qp_success=mean(exitflag>0); out.fallback=fallback;
out.max_qp=max(maxcc,[],'omitnan'); out.violation_rate=mean(s>limit,2);
end
