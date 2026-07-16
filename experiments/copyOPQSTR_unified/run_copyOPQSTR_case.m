function out=run_copyOPQSTR_case(method,base,D,cfg)
% Common closed-loop evaluator; only the method-specific identification policy differs.
T=cfg.T_cl; x=zeros(cfg.n,1); y=zeros(cfg.p,T); u=zeros(cfg.m,T); ef=zeros(1,T); cc=nan(1,T); J=nan(1,T); fallback=0; pred_sq=[];
model=base; model.Sigma_obs=D.sigma_e_profile(1)^2*eye(cfg.p); eb=zeros(cfg.ell,cfg.noise_window); ec=0; ob=zeros(cfg.p,cfg.noise_window); oc=0; zprev=[]; yhist=D.y_off; uhist=D.u_off; reid_count=0;
opt.N=cfg.N; opt.Q=cfg.Q; opt.Ru=cfg.Ru; opt.u_min=cfg.u_min; opt.u_max=cfg.u_max; opt.H=zeros(numel(cfg.tracked),cfg.p); opt.H(:,cfg.tracked)=eye(numel(cfg.tracked)); opt.h=cfg.y_max*ones(numel(cfg.tracked),1); opt.alpha_joint=cfg.alpha_joint;
for k=1:T
 vk=D.sigma_e_profile(k)*D.Vstd(:,k); yk=D.C*x+vk; y(:,k)=yk; z=model.R'*(yk-model.y_mean); IPR=eye(cfg.p)-model.P*model.R'; ores=IPR*(yk-model.y_mean); io=mod(k-1,cfg.noise_window)+1; ob(:,io)=ores; oc=min(oc+1,cfg.noise_window);
 if k>=2, er=z-model.A*zprev-model.B*(u(:,k-1)-model.u_mean); ie=mod(k-2,cfg.noise_window)+1; eb(:,ie)=er; ec=min(ec+1,cfg.noise_window); end
 if ec>=5, Qe=eb(:,1:ec)-mean(eb(:,1:ec),2); model.Sigma_eps=(Qe*Qe')/max(ec-1,1)+1e-8*eye(cfg.ell); model.Sigma_eps=(model.Sigma_eps+model.Sigma_eps')/2; end
 if oc>=5, Qo=ob(:,1:oc)-mean(ob(:,1:oc),2); so=norm(Qo,'fro')/sqrt(max((cfg.p-cfg.ell)*(oc-1),1)); model.Sigma_obs=max(so^2,1e-8)*eye(cfg.p); end
 rk=D.Rf(:,min(k+1,T));
 try, [~,yp,U,info]=centered_smpc_step(yk,rk,model,opt); uk=U(1:cfg.m); ef(k)=info.exitflag; cc(k)=max(info.A_ch*U-info.b_ch); J(k)=info.cost; catch, fallback=fallback+1; uk=min(max(model.u_mean,cfg.u_min),cfg.u_max); ef(k)=-1; yp=nan(cfg.p,1); end %#ok<NASGU>
 if k<T, xnext=cfg.A*x+cfg.B*uk+D.sigma_w_profile(k)*D.Wstd(:,k); ynext=D.C*xnext+D.sigma_e_profile(k+1)*D.Vstd(:,k+1); pred=model.y_mean+model.P*(model.A*z+model.B*(uk-model.u_mean)); pred_sq(end+1)=mean((ynext-pred).^2); x=xnext; end
 u(:,k)=uk; zprev=z; yhist=[yhist,yk]; uhist=[uhist,uk]; %#ok<AGROW>
 if base.online_reidentify && mod(k,base.reidentify_period)==0
   first=max(1,size(yhist,2)-cfg.reidentify_recent+1); yr=[D.y_off,yhist(:,first:end)]; ur=[D.u_off,uhist(:,first:end)];
   [Anew,Bnew,Pnew,Rnew,Snew,stnew]=control_aware_iterative_ivr_varx(yr,ur,cfg.ell,cfg.tracked);
   model.A=Anew; model.B=Bnew; model.P=Pnew; model.R=Rnew; model.Sigma_eps=Snew; model.y_mean=stnew.y_mean; model.u_mean=stnew.u_mean; model.Sigma_obs=max(trace(model.Sigma_obs)/cfg.p,1e-8)*eye(cfg.p); eb=zeros(cfg.ell,cfg.noise_window); ec=0; ob=zeros(cfg.p,cfg.noise_window); oc=0; zprev=model.R'*(yk-model.y_mean); reid_count=reid_count+1;
 end
end
warm=cfg.warm_start:T; e=y(cfg.tracked,warm)-D.Rf(cfg.tracked,warm); E=zeros(cfg.p,numel(cfg.tracked)); E(cfg.tracked,:)=eye(numel(cfg.tracked));
out.method=method; out.MAE=mean(abs(e),2); out.RMSE=sqrt(mean(e.^2,2)); out.prediction_rmse=sqrt(mean(pred_sq)); out.reconstruction_residual=base.stats.reconstruction_residual; out.coverage_error=norm(base.P*base.R'*E-E,'fro'); out.dual_error=norm(base.R'*base.P-eye(cfg.ell),'fro'); out.avg_cost=mean(J(warm),'omitnan'); out.qp_success_rate=mean(ef>0); out.fallbacks=fallback; out.upper_violation_rate=sum(y(cfg.tracked,:)>cfg.y_max,2)/T; out.abs_violation_rate=sum(abs(y(cfg.tracked,:))>cfg.y_max,2)/T; out.u_rms=sqrt(mean(u.^2,2)); out.max_constraint=max(cc,[],'omitnan'); out.reidentify_count=reid_count; out.y=y; out.u=u;
end
