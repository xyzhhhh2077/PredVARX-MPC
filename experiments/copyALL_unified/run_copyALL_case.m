function out=run_copyALL_case(method,base,D,cfg)
% Common evaluator. All methods share plant/data/SMPC; only identification/adaptation differs.
T=cfg.T_cl; x=zeros(cfg.n,1); y=zeros(cfg.p,T); u=zeros(cfg.m,T); ef=zeros(1,T); cc=nan(1,T); J=nan(1,T); fallback=0; pred_sq=[];
model=base; model.Sigma_obs=D.sigma_e_profile(1)^2*eye(cfg.p); zprev=[]; yonline=zeros(cfg.p,T); uonline=zeros(cfg.m,T); reid_count=0; bias=zeros(cfg.p,1);
win=max(cfg.noise_window,base.dinkla_window); eb=zeros(cfg.ell,win); ec=0; ob=zeros(cfg.p,cfg.noise_window); oc=0;
opt.N=cfg.N; opt.Q=cfg.Q; opt.Ru=cfg.Ru; opt.u_min=cfg.u_min; opt.u_max=cfg.u_max; opt.H=zeros(numel(cfg.tracked),cfg.p); opt.H(:,cfg.tracked)=eye(numel(cfg.tracked)); opt.h=cfg.y_max*ones(numel(cfg.tracked),1); opt.alpha_joint=cfg.alpha_joint;
for k=1:T
 vk=D.sigma_e_profile(k)*D.Vstd(:,k); yk=D.C*x+vk; y(:,k)=yk; z=model.R'*(yk-model.y_mean); IPR=eye(cfg.p)-model.P*model.R'; ores=IPR*(yk-model.y_mean); io=mod(k-1,cfg.noise_window)+1; ob(:,io)=ores; oc=min(oc+1,cfg.noise_window);
 if k>=2, er=z-model.A*zprev-model.B*(u(:,k-1)-model.u_mean); ie=mod(k-2,win)+1; eb(:,ie)=er; ec=min(ec+1,win); end
 cov_window=cfg.noise_window; if strcmp(base.covariance_mode,'dinkla'), cov_window=base.dinkla_window; end
 if ec>=5, ids=1:min(ec,cov_window); Qe=eb(:,ids)-mean(eb(:,ids),2); model.Sigma_eps=(Qe*Qe')/max(numel(ids)-1,1)+1e-8*eye(cfg.ell); model.Sigma_eps=(model.Sigma_eps+model.Sigma_eps')/2; end
 if oc>=5, Qo=ob(:,1:oc)-mean(ob(:,1:oc),2); so=norm(Qo,'fro')/sqrt(max((cfg.p-cfg.ell)*(oc-1),1)); model.Sigma_obs=max(so^2,1e-8)*eye(cfg.p); end
 rnow=D.Rf(:,k); if base.bias_correction, bias=(1-base.bias_gamma)*bias+base.bias_gamma*(yk-rnow); end
 rk=D.Rf(:,min(k+1,T)); if base.bias_correction, rk=rk-bias; end
 try, [~,~,U,info]=centered_smpc_step(yk,rk,model,opt); uk=U(1:cfg.m); ef(k)=info.exitflag; cc(k)=max(info.A_ch*U-info.b_ch); J(k)=info.cost; catch, fallback=fallback+1; uk=min(max(model.u_mean,cfg.u_min),cfg.u_max); ef(k)=-1; end
 if k<T, xnext=cfg.A*x+cfg.B*uk+D.sigma_w_profile(k)*D.Wstd(:,k); ynext=D.C*xnext+D.sigma_e_profile(k+1)*D.Vstd(:,k+1); pred=model.y_mean+model.P*(model.A*z+model.B*(uk-model.u_mean)); pred_sq(end+1)=mean((ynext-pred).^2); x=xnext; end %#ok<AGROW>
 u(:,k)=uk; yonline(:,k)=yk; uonline(:,k)=uk; zprev=z;
 if base.online_reidentify && mod(k,base.reidentify_period)==0
   switch base.reidentify_mode
     case 'sliding'
       first=max(1,k-base.reidentify_recent+1); yr=yonline(:,first:k); ur=uonline(:,first:k);
     case 'cumulative'
       yr=[D.y_off,yonline(:,1:k)]; ur=[D.u_off,uonline(:,1:k)];
     case 'offline_plus_recent'
       first=max(1,k-base.reidentify_recent+1); yr=[D.y_off,yonline(:,first:k)]; ur=[D.u_off,uonline(:,first:k)];
     otherwise
       error('Unknown reidentify mode: %s',base.reidentify_mode);
   end
   [Anew,Bnew,Pnew,Rnew,Snew,stnew]=control_aware_iterative_ivr_varx(yr,ur,cfg.ell,cfg.tracked);
   model.A=Anew; model.B=Bnew; model.P=Pnew; model.R=Rnew; model.Sigma_eps=Snew; model.y_mean=stnew.y_mean; model.u_mean=stnew.u_mean; model.Sigma_obs=max(trace(model.Sigma_obs)/cfg.p,1e-8)*eye(cfg.p);
   eb=zeros(cfg.ell,win); ec=0; ob=zeros(cfg.p,cfg.noise_window); oc=0; zprev=model.R'*(yk-model.y_mean); reid_count=reid_count+1;
 end
end
warm=cfg.warm_start:T; e=y(cfg.tracked,warm)-D.Rf(cfg.tracked,warm); E=zeros(cfg.p,numel(cfg.tracked)); E(cfg.tracked,:)=eye(numel(cfg.tracked));
out.method=method; out.MAE=mean(abs(e),2); out.RMSE=sqrt(mean(e.^2,2)); out.prediction_rmse=sqrt(mean(pred_sq)); out.reconstruction_residual=base.stats.reconstruction_residual; out.coverage_error=norm(base.P*base.R'*E-E,'fro'); out.dual_error=norm(base.R'*base.P-eye(cfg.ell),'fro'); out.avg_cost=mean(J(warm),'omitnan'); out.qp_success_rate=mean(ef>0); out.fallbacks=fallback; out.upper_violation_rate=sum(y(cfg.tracked,:)>cfg.y_max,2)/T; out.abs_violation_rate=sum(abs(y(cfg.tracked,:))>cfg.y_max,2)/T; out.u_rms=sqrt(mean(u.^2,2)); out.max_constraint=max(cc,[],'omitnan'); out.reidentify_count=reid_count; out.bias_final=bias(cfg.tracked); out.y=y; out.u=u;
end
