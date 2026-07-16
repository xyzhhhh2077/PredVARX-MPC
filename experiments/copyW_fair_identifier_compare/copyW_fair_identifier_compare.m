%% copyW_fair_identifier_compare -- copyV conditions, identifier-only ablation
clear; clc; close all; addpath(fileparts(mfilename('fullpath')));
rng(20260710,'twister');

n=6; m=3; p=30; ell=5; tracked=[1 2]; T_off=1500; T_cl=1200; N=18;
sw=0.045; se=0.055; noise_cycle=400;
sw_min=0.020; sw_max=0.090; se_min=0.025; se_max=0.100; noise_phase_e=pi/3;
u_min=-3; u_max=3; y_max=2.00; alpha_joint=0.10;
A=diag([0.94 0.88 0.78 0.64 0.50 0.35]);
A(1,2)=0.10; A(2,3)=-0.06; A(3,4)=0.05; A(4,5)=0.04;
B=[0.34 -0.10 0.05;0.12 0.28 -0.06;0.05 0.12 0.24;-0.05 0.06 0.18;0.02 -0.10 0.14;0.08 0.02 -0.08];
C=zeros(p,n); C(1,1)=1; C(1,3)=0.16; C(2,2)=1; C(2,4)=-0.12;
for i=3:p, C(i,:)=0.45*randn(1,n); C(i,:)=C(i,:)/max(norm(C(i,:)),1e-12); end

u_off=1.20*randn(m,T_off); x_off=zeros(n,T_off+1); y_off=zeros(p,T_off);
for k=1:T_off
    y_off(:,k)=C*x_off(:,k)+se*randn(p,1);
    x_off(:,k+1)=A*x_off(:,k)+B*u_off(:,k)+sw*randn(n,1);
end

phase=2*pi*(0:T_cl-1)/noise_cycle;
sigma_w_profile=sw_min+(sw_max-sw_min)*0.5*(1-cos(phase));
sigma_e_profile=se_min+(se_max-se_min)*0.5*(1-cos(phase+noise_phase_e));
Wstd=randn(n,T_cl); Vstd=randn(p,T_cl); % common random numbers for all cases
Rf=zeros(p,T_cl); levels=[0.25 1.50 0.65 1.85 0.45;0.35 1.25 1.75 0.80 1.55]; seg_len=floor(T_cl/5);
for s=1:5, ix=(s-1)*seg_len+1:min(s*seg_len,T_cl); Rf(1,ix)=levels(1,s); Rf(2,ix)=levels(2,s); end

names={'copyP_svd_ols','main_qr_ivr','copyO_oblique_ivr','copyV_control_aware_ivr'};
cases=cell(1,4);
[Ai,Bi,Pi,Ri,Si,sti]=control_ready_subspace_varx(y_off,u_off,ell); sti.tracked_projection_error=norm(Pi*Pi'*eye_cols(p,tracked)-eye_cols(p,tracked),'fro'); sti.reconstruction_residual=sti.subspace_residual;
cases{1}=pack_model(Ai,Bi,Pi,Ri,Si,sti,'copyP-style global SVD+OLS');
[Ai,Bi,Pi,~,Ri,~,~,Si]=predvarx_identify_main_qr(y_off,u_off,ell,1,1,A,B,C,n,m,p); sti=base_stats(y_off,u_off,Pi,Ri,tracked); cases{2}=pack_model(Ai,Bi,Pi,Ri,Si,sti,'main full QR-IVR');
[Ai,Bi,Pi,~,Ri,~,~,Si]=predvarx_identify_oblique(y_off,u_off,ell,1,1,A,B,C,n,m,p); sti=base_stats(y_off,u_off,Pi,Ri,tracked); cases{3}=pack_model(Ai,Bi,Pi,Ri,Si,sti,'copyO oblique dual-basis IVR');
[Ai,Bi,Pi,Ri,Si,sti]=control_aware_iterative_ivr_varx(y_off,u_off,ell,tracked); cases{4}=pack_model(Ai,Bi,Pi,Ri,Si,sti,'copyV tracked-axis iterative IVR');

opt.N=N; opt.Q=zeros(p); opt.Q(1,1)=80; opt.Q(2,2)=80; opt.Ru=0.18*eye(m); opt.u_min=u_min; opt.u_max=u_max; opt.H=zeros(numel(tracked),p); opt.H(:,tracked)=eye(numel(tracked)); opt.h=y_max*ones(numel(tracked),1); opt.alpha_joint=alpha_joint;
results=cell(1,4);
for ci=1:4
    results{ci}=run_case(cases{ci},A,B,C,Rf,Wstd,Vstd,sigma_w_profile,sigma_e_profile,opt,tracked,y_max,u_min,u_max);
    fprintf('%s: MAE=[%.4f %.4f] RMSE=[%.4f %.4f] predRMSE=%.4f coverage=%.3e QP=%.3f fallbacks=%d\n',names{ci},results{ci}.MAE,results{ci}.RMSE,results{ci}.prediction_rmse,cases{ci}.stats.tracked_projection_error,results{ci}.qp_success_rate,results{ci}.infeasible_count);
end

results_dir=fullfile(fileparts(mfilename('fullpath')),'results'); if ~exist(results_dir,'dir'), mkdir(results_dir); end
save(fullfile(results_dir,'copyW_fair_identifier_compare_data.mat'),'names','cases','results','A','B','C','Rf','u_off','x_off','y_off','Wstd','Vstd','sigma_w_profile','sigma_e_profile','opt','tracked','-v7.3');
fid=fopen(fullfile(results_dir,'copyW_fair_identifier_compare_metrics.csv'),'w');
fprintf(fid,'name,mae_y1,mae_y2,rmse_y1,rmse_y2,prediction_rmse,reconstruction_residual,tracked_projection_error,dual_error,avg_cost,qp_success,fallbacks,upper_viol_y1,upper_viol_y2,u_rms_1,u_rms_2,u_rms_3\n');
for ci=1:4
 r=results{ci}; st=cases{ci}.stats; fprintf(fid,'%s,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%d,%.12g,%.12g,%.12g,%.12g,%.12g\n',names{ci},r.MAE,r.RMSE,r.prediction_rmse,st.reconstruction_residual,st.tracked_projection_error,st.dual_error,r.avg_cost,r.qp_success_rate,r.infeasible_count,r.upper_violation_rate,r.u_rms); end
fclose(fid);

fig=figure('Position',[80 80 1800 1100],'Color','w'); t=1:T_cl; tl=tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
ax1=nexttile; plot(t,Rf(1,:),'k--','LineWidth',1.2); hold on; for ci=1:4, plot(t,results{ci}.y(1,:),'LineWidth',0.7); end; grid on; ylabel('y_1'); legend([{'r_1'},names],'Interpreter','none','Location','eastoutside');
ax2=nexttile; plot(t,Rf(2,:),'k--','LineWidth',1.2); hold on; for ci=1:4, plot(t,results{ci}.y(2,:),'LineWidth',0.7); end; grid on; ylabel('y_2'); legend([{'r_2'},names],'Interpreter','none','Location','eastoutside');
ax3=nexttile; vals=zeros(4,3); for ci=1:4, vals(ci,:)=[mean(results{ci}.MAE),results{ci}.prediction_rmse,results{ci}.avg_cost/1000]; end; bar(vals); grid on; set(gca,'XTickLabel',names,'XTickLabelRotation',15); legend({'mean MAE','one-step prediction RMSE','avg cost / 1000'}); title('Common-data identifier comparison');
linkaxes([ax1 ax2],'x'); title(tl,'copyW fair comparison: copyV plant, noise, reference, and SMPC'); print(fig,fullfile(results_dir,'copyW_fair_identifier_compare_fig'),'-dpng','-r160');

function E=eye_cols(p,tracked), E=zeros(p,numel(tracked)); E(tracked,:)=eye(numel(tracked)); end
function st=base_stats(y,u,P,R,tracked)
yc=y-mean(y,2); E=eye_cols(size(y,1),tracked); st.y_mean=mean(y,2); st.u_mean=mean(u,2); st.tracked_projection_error=norm(P*R'*E-E,'fro'); st.reconstruction_residual=norm(yc-P*(R'*yc),'fro')/norm(yc,'fro'); st.dual_error=norm(R'*P-eye(size(P,2)),'fro');
end
function c=pack_model(A,B,P,R,S,st,description)
if ~isfield(st,'dual_error'), st.dual_error=norm(R'*P-eye(size(P,2)),'fro'); end
c.A=A; c.B=B; c.P=P; c.R=R; c.Sigma_eps=(S+S')/2; c.y_mean=st.y_mean; c.u_mean=st.u_mean; c.stats=st; c.description=description;
end
function out=run_case(c,A,B,C,Rf,Wstd,Vstd,swp,sep,opt,tracked,ymax,umin,umax)
T=size(Rf,2); n=size(A,1); p=size(C,1); m=size(B,2); ell=size(c.A,1); x=zeros(n,1); y=zeros(p,T); yhat=zeros(p,T); u=zeros(m,T); ef=zeros(1,T); cc=nan(1,T); J=nan(1,T); fallback=0; pred_sq=[];
model=c; model.Sigma_obs=sep(1)^2*eye(p); noise_window=40; eb=zeros(ell,noise_window); ec=0; ob=zeros(p,noise_window); oc=0; IPR=eye(p)-c.P*c.R';
for k=1:T
 vk=sep(k)*Vstd(:,k); yk=C*x+vk; y(:,k)=yk; z=c.R'*(yk-c.y_mean); ores=IPR*(yk-c.y_mean); io=mod(k-1,noise_window)+1; ob(:,io)=ores; oc=min(oc+1,noise_window);
 if k>=2, er=z-model.A*zprev-model.B*(u(:,k-1)-model.u_mean); ie=mod(k-2,noise_window)+1; eb(:,ie)=er; ec=min(ec+1,noise_window); end
 if ec>=5, Q=eb(:,1:ec)-mean(eb(:,1:ec),2); model.Sigma_eps=(Q*Q')/max(ec-1,1)+1e-8*eye(ell); model.Sigma_eps=(model.Sigma_eps+model.Sigma_eps')/2; end
 if oc>=5, Q=ob(:,1:oc)-mean(ob(:,1:oc),2); so=norm(Q,'fro')/sqrt(max((p-ell)*(oc-1),1)); model.Sigma_obs=max(so^2,1e-8)*eye(p); end
 rk=Rf(:,min(k+1,T));
 try, [~,yp,U,info]=centered_smpc_step(yk,rk,model,opt); uk=U(1:m); yhat(:,k)=yp; ef(k)=info.exitflag; cc(k)=max(info.A_ch*U-info.b_ch); J(k)=info.cost; catch, fallback=fallback+1; uk=min(max(model.u_mean,umin),umax); ef(k)=-1; end
 if k<T, xnext=A*x+B*uk+swp(k)*Wstd(:,k); ynext=C*xnext+sep(k+1)*Vstd(:,k+1); pred_sq(end+1)=mean((ynext-(c.y_mean+c.P*(model.A*z+model.B*(uk-model.u_mean)))).^2); x=xnext; end
 u(:,k)=uk; zprev=z;
end
warm=151:T; e=y(tracked,warm)-Rf(tracked,warm); out.y=y; out.yhat=yhat; out.u=u; out.MAE=mean(abs(e),2); out.RMSE=sqrt(mean(e.^2,2)); out.prediction_rmse=sqrt(mean(pred_sq)); out.upper_violation_rate=sum(y(tracked,:)>ymax,2)/T; out.qp_success_rate=mean(ef>0); out.infeasible_count=fallback; out.avg_cost=mean(J(warm),'omitnan'); out.u_rms=sqrt(mean(u.^2,2)); out.max_constraint=max(cc,[],'omitnan');
end
