function report = run_one_config(cfg, plant, data)
% RUN_ONE_CONFIG Single closed-loop config for copyAD ladder studies.
% cfg fields:
%   name, sigma_eps_mode, use_cross_cov, sigma_obs_mode,
%   enable_soft_recovery, use_terminal_cost, input_residualize,
%   y_max, alpha_joint, T_cl, N, stress_label

root = fileparts(mfilename('fullpath'));
addpath(root); addpath(fullfile(root,'lib'));

n=plant.n; m=plant.m; p=plant.p; ell=plant.ell; tracked=plant.tracked;
A=plant.A; B=plant.B; C=plant.C; Sigma_n=plant.Sigma_n; L_n=plant.L_n;
se=plant.se; sw=plant.sw;
u_min=plant.u_min; u_max=plant.u_max;
noise_cycle=plant.noise_cycle;
sw_min=plant.sw_min; sw_max=plant.sw_max;
se_min=plant.se_min; se_max=plant.se_max; noise_phase_e=plant.noise_phase_e;

y_off=data.y_off; u_off=data.u_off;
T_cl=cfg.T_cl; N=cfg.N; y_max=cfg.y_max; alpha_joint=cfg.alpha_joint;

[Ahat,Bhat,Phat,Rhat,Sigma_eps,stats] = split_control_free_ivr_varx( ...
    y_off,u_off,ell,tracked,Sigma_n, ...
    'input_residualize', cfg.input_residualize);

switch lower(cfg.sigma_eps_mode)
    case 'ml', Sigma_eps_use = stats.Sigma_eps_ml;
    case 'ols', Sigma_eps_use = stats.Sigma_eps_ols;
    otherwise, Sigma_eps_use = Sigma_eps; cfg.sigma_eps_mode='t2';
end

[~,~,~,~,Szo,~] = cross_cov_diagnostics(y_off, Phat, Rhat);

model.A=Ahat; model.B=Bhat; model.P=Phat; model.R=Rhat;
model.y_mean=stats.y_mean; model.u_mean=stats.u_mean;
model.Sigma_eps=(Sigma_eps_use+Sigma_eps_use')/2;
model.Sigma_obs=Sigma_n;
model.Sigma_zo=Szo;

opt.N=N;
opt.Q=zeros(p); opt.Q(1,1)=80; opt.Q(2,2)=80;
opt.Ru=0.25*eye(m);
opt.H=zeros(2,p); opt.H(1,1)=1; opt.H(2,2)=1;
opt.h=y_max*[1;1];
opt.u_min=u_min; opt.u_max=u_max;
opt.alpha_joint=alpha_joint;
opt.use_cross_cov=cfg.use_cross_cov;
opt.use_terminal_cost=cfg.use_terminal_cost;

seg_len=T_cl/5;
Rf=zeros(p,T_cl);
rset=[0.6 -0.4; 1.2 0.3; -0.5 1.0; 0.9 -0.8; 0.2 0.5];
for s=1:5
    idx=(s-1)*seg_len+1:s*seg_len;
    Rf(1,idx)=rset(s,1); Rf(2,idx)=rset(s,2);
end
[sigma_w_profile, sigma_e_profile] = smooth_noise_profile( ...
    T_cl, noise_cycle, sw_min, sw_max, se_min, se_max, noise_phase_e);

% independent stream for CL but fixed seed offset per config name
rng(cfg.cl_seed,'twister');

x=zeros(n,1);
y=zeros(p,T_cl); u=zeros(m,T_cl);
exitflag=nan(1,T_cl); max_cc=nan(1,T_cl);
cc_level=repmat({''},1,T_cl);
soft_mode=repmat({''},1,T_cl);
noise_window=40;
eps_buffer=zeros(ell,noise_window); eps_count=0;
obs_buffer=zeros(p,noise_window); obs_count=0;
IminusPR=eye(p)-Phat*Rhat';
z_prev=[];
primary_fail=0; soft_ok=0; uncert=0;

for k=1:T_cl
    vk=(sigma_e_profile(k)/se)*L_n*randn(p,1);
    yk=C*x+vk;
    y(:,k)=yk;
    zk=model.R'*(yk-model.y_mean);
    obs_res=IminusPR*(yk-model.y_mean);
    io=mod(k-1,noise_window)+1;
    obs_buffer(:,io)=obs_res;
    obs_count=min(obs_count+1,noise_window);
    if k>=2
        eps_res=zk-model.A*z_prev-model.B*(u(:,k-1)-model.u_mean);
        ie=mod(k-2,noise_window)+1;
        eps_buffer(:,ie)=eps_res;
        eps_count=min(eps_count+1,noise_window);
    end
    if eps_count>=5
        Ewin=eps_buffer(:,1:eps_count)-mean(eps_buffer(:,1:eps_count),2);
        Geps=Ewin*Ewin'; Nres=eps_count;
        switch lower(cfg.sigma_eps_mode)
            case 'ml', den=max(Nres,1);
            case 'ols', den=max(Nres-(ell+m),1);
            otherwise, den=max(Nres-1,1);
        end
        model.Sigma_eps=(Geps/den)+1e-8*eye(ell);
        model.Sigma_eps=(model.Sigma_eps+model.Sigma_eps')/2;
    end
    if obs_count>=5
        Owin=obs_buffer(:,1:obs_count);
        model.Sigma_obs=build_sigma_obs_trial(cfg.sigma_obs_mode,Owin,Sigma_n,p,ell,1e-8);
    end
    rk=Rf(:,min(k+1,T_cl));
    try
        [~,~,U,info]=centered_smpc_step(yk,rk,model,opt);
        uk=U(1:m);
        exitflag(k)=info.exitflag;
        max_cc(k)=max(info.A_ch*U-info.b_ch);
        cc_level{k}='qp';
        soft_mode{k}='primary';
    catch
        primary_fail=primary_fail+1;
        recovered=false;
        if cfg.enable_soft_recovery
            fb2=soft_recovery_smpc(yk,rk,model,opt);
            if fb2.success
                recovered=true; soft_ok=soft_ok+1;
                uk=fb2.uk; exitflag(k)=fb2.exitflag;
                max_cc(k)=fb2.max_cc_violation;
                cc_level{k}='soft_recovery'; soft_mode{k}=fb2.mode;
            end
        end
        if ~recovered
            uncert=uncert+1;
            uk=min(max(model.u_mean,u_min),u_max);
            exitflag(k)=-1; max_cc(k)=NaN;
            cc_level{k}='uncertified_fallback'; soft_mode{k}='none';
        end
    end
    u(:,k)=uk;
    wk=sigma_w_profile(k)*randn(n,1);
    x=A*x+B*uk+wk;
    z_prev=zk;
end

warm=151:T_cl;
err=y(tracked,warm)-Rf(tracked,warm);
report=struct();
report.name=cfg.name;
report.cfg=cfg;
report.stats_geom_left=stats.tracked_left_error;
report.stats_geom_dual=stats.dual_error;
report.ivr_input_conditional=stats.ivr_input_conditional;
report.MAE=mean(abs(err),2);
report.RMSE=sqrt(mean(err.^2,2));
report.Bias=mean(err,2);
% empirical one-sided upper violation on tracked outputs over full horizon
report.upper_viol_rate=sum(y(tracked,:)>y_max,2)/T_cl;
report.abs_viol_rate=sum(abs(y(tracked,:))>y_max,2)/T_cl;
% joint any-tracked upper violation rate (coverage proxy)
report.joint_upper_viol_rate=mean(any(y(tracked,:)>y_max,1));
report.joint_upper_cover_rate=1-report.joint_upper_viol_rate;
report.nominal_cover_target=1-alpha_joint; % rough Boole joint target
report.qp_rate=mean(strcmp(cc_level,'qp'));
report.soft_rate=mean(strcmp(cc_level,'soft_recovery'));
report.uncert_rate=mean(strcmp(cc_level,'uncertified_fallback'));
report.primary_fail_count=primary_fail;
report.soft_ok_count=soft_ok;
report.uncert_count=uncert;
report.soft_success_given_fail=soft_ok/max(primary_fail,1);
qp_mask=strcmp(cc_level,'qp');
if any(qp_mask)
    report.max_qp_resid=max(max_cc(qp_mask & ~isnan(max_cc)));
    warm_qp = false(1,T_cl); warm_qp(warm)=true;
    msk = qp_mask & warm_qp & ~isnan(max_cc);
    if any(msk)
        report.active_rate=mean(max_cc(msk)>-1e-3);
    else
        report.active_rate=NaN;
    end
else
    report.max_qp_resid=NaN; report.active_rate=NaN;
end
% soft-step empirical safety: among soft steps, fraction still under y_max next? use same-step y
soft_mask=strcmp(cc_level,'soft_recovery');
if any(soft_mask)
    report.soft_step_joint_viol=mean(any(y(tracked,soft_mask)>y_max,1));
else
    report.soft_step_joint_viol=NaN;
end
unc_mask=strcmp(cc_level,'uncertified_fallback');
if any(unc_mask)
    report.uncert_step_joint_viol=mean(any(y(tracked,unc_mask)>y_max,1));
else
    report.uncert_step_joint_viol=NaN;
end
end
