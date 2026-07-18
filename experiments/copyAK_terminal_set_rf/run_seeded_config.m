function report = run_seeded_config(cfg, plant, data)
% RUN_SEEDED_CONFIG Closed-loop one config; optional terminal cost/set probes.
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

if ~isfield(cfg,'input_residualize'), cfg.input_residualize=false; end
if ~isfield(cfg,'use_cross_cov'), cfg.use_cross_cov=false; end
if ~isfield(cfg,'use_terminal_cost'), cfg.use_terminal_cost=false; end
if ~isfield(cfg,'use_terminal_set'), cfg.use_terminal_set=false; end
if ~isfield(cfg,'alpha_term'), cfg.alpha_term=NaN; end
if ~isfield(cfg,'terminal_box_scale'), cfg.terminal_box_scale=1.0; end
if ~isfield(cfg,'enable_soft_recovery'), cfg.enable_soft_recovery=true; end
if ~isfield(cfg,'sigma_eps_mode'), cfg.sigma_eps_mode='t2'; end
if ~isfield(cfg,'sigma_obs_mode'), cfg.sigma_obs_mode='declared_shape'; end

[Ahat,Bhat,Phat,Rhat,Sigma_eps,stats] = split_control_free_ivr_varx( ...
    y_off,u_off,ell,tracked,Sigma_n,'input_residualize',cfg.input_residualize);

switch lower(cfg.sigma_eps_mode)
    case 'ml', SigE = stats.Sigma_eps_ml;
    case 'ols', SigE = stats.Sigma_eps_ols;
    otherwise, SigE = Sigma_eps; cfg.sigma_eps_mode='t2';
end
[~,~,~,~,Szo,~] = cross_cov_diagnostics(y_off,Phat,Rhat);

model.A=Ahat; model.B=Bhat; model.P=Phat; model.R=Rhat;
model.y_mean=stats.y_mean; model.u_mean=stats.u_mean;
model.Sigma_eps=(SigE+SigE')/2;
model.Sigma_obs=Sigma_n;
model.Sigma_zo=Szo;

opt.N=N; opt.Q=zeros(p); opt.Q(1,1)=80; opt.Q(2,2)=80; opt.Ru=0.25*eye(m);
opt.H=zeros(2,p); opt.H(1,1)=1; opt.H(2,2)=1; opt.h=y_max*[1;1];
opt.u_min=u_min; opt.u_max=u_max; opt.alpha_joint=alpha_joint;
opt.use_cross_cov=cfg.use_cross_cov;
opt.use_terminal_cost=cfg.use_terminal_cost;
opt.use_terminal_set=cfg.use_terminal_set;
opt.alpha_term=cfg.alpha_term;
opt.terminal_box_scale=cfg.terminal_box_scale;

% Offline Lyapunov diagnostics (same Ahat as controller)
rho_A = max(abs(eig(Ahat)));
lyap_res = NaN; Vdec = NaN; Pterm_off = [];
try
    Qf = eye(ell);
    Pt = dlyap(Ahat', Qf); Pt=(Pt+Pt')/2;
    Pterm_off = Pt;
    lyap_res = norm(Ahat'*Pt*Ahat - Pt + Qf, 'fro');
    rng(11,'twister');
    ok=0; Nmc=200;
    for ii=1:Nmc
        zz=randn(ell,1); zz=zz/(norm(zz)+eps);
        if (Ahat*zz)'*Pt*(Ahat*zz) <= zz'*Pt*zz - 0.5*(zz'*Qf*zz) + 1e-9
            ok=ok+1;
        end
    end
    Vdec = ok/Nmc;
catch
end

% If terminal set requested without alpha, calibrate from free unit V quantile
if cfg.use_terminal_set && (~isfinite(cfg.alpha_term) || cfg.alpha_term<=0)
    if ~isempty(Pterm_off)
        cal = calibrate_alpha_term(Ahat, Pterm_off, 'Nmc', 600, 'seed', 19);
        % Scale by typical offline latent energy
        Zc = Rhat' * (y_off - stats.y_mean);
        vtyp = median(sum((Pterm_off*Zc).*Zc,1));  % median z'P z offline
        if ~isfinite(vtyp) || vtyp<=0
            vtyp = cal.alpha_recommend;
        end
        cfg.alpha_term = max(cal.alpha_recommend, 0.5*vtyp);
        opt.alpha_term = cfg.alpha_term;
    else
        cfg.use_terminal_set = false;
        opt.use_terminal_set = false;
    end
end

seg=T_cl/5; Rf=zeros(p,T_cl);
rset=[0.6 -0.4; 1.2 0.3; -0.5 1.0; 0.9 -0.8; 0.2 0.5];
for s=1:5
    idx=(s-1)*seg+1:s*seg; Rf(1,idx)=rset(s,1); Rf(2,idx)=rset(s,2);
end
[swp,sep]=smooth_noise_profile(T_cl,noise_cycle,sw_min,sw_max,se_min,se_max,noise_phase_e);
rng(cfg.cl_seed,'twister');

x=zeros(n,1); y=zeros(p,T_cl); u=zeros(m,T_cl);
cc=repmat({''},1,T_cl); stage=repmat({''},1,T_cl);
nw=40; eb=zeros(ell,nw); ec=0; ob=zeros(p,nw); oc=0;
Ipr=eye(p)-Phat*Rhat'; zprev=[];
pf=0; so=0; uc=0;
n_risk=0; n_short=0; n_bound=0;
Vterm_hist = nan(1,T_cl);
term_cost_flag = false; term_set_flag = false;

for k=1:T_cl
    vk=(sep(k)/se)*L_n*randn(p,1); yk=C*x+vk; y(:,k)=yk;
    zk=model.R'*(yk-model.y_mean);
    ores=Ipr*(yk-model.y_mean);
    io=mod(k-1,nw)+1; ob(:,io)=ores; oc=min(oc+1,nw);
    if k>=2
        er=zk-model.A*zprev-model.B*(u(:,k-1)-model.u_mean);
        ie=mod(k-2,nw)+1; eb(:,ie)=er; ec=min(ec+1,nw);
    end
    if ec>=5
        E=eb(:,1:ec)-mean(eb(:,1:ec),2); G=E*E'; Nr=ec;
        switch lower(cfg.sigma_eps_mode)
            case 'ml', d=max(Nr,1);
            case 'ols', d=max(Nr-(ell+m),1);
            otherwise, d=max(Nr-1,1);
        end
        model.Sigma_eps=(G/d)+1e-8*eye(ell);
        model.Sigma_eps=(model.Sigma_eps+model.Sigma_eps')/2;
    end
    if oc>=5
        model.Sigma_obs=build_sigma_obs_trial(cfg.sigma_obs_mode,ob(:,1:oc),Sigma_n,p,ell,1e-8);
    end
    rk=Rf(:,min(k+1,T_cl));
    try
        [~,~,U,info]=centered_smpc_step(yk,rk,model,opt);
        uk=U(1:m); cc{k}='qp'; stage{k}='primary';
        if isfield(info,'V_term'), Vterm_hist(k)=info.V_term; end
        if isfield(info,'terminal_cost_applied')
            term_cost_flag = term_cost_flag || info.terminal_cost_applied;
        end
        if isfield(info,'terminal_set_applied')
            term_set_flag = term_set_flag || info.terminal_set_applied;
        end
    catch
        pf=pf+1; ok=false;
        if cfg.enable_soft_recovery
            fb=soft_recovery_smpc(yk,rk,model,opt);
            if fb.success
                ok=true; so=so+1; uk=fb.uk; cc{k}='soft_recovery';
                if isfield(fb,'stage_class'), stage{k}=fb.stage_class; else, stage{k}=fb.mode; end
                switch stage{k}
                    case 'risk_inflate', n_risk=n_risk+1;
                    case 'short_horizon', n_short=n_short+1;
                    case 'bound_only', n_bound=n_bound+1;
                end
            end
        end
        if ~ok
            uc=uc+1; uk=min(max(model.u_mean,u_min),u_max);
            cc{k}='uncertified_fallback'; stage{k}='none';
        end
    end
    u(:,k)=uk; x=A*x+B*uk+swp(k)*randn(n,1); zprev=zk;
end

warm=max(1,floor(0.2*T_cl)):T_cl;
err=y(tracked,warm)-Rf(tracked,warm);
report=struct();
report.name=cfg.name; report.seed=cfg.cl_seed; report.cfg=cfg;
report.MAE=mean(abs(err),2); report.RMSE=sqrt(mean(err.^2,2));
report.joint_viol=mean(any(y(tracked,:)>y_max,1));
report.joint_cover=1-report.joint_viol;
report.target=1-alpha_joint;
report.qp_rate=mean(strcmp(cc,'qp'));
report.soft_rate=mean(strcmp(cc,'soft_recovery'));
report.uncert_rate=mean(strcmp(cc,'uncertified_fallback'));
report.qp_success=report.qp_rate;
report.primary_fail=pf; report.soft_ok=so; report.uncert=uc;
report.soft_success_given_fail=so/max(pf,1);
report.n_risk=n_risk; report.n_short=n_short; report.n_bound=n_bound;
sm=strcmp(cc,'soft_recovery');
if any(sm)
    report.soft_step_viol=mean(any(y(tracked,sm)>y_max,1));
else
    report.soft_step_viol=NaN;
end
um=strcmp(cc,'uncertified_fallback');
if any(um)
    report.uncert_step_viol=mean(any(y(tracked,um)>y_max,1));
else
    report.uncert_step_viol=NaN;
end
report.left_err=stats.tracked_left_error;
report.rho_A=rho_A;
report.lyap_residual_fro=lyap_res;
report.V_decrease_ok_ratio=Vdec;
report.alpha_term=cfg.alpha_term;
report.use_terminal_cost=cfg.use_terminal_cost;
report.use_terminal_set=cfg.use_terminal_set;
report.terminal_cost_applied=term_cost_flag;
report.terminal_set_applied=term_set_flag;
vv=Vterm_hist(isfinite(Vterm_hist));
if isempty(vv)
    report.V_term_median=NaN; report.V_term_p90=NaN; report.V_term_max=NaN;
else
    report.V_term_median=median(vv);
    report.V_term_p90=quantile(vv,0.90);
    report.V_term_max=max(vv);
end
report.Pterm_offline=Pterm_off;
end
