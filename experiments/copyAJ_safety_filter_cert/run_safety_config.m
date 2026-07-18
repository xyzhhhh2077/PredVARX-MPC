function report = run_safety_config(cfg, plant, data)
% RUN_SAFETY_CONFIG Closed-loop recovery modes for copyAJ:
%   'none'                 - uncertified sat(u_mean) on primary fail
%   'soft_relaxed'         - old soft recovery (NOT original alpha)
%   'mean_safety_filter'   - L2' one-step deterministic mean hard-bound filter
%
% Certificate levels used (MUST NOT be named qp_original / original_alpha):
%   qp_primary           - primary N-step chance QP success
%   soft_relaxed         - soft recovery success
%   mean_safety_filter   - mean hard-bound safety filter success
%   uncertified_fallback - no certified recovery

root = fileparts(mfilename('fullpath'));
addpath(root); addpath(fullfile(root,'lib'));

n=plant.n; m=plant.m; p=plant.p; ell=plant.ell; tracked=plant.tracked;
A=plant.A; B=plant.B; C=plant.C; Sigma_n=plant.Sigma_n; L_n=plant.L_n;
se=plant.se; u_min=plant.u_min; u_max=plant.u_max;
noise_cycle=plant.noise_cycle;
sw_min=plant.sw_min; sw_max=plant.sw_max;
se_min=plant.se_min; se_max=plant.se_max; noise_phase_e=plant.noise_phase_e;

y_off=data.y_off; u_off=data.u_off;
T_cl=cfg.T_cl; N=cfg.N; y_max=cfg.y_max; alpha_joint=cfg.alpha_joint;

[Ahat,Bhat,Phat,Rhat,Sigma_eps,stats] = split_control_free_ivr_varx( ...
    y_off,u_off,ell,tracked,Sigma_n,'input_residualize',false);
[~,~,~,~,Szo,~] = cross_cov_diagnostics(y_off,Phat,Rhat);

model.A=Ahat; model.B=Bhat; model.P=Phat; model.R=Rhat;
model.y_mean=stats.y_mean; model.u_mean=stats.u_mean;
model.Sigma_eps=(Sigma_eps+Sigma_eps')/2;
model.Sigma_obs=Sigma_n; model.Sigma_zo=Szo;

opt.N=N; opt.Q=zeros(p); opt.Q(1,1)=80; opt.Q(2,2)=80; opt.Ru=0.25*eye(m);
opt.H=zeros(2,p); opt.H(1,1)=1; opt.H(2,2)=1; opt.h=y_max*[1;1];
opt.u_min=u_min; opt.u_max=u_max; opt.alpha_joint=alpha_joint;
opt.use_cross_cov=false; opt.use_terminal_cost=false;
opt.mean_bound_margin=0;

seg=T_cl/5; Rf=zeros(p,T_cl);
rset=[0.6 -0.4; 1.2 0.3; -0.5 1.0; 0.9 -0.8; 0.2 0.5];
for s=1:5
    idx=(s-1)*seg+1:s*seg; Rf(1,idx)=rset(s,1); Rf(2,idx)=rset(s,2);
end
[swp,sep]=smooth_noise_profile(T_cl,noise_cycle,sw_min,sw_max,se_min,se_max,noise_phase_e);
rng(cfg.cl_seed,'twister');

x=zeros(n,1); y=zeros(p,T_cl); u=zeros(m,T_cl);
level=repmat({''},1,T_cl); stage=repmat({''},1,T_cl);
nw=40; eb=zeros(ell,nw); ec=0; ob=zeros(p,nw); oc=0;
Ipr=eye(p)-Phat*Rhat'; zprev=[]; u_prev=[];
pf=0; n_primary=0; n_soft_rel=0; n_msf=0; n_unc=0;
n_msf_qp=0; n_msf_prev=0; n_msf_umean=0;
n_soft_also_mean=0;
n_msf_also_orig=0;  % diagnostic only: does mean filter uk pass original-alpha hold?

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
        model.Sigma_eps=(G/max(Nr-1,1))+1e-8*eye(ell);
        model.Sigma_eps=(model.Sigma_eps+model.Sigma_eps')/2;
    end
    if oc>=5
        model.Sigma_obs=build_sigma_obs_trial('declared_shape',ob(:,1:oc),Sigma_n,p,ell,1e-8);
    end
    rk=Rf(:,min(k+1,T_cl));
    try
        [~,~,U,info]=centered_smpc_step(yk,rk,model,opt);
        uk=U(1:m);
        mv=max(info.A_ch*U-info.b_ch);
        if ~(info.exitflag>0 && mv<=1e-7)
            error('primary row fail');
        end
        level{k}='qp_primary'; stage{k}='primary';
        n_primary=n_primary+1;
    catch
        pf=pf+1;
        switch cfg.recovery_mode
            case 'soft_relaxed'
                fb=soft_recovery_smpc(yk,rk,model,opt);
                if fb.success
                    uk=fb.uk; level{k}='soft_relaxed'; stage{k}=fb.stage_class;
                    n_soft_rel=n_soft_rel+1;
                    % diagnostic: does soft uk also pass mean hard bound?
                    ztmp=model.R'*(yk-model.y_mean);
                    G1=model.P*model.B;
                    mu0=model.y_mean+model.P*(model.A*ztmp)-G1*model.u_mean;
                    yhat=mu0+G1*uk;
                    H=opt.H; he=opt.h(:);
                    mv_mean=max([H*yhat-he; -H*yhat-he]);
                    if mv_mean<=1e-7
                        n_soft_also_mean=n_soft_also_mean+1;
                    end
                else
                    uk=min(max(model.u_mean,u_min),u_max);
                    level{k}='uncertified_fallback'; stage{k}='uncertified';
                    n_unc=n_unc+1;
                end
            case 'mean_safety_filter'
                fb=mean_safety_filter(yk,rk,model,opt,u_prev);
                uk=fb.uk;
                if fb.success
                    level{k}='mean_safety_filter'; stage{k}=fb.stage_class;
                    n_msf=n_msf+1;
                    switch fb.stage_class
                        case 'one_step_mean_qp', n_msf_qp=n_msf_qp+1;
                        case 'backup_u_prev_mean', n_msf_prev=n_msf_prev+1;
                        case 'backup_u_mean_mean', n_msf_umean=n_msf_umean+1;
                    end
                    % diagnostic only: original-alpha constant-hold (expect often fail)
                    try
                        cert=certify_input_original_alpha(uk,yk,model,opt);
                        if cert.pass, n_msf_also_orig=n_msf_also_orig+1; end
                    catch
                    end
                else
                    uk=min(max(model.u_mean,u_min),u_max);
                    level{k}='uncertified_fallback'; stage{k}='uncertified';
                    n_unc=n_unc+1;
                end
            otherwise % none
                uk=min(max(model.u_mean,u_min),u_max);
                level{k}='uncertified_fallback'; stage{k}='uncertified';
                n_unc=n_unc+1;
        end
    end
    u(:,k)=uk; u_prev=uk;
    x=A*x+B*uk+swp(k)*randn(n,1); zprev=zk;
end

warm=max(1,floor(0.2*T_cl)):T_cl;
err=y(tracked,warm)-Rf(tracked,warm);
report=struct();
report.name=cfg.name; report.seed=cfg.cl_seed; report.cfg=cfg;
report.MAE=mean(abs(err),2); report.RMSE=sqrt(mean(err.^2,2));
report.joint_viol=mean(any(y(tracked,:)>y_max,1));
report.joint_cover=1-report.joint_viol;
report.primary_fail=pf; report.n_primary=n_primary;
report.n_soft_rel=n_soft_rel; report.n_msf=n_msf; report.n_unc=n_unc;
report.primary_rate=n_primary/T_cl;
report.uncert_rate=n_unc/T_cl;
report.soft_relaxed_rate=n_soft_rel/T_cl;
report.msf_rate=n_msf/T_cl;
report.msf_success_given_fail=n_msf/max(pf,1);
report.soft_success_given_fail=n_soft_rel/max(pf,1);
report.soft_also_mean_rate=n_soft_also_mean/max(n_soft_rel,1);
report.msf_also_orig_rate=n_msf_also_orig/max(n_msf,1);
report.n_msf_qp=n_msf_qp; report.n_msf_prev=n_msf_prev; report.n_msf_umean=n_msf_umean;
report.n_soft_also_mean=n_soft_also_mean; report.n_msf_also_orig=n_msf_also_orig;
report.left_err=stats.tracked_left_error;
report.level=level; report.stage=stage;
report.cert_levels_used={'qp_primary','soft_relaxed','mean_safety_filter','uncertified_fallback'};
report.note=['L2_prime mean_safety_filter is one-step deterministic mean hard-bound; ' ...
    'NOT original-alpha chance recovery; NOT recursive feasibility.'];
end
