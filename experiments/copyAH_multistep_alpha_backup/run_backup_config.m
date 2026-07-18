function report = run_backup_config(cfg, plant, data)
% RUN_BACKUP_CONFIG Recovery modes:
%   none | soft_relaxed | original_alpha_multistep

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

seg=T_cl/5; Rf=zeros(p,T_cl);
rset=[0.6 -0.4; 1.2 0.3; -0.5 1.0; 0.9 -0.8; 0.2 0.5];
for s=1:5
    idx=(s-1)*seg+1:s*seg; Rf(1,idx)=rset(s,1); Rf(2,idx)=rset(s,2);
end
[swp,sep]=smooth_noise_profile(T_cl,noise_cycle,sw_min,sw_max,se_min,se_max,noise_phase_e);
rng(cfg.cl_seed,'twister');

x=zeros(n,1); y=zeros(p,T_cl); u=zeros(m,T_cl);
level=repmat({''},1,T_cl); stage=repmat({''},1,T_cl);
orig_cert_flag=false(1,T_cl);
nw=40; eb=zeros(ell,nw); ec=0; ob=zeros(p,nw); oc=0;
Ipr=eye(p)-Phat*Rhat'; zprev=[]; u_prev=[];
pf=0; n_primary=0; n_orig_rec=0; n_soft_rel=0; n_unc=0;
n_short=0; n_one=0; n_feas=0; n_redQ=0; n_prev=0; n_umean=0;
n_soft_also_orig=0;

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
        n_primary=n_primary+1; orig_cert_flag(k)=true;
    catch
        pf=pf+1;
        switch cfg.recovery_mode
            case 'soft_relaxed'
                fb=soft_recovery_smpc(yk,rk,model,opt);
                if fb.success
                    uk=fb.uk; level{k}='soft_relaxed'; stage{k}=fb.stage_class;
                    n_soft_rel=n_soft_rel+1;
                    % check if recovered uk constant-hold passes original alpha
                    cert=certify_input_original_alpha(uk,yk,model,opt);
                    orig_cert_flag(k)=cert.pass;
                    if cert.pass, n_soft_also_orig=n_soft_also_orig+1; end
                else
                    uk=min(max(model.u_mean,u_min),u_max);
                    level{k}='uncertified_fallback'; stage{k}='uncertified';
                    n_unc=n_unc+1; orig_cert_flag(k)=false;
                end
            case 'original_alpha_multistep'
                fb=multistep_original_alpha_recovery(yk,rk,model,opt,u_prev);
                uk=fb.uk; level{k}=fb.cc_cert_level; stage{k}=fb.stage_class;
                if fb.success
                    n_orig_rec=n_orig_rec+1; orig_cert_flag(k)=true;
                    switch fb.stage_class
                        case 'short_horizon_fullU', n_short=n_short+1;
                        case 'one_step_fullU', n_one=n_one+1;
                        case 'feasibility_fullU', n_feas=n_feas+1;
                        case 'reduced_Q_fullU', n_redQ=n_redQ+1;
                        case 'backup_u_prev', n_prev=n_prev+1;
                        case 'backup_u_mean', n_umean=n_umean+1;
                    end
                else
                    n_unc=n_unc+1; orig_cert_flag(k)=false;
                end
            otherwise
                uk=min(max(model.u_mean,u_min),u_max);
                level{k}='uncertified_fallback'; stage{k}='uncertified';
                n_unc=n_unc+1; orig_cert_flag(k)=false;
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
report.n_orig_rec=n_orig_rec; report.n_soft_rel=n_soft_rel; report.n_unc=n_unc;
report.orig_cert_rate=mean(orig_cert_flag);
report.uncert_rate=n_unc/T_cl;
report.soft_relaxed_rate=n_soft_rel/T_cl;
report.orig_rec_rate=n_orig_rec/T_cl;
report.orig_success_given_fail=n_orig_rec/max(pf,1);
report.soft_also_orig_rate=n_soft_also_orig/max(n_soft_rel,1);
report.n_short=n_short; report.n_one=n_one; report.n_feas=n_feas;
report.n_redQ=n_redQ; report.n_prev=n_prev; report.n_umean=n_umean;
report.left_err=stats.tracked_left_error;
end
