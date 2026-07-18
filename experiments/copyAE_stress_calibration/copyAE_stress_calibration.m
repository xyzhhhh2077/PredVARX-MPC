%% copyAE_stress_calibration
% P0+ multi-seed soft-stage breakdown; P1+ tight-limit denom/obs coverage.
% Does not modify main/copyX/copyAA/copyAB/copyAC/copyAD.
clear; clc; close all;
root = fileparts(mfilename('fullpath'));
addpath(root); addpath(fullfile(root,'lib'));
rng(20260710,'twister');

%% Plant (same family)
plant=struct('n',6,'m',3,'p',30,'ell',5,'tracked',[1 2], ...
    'sw',0.045,'se',0.055,'noise_cycle',400, ...
    'sw_min',0.02,'sw_max',0.09,'se_min',0.025,'se_max',0.10,'noise_phase_e',pi/3, ...
    'u_min',-3,'u_max',3);
A=diag([0.94 0.88 0.78 0.64 0.50 0.35]);
A(1,2)=0.10; A(2,3)=-0.06; A(3,4)=0.05; A(4,5)=0.04;
B=[0.34 -0.10 0.05; 0.12 0.28 -0.06; 0.05 0.12 0.24; -0.05 0.06 0.18; 0.02 -0.10 0.14; 0.08 0.02 -0.08];
C=zeros(plant.p,plant.n); C(1,1)=1; C(1,3)=0.16; C(2,2)=1; C(2,4)=-0.12;
for i=3:plant.p
    C(i,:)=0.45*randn(1,plant.n); C(i,:)=C(i,:)/max(norm(C(i,:)),1e-12);
end
rel=linspace(0.55,1.65,plant.p)'; Corr=eye(plant.p);
Corr(3,4)=0.3; Corr(4,3)=0.3; Corr(5,6)=-0.22; Corr(6,5)=-0.22; Corr(8,9)=0.18; Corr(9,8)=0.18;
Sn=diag(plant.se*rel)*Corr*diag(plant.se*rel);
plant.A=A; plant.B=B; plant.C=C; plant.Sigma_n=Sn; plant.L_n=chol(Sn,'lower');

T_off=1500; u_off=1.2*randn(plant.m,T_off); x=zeros(plant.n,1); y_off=zeros(plant.p,T_off);
for k=1:T_off
    y_off(:,k)=C*x+plant.L_n*randn(plant.p,1);
    x=A*x+B*u_off(:,k)+plant.sw*randn(plant.n,1);
end
data.y_off=y_off; data.u_off=u_off;

base=struct('sigma_eps_mode','t2','use_cross_cov',false,'sigma_obs_mode','declared_shape', ...
    'enable_soft_recovery',true,'use_terminal_cost',false,'input_residualize',false, ...
    'y_max',0.55,'alpha_joint',0.10,'T_cl',600,'N',18,'cl_seed',1,'name','');

jobs={}; %#ok<*SAGROW>
seeds=[71001 71002 71003 71004 71005];

%% P0+: multi-seed soft OFF vs ON
for s=seeds
    c=base; c.name=sprintf('P0_no_soft_s%d',s); c.enable_soft_recovery=false; c.cl_seed=s; jobs{end+1}=c;
    c=base; c.name=sprintf('P0_soft_s%d',s); c.enable_soft_recovery=true; c.cl_seed=s; jobs{end+1}=c;
end

%% P1+: tight y_max coverage for eps modes (soft ON to keep runnable)
ym_list=[0.80 0.65 0.55];
for ym=ym_list
    for mode={'t2','ml','ols'}
        c=base; c.name=sprintf('P1_eps_%s_y%.2f',mode{1},ym);
        c.y_max=ym; c.sigma_eps_mode=mode{1}; c.enable_soft_recovery=true; c.cl_seed=82000+round(100*ym);
        jobs{end+1}=c;
    end
end

%% P1+: obs modes at y_max=0.65
for mode={'declared_shape','residual_support','additive'}
    c=base; c.name=sprintf('P1_obs_%s_y0.65',mode{1});
    c.y_max=0.65; c.sigma_obs_mode=mode{1}; c.enable_soft_recovery=true; c.cl_seed=83065;
    jobs{end+1}=c;
end

nJ=numel(jobs);
R=cell(nJ,1);
for i=1:nJ
    fprintf('[%d/%d] %s\n',i,nJ,jobs{i}.name);
    R{i}=run_seeded_config(jobs{i},plant,data);
    r=R{i};
    fprintf('  cover=%.4f MAE=%.3f/%.3f pf=%d soft=%d unc=%d stages risk/short/bound=%d/%d/%d soft_viol=%.4f\n', ...
        r.joint_cover,r.MAE(1),r.MAE(2),r.primary_fail,r.soft_ok,r.uncert, ...
        r.n_risk,r.n_short,r.n_bound,r.soft_step_viol);
end

resdir=fullfile(root,'results'); if ~exist(resdir,'dir'), mkdir(resdir); end
fid=fopen(fullfile(resdir,'copyAE_stress_calibration_metrics.csv'),'w');
fprintf(fid,['name,seed,y_max,eps,obs,soft_on,cover,viol,target,MAE1,MAE2,' ...
    'qp_rate,soft_rate,uncert_rate,primary_fail,soft_ok,uncert,' ...
    'soft_success_given_fail,n_risk,n_short,n_bound,soft_step_viol,uncert_step_viol\n']);
for i=1:nJ
    r=R{i}; c=r.cfg;
    fprintf(fid,['%s,%d,%.3f,%s,%s,%d,%.6f,%.6f,%.3f,%.6f,%.6f,' ...
        '%.6f,%.6f,%.6f,%d,%d,%d,%.6f,%d,%d,%d,%.6f,%.6f\n'], ...
        r.name,r.seed,c.y_max,c.sigma_eps_mode,c.sigma_obs_mode,c.enable_soft_recovery, ...
        r.joint_cover,r.joint_viol,r.target,r.MAE(1),r.MAE(2), ...
        r.qp_rate,r.soft_rate,r.uncert_rate,r.primary_fail,r.soft_ok,r.uncert, ...
        r.soft_success_given_fail,r.n_risk,r.n_short,r.n_bound,r.soft_step_viol,r.uncert_step_viol);
end
fclose(fid);

%% Aggregate P0
fmd=fopen(fullfile(resdir,'copyAE_stress_calibration_report.md'),'w');
fprintf(fmd,'# copyAE stress calibration (P0+ / P1+)\n\n');
fprintf(fmd,'Multi-seed soft-stage breakdown and tight-limit coverage probes. **Not theorem closures.**\n\n');
fprintf(fmd,'`T_cl=600`. Offline data `rng(20260710)`. Soft stages: risk_inflate / short_horizon / bound_only.\n\n');

fprintf(fmd,'## P0+ Soft multi-seed (y_max=0.55)\n\n');
fprintf(fmd,'| seed | soft | primary_fail | soft_ok | uncert | succ|fail | cover | soft_step_viol | risk | shortN | bound |\n');
fprintf(fmd,'|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
soft_succ=[]; no_unc=[];
for i=1:nJ
    if startsWith(R{i}.name,'P0_')
        r=R{i};
        fprintf(fmd,'| %d | %d | %d | %d | %d | %.3f | %.4f | %.4f | %d | %d | %d |\n', ...
            r.seed, r.cfg.enable_soft_recovery, r.primary_fail, r.soft_ok, r.uncert, ...
            r.soft_success_given_fail, r.joint_cover, r.soft_step_viol, r.n_risk, r.n_short, r.n_bound);
        if r.cfg.enable_soft_recovery
            soft_succ(end+1)=r.soft_success_given_fail; %#ok<AGROW>
            no_unc(end+1)=(r.uncert==0);
        end
    end
end
if ~isempty(soft_succ)
    fprintf(fmd,'\n**Aggregate soft-ON:** mean succ|fail = %.3f; fraction seeds with uncert=0: %.2f (%d/%d).\n', ...
        mean(soft_succ), mean(no_unc), sum(no_unc), numel(no_unc));
end
fprintf(fmd,'\nNon-claim: stage mix shows *how* soft recovers, not original-alpha validity.\n\n');

fprintf(fmd,'## P1+ Sigma_eps under tight y_max (soft ON)\n\n');
fprintf(fmd,'| config | y_max | cover | MAE1 | MAE2 | soft_rate | uncert_rate |\n|---|---:|---:|---:|---:|---:|---:|\n');
for i=1:nJ
    if startsWith(R{i}.name,'P1_eps_')
        r=R{i};
        fprintf(fmd,'| %s | %.2f | %.4f | %.4f | %.4f | %.3f | %.3f |\n', ...
            r.name,r.cfg.y_max,r.joint_cover,r.MAE(1),r.MAE(2),r.soft_rate,r.uncert_rate);
    end
end
fprintf(fmd,'\n');

fprintf(fmd,'## P1+ Sigma_obs at y_max=0.65\n\n');
fprintf(fmd,'| config | cover | MAE1 | MAE2 | soft_rate |\n|---|---:|---:|---:|---:|\n');
for i=1:nJ
    if startsWith(R{i}.name,'P1_obs_')
        r=R{i};
        fprintf(fmd,'| %s | %.4f | %.4f | %.4f | %.3f |\n', ...
            r.name,r.joint_cover,r.MAE(1),r.MAE(2),r.soft_rate);
    end
end
fprintf(fmd,'\n## Takeaways\n\n');
fprintf(fmd,'1. P0+: multi-seed confirms soft path reliability under stress when enabled.\n');
fprintf(fmd,'2. Stage histogram (risk/short/bound) documents recovery mix.\n');
fprintf(fmd,'3. P1+: if cover separates across denoms/modes at tight limits, use as empirical ranking only.\n');
fprintf(fmd,'4. Still open: original-alpha certificate; DOF theorem; Cov(o) identity; opinion 9 stability.\n');
fclose(fmd);

save(fullfile(resdir,'copyAE_stress_calibration_data.mat'),'R','plant','-v7.3');
fprintf('COPYAE_DONE\n');
