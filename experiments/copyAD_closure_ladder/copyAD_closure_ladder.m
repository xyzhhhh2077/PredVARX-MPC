%% copyAD_closure_ladder
% Sequential closure trials for remaining open opinions (P0->P1->P2).
% Does not modify main/, copyX, copyAA, copyAB, copyAC.
clear; clc; close all;
root_dir = fileparts(mfilename('fullpath'));
addpath(root_dir); addpath(fullfile(root_dir,'lib'));
rng(20260710,'twister');

%% Shared plant (same family as copyAB/AC)
plant = struct();
plant.n=6; plant.m=3; plant.p=30; plant.ell=5; plant.tracked=[1 2];
plant.sw=0.045; plant.se=0.055;
plant.noise_cycle=400;
plant.sw_min=0.020; plant.sw_max=0.090;
plant.se_min=0.025; plant.se_max=0.100; plant.noise_phase_e=pi/3;
plant.u_min=-3; plant.u_max=3;
A=diag([0.94 0.88 0.78 0.64 0.50 0.35]);
A(1,2)=0.10; A(2,3)=-0.06; A(3,4)=0.05; A(4,5)=0.04;
B=[0.34 -0.10 0.05; 0.12 0.28 -0.06; 0.05 0.12 0.24; -0.05 0.06 0.18; 0.02 -0.10 0.14; 0.08 0.02 -0.08];
C=zeros(plant.p,plant.n);
C(1,1)=1; C(1,3)=0.16; C(2,2)=1; C(2,4)=-0.12;
for i=3:plant.p
    C(i,:)=0.45*randn(1,plant.n);
    C(i,:)=C(i,:)/max(norm(C(i,:)),1e-12);
end
sensor_rel=linspace(0.55,1.65,plant.p)';
Corr=eye(plant.p);
Corr(3,4)=0.30; Corr(4,3)=0.30; Corr(5,6)=-0.22; Corr(6,5)=-0.22; Corr(8,9)=0.18; Corr(9,8)=0.18;
Dn=diag(plant.se*sensor_rel);
Sigma_n=Dn*Corr*Dn;
plant.A=A; plant.B=B; plant.C=C; plant.Sigma_n=Sigma_n; plant.L_n=chol(Sigma_n,'lower');

%% One offline dataset
T_off=1500;
u_off=1.20*randn(plant.m,T_off);
x=zeros(plant.n,1); y_off=zeros(plant.p,T_off);
for k=1:T_off
    y_off(:,k)=C*x+plant.L_n*randn(plant.p,1);
    x=A*x+B*u_off(:,k)+plant.sw*randn(plant.n,1);
end
data.y_off=y_off; data.u_off=u_off;

%% Config factory
base = struct('sigma_eps_mode','t2','use_cross_cov',false,'sigma_obs_mode','declared_shape', ...
    'enable_soft_recovery',false,'use_terminal_cost',false,'input_residualize',false, ...
    'y_max',2.00,'alpha_joint',0.10,'T_cl',800,'N',18,'cl_seed',20260718);

cfgs = {};

% ---- P0: soft recovery stress (tight y_max) ----
c=base; c.name='P0_stress_no_soft'; c.y_max=0.55; c.enable_soft_recovery=false; c.cl_seed=71001;
cfgs{end+1}=c; %#ok<*SAGROW>
c=base; c.name='P0_stress_with_soft'; c.y_max=0.55; c.enable_soft_recovery=true; c.cl_seed=71001;
cfgs{end+1}=c;
c=base; c.name='P0_stress_soft_plus_cross'; c.y_max=0.55; c.enable_soft_recovery=true; c.use_cross_cov=true; c.cl_seed=71001;
cfgs{end+1}=c;

% ---- P1: Sigma_eps denom coverage (nominal y_max) ----
for mode = {'t2','ml','ols'}
    c=base; c.name=['P1_eps_' mode{1}]; c.sigma_eps_mode=mode{1}; c.cl_seed=72001;
    cfgs{end+1}=c;
end

% ---- P1: Sigma_obs mode coverage ----
for mode = {'declared_shape','residual_support','additive'}
    c=base; c.name=['P1_obs_' mode{1}]; c.sigma_obs_mode=mode{1}; c.cl_seed=73001;
    cfgs{end+1}=c;
end

% ---- P2: cross on/off ----
c=base; c.name='P2_cross_off'; c.use_cross_cov=false; c.cl_seed=74001; cfgs{end+1}=c;
c=base; c.name='P2_cross_on'; c.use_cross_cov=true; c.cl_seed=74001; cfgs{end+1}=c;

% ---- P2: residualize on/off ----
c=base; c.name='P2_resid_off'; c.input_residualize=false; c.cl_seed=75001; cfgs{end+1}=c;
c=base; c.name='P2_resid_on'; c.input_residualize=true; c.cl_seed=75001; cfgs{end+1}=c;

% ---- Reference AB-like and AC-like short runs ----
c=base; c.name='REF_AB_like'; c.cl_seed=76001; cfgs{end+1}=c;
c=base; c.name='REF_AC_like'; c.sigma_eps_mode='ols'; c.use_cross_cov=true; ...
    c.sigma_obs_mode='additive'; c.enable_soft_recovery=true; c.use_terminal_cost=true; ...
    c.input_residualize=true; c.cl_seed=76001; cfgs{end+1}=c;

%% Run all
nC=numel(cfgs);
reports=cell(nC,1);
for i=1:nC
    fprintf('[%d/%d] %s ...\n', i, nC, cfgs{i}.name);
    reports{i}=run_one_config(cfgs{i}, plant, data);
    r=reports{i};
    fprintf('  MAE=[%.4f %.4f] joint_cover=%.4f (target~%.2f) qp=%.3f soft=%.3f unc=%.3f prim_fail=%d soft_ok=%d\n', ...
        r.MAE(1), r.MAE(2), r.joint_upper_cover_rate, r.nominal_cover_target, ...
        r.qp_rate, r.soft_rate, r.uncert_rate, r.primary_fail_count, r.soft_ok_count);
end

%% Write tables
resdir=fullfile(root_dir,'results');
if ~exist(resdir,'dir'), mkdir(resdir); end

fid=fopen(fullfile(resdir,'copyAD_closure_ladder_metrics.csv'),'w');
fprintf(fid,['name,y_max,eps,cross,obs,soft,terminal,resid,MAE1,MAE2,RMSE1,RMSE2,' ...
    'joint_cover,joint_viol,target_cover,qp_rate,soft_rate,uncert_rate,' ...
    'primary_fail,soft_ok,uncert,soft_success_given_fail,soft_step_viol,uncert_step_viol,active_rate,left_err\n']);
for i=1:nC
    r=reports{i}; c=r.cfg;
    fprintf(fid,['%s,%.4f,%s,%d,%s,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,' ...
        '%.6f,%.6f,%.4f,%.6f,%.6f,%.6f,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.3e\n'], ...
        r.name, c.y_max, c.sigma_eps_mode, c.use_cross_cov, c.sigma_obs_mode, ...
        c.enable_soft_recovery, c.use_terminal_cost, c.input_residualize, ...
        r.MAE(1), r.MAE(2), r.RMSE(1), r.RMSE(2), ...
        r.joint_upper_cover_rate, r.joint_upper_viol_rate, r.nominal_cover_target, ...
        r.qp_rate, r.soft_rate, r.uncert_rate, ...
        r.primary_fail_count, r.soft_ok_count, r.uncert_count, r.soft_success_given_fail, ...
        r.soft_step_joint_viol, r.uncert_step_joint_viol, r.active_rate, r.stats_geom_left);
end
fclose(fid);

% Markdown summary
fmd=fopen(fullfile(resdir,'copyAD_closure_ladder_report.md'),'w');
fprintf(fmd,'# copyAD closure ladder report\n\n');
fprintf(fmd,'Sequential empirical trials for remaining open opinions. **Not theorem closures.**\n\n');
fprintf(fmd,'Shared plant/seed family with copyAB/AC offline data `rng(20260710)`; each config has its own `cl_seed`.\n');
fprintf(fmd,'Horizon here `T_cl=800` (shorter ladder runs).\n\n');

fprintf(fmd,'## P0 Soft recovery under tight y_max=0.55\n\n');
fprintf(fmd,'| config | primary_fail | soft_ok | uncert | soft_success|fail | joint_cover | MAE |\n|---|---:|---:|---:|---:|---:|---:|\n');
for i=1:nC
    if startsWith(reports{i}.name,'P0_')
        r=reports{i};
        fprintf(fmd,'| %s | %d | %d | %d | %.3f | %.4f | %.4f/%.4f |\n', ...
            r.name, r.primary_fail_count, r.soft_ok_count, r.uncert_count, ...
            r.soft_success_given_fail, r.joint_upper_cover_rate, r.MAE(1), r.MAE(2));
    end
end
fprintf(fmd,'\nNon-claim: soft success is a **different** certificate (risk inflate / short N / bounds-only), not original alpha.\n\n');

fprintf(fmd,'## P1 Sigma_eps denominators (coverage proxy)\n\n');
fprintf(fmd,'| config | joint_cover | target | MAE1 | MAE2 | active |\n|---|---:|---:|---:|---:|---:|\n');
for i=1:nC
    if startsWith(reports{i}.name,'P1_eps_')
        r=reports{i};
        fprintf(fmd,'| %s | %.4f | %.2f | %.4f | %.4f | %.3f |\n', ...
            r.name, r.joint_upper_cover_rate, r.nominal_cover_target, r.MAE(1), r.MAE(2), r.active_rate);
    end
end
fprintf(fmd,'\nNon-claim: joint empirical cover is only a proxy; not a proof of correct residual DOF.\n\n');

fprintf(fmd,'## P1 Sigma_obs modes\n\n');
fprintf(fmd,'| config | joint_cover | MAE1 | MAE2 | active |\n|---|---:|---:|---:|---:|\n');
for i=1:nC
    if startsWith(reports{i}.name,'P1_obs_')
        r=reports{i};
        fprintf(fmd,'| %s | %.4f | %.4f | %.4f | %.3f |\n', ...
            r.name, r.joint_upper_cover_rate, r.MAE(1), r.MAE(2), r.active_rate);
    end
end
fprintf(fmd,'\nNon-claim: modes are engineering objects; none identified with Cov(o) theorem.\n\n');

fprintf(fmd,'## P2 Cross-cov and residualize\n\n');
fprintf(fmd,'| config | joint_cover | MAE1 | MAE2 | qp | active |\n|---|---:|---:|---:|---:|---:|\n');
for i=1:nC
    if startsWith(reports{i}.name,'P2_') || startsWith(reports{i}.name,'REF_')
        r=reports{i};
        fprintf(fmd,'| %s | %.4f | %.4f | %.4f | %.3f | %.3f |\n', ...
            r.name, r.joint_upper_cover_rate, r.MAE(1), r.MAE(2), r.qp_rate, r.active_rate);
    end
end
fprintf(fmd,'\n## Still open (theorem level)\n\n');
fprintf(fmd,'1. Recursive feasibility / stability (opinion 9) — not attempted beyond soft terminal cost elsewhere.\n');
fprintf(fmd,'2. Soft recovery original-risk certificate.\n');
fprintf(fmd,'3. Statistically optimal Sigma_eps denom; Sigma_obs = Cov(o); Boole with cross terms; input-conditional PredVARX optimum.\n');
fclose(fmd);

save(fullfile(resdir,'copyAD_closure_ladder_data.mat'),'reports','plant','-v7.3');
fprintf('\nWROTE results to %s\n', resdir);
fprintf('COPYAD_LADDER_DONE\n');
