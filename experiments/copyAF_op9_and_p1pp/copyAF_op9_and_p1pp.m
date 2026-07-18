%% copyAF_op9_and_p1pp
% 1) Opinion 9 Lyapunov/Schur numerical probes (NOT RF/stability proofs)
% 2) P1++ multi-seed tight y_max coverage for eps/obs modes
% Does not modify main/copyX/copyAA/copyAB/copyAC/copyAD/copyAE.
clear; clc; close all;
root = fileparts(mfilename('fullpath'));
addpath(root); addpath(fullfile(root,'lib'));
rng(20260710,'twister');

%% Plant
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

resdir=fullfile(root,'results'); if ~exist(resdir,'dir'), mkdir(resdir); end

%% ========== Opinion 9 probe ==========
fprintf('=== Opinion 9 Lyapunov/Schur probe ===\n');
op9 = opinion9_lyapunov_probe(y_off, u_off, plant.ell, plant.tracked, Sn);
fprintf('rho(A)=%.6f Schur=%d dlyap_ok=%d lyap_res=%.3e Vdec_ok=%.3f ctrb_rank=%d left=%.2e\n', ...
    op9.rho_A, op9.is_schur, op9.dlyap_ok, op9.lyap_residual_fro, ...
    op9.V_decrease_ok_ratio, op9.ctrb_rank, op9.left_error);

% Closed-loop terminal ON vs OFF under moderate stress (not RF proof)
base=struct('sigma_eps_mode','t2','use_cross_cov',false,'sigma_obs_mode','declared_shape', ...
    'enable_soft_recovery',true,'use_terminal_cost',false,'input_residualize',false, ...
    'y_max',0.70,'alpha_joint',0.10,'T_cl',500,'N',18,'cl_seed',91001,'name','');
c0=base; c0.name='OP9_term_OFF'; c0.use_terminal_cost=false;
c1=base; c1.name='OP9_term_ON'; c1.use_terminal_cost=true;
fprintf('OP9 CL term OFF...\n'); r0=run_seeded_config(c0,plant,data);
fprintf('OP9 CL term ON...\n'); r1=run_seeded_config(c1,plant,data);
fprintf('term OFF: cover=%.4f MAE=%.4f/%.4f soft=%.3f\n', r0.joint_cover,r0.MAE(1),r0.MAE(2),r0.soft_rate);
fprintf('term ON : cover=%.4f MAE=%.4f/%.4f soft=%.3f\n', r1.joint_cover,r1.MAE(1),r1.MAE(2),r1.soft_rate);

%% ========== P1++ multi-seed tight coverage ==========
fprintf('\n=== P1++ multi-seed tight coverage ===\n');
jobs={}; %#ok<*SAGROW>
seeds=[92001 92002 92003 92004 92005];
ym_list=[0.50 0.45];
for ym=ym_list
    for mode={'t2','ml','ols'}
        for s=seeds
            c=base; c.name=sprintf('P1pp_eps_%s_y%.2f_s%d',mode{1},ym,s);
            c.y_max=ym; c.sigma_eps_mode=mode{1}; c.cl_seed=s; c.T_cl=500;
            c.enable_soft_recovery=true; jobs{end+1}=c;
        end
    end
end
% obs modes multi-seed at y=0.50
for mode={'declared_shape','residual_support','additive'}
    for s=seeds
        c=base; c.name=sprintf('P1pp_obs_%s_y0.50_s%d',mode{1},s);
        c.y_max=0.50; c.sigma_obs_mode=mode{1}; c.cl_seed=s; c.T_cl=500;
        c.enable_soft_recovery=true; jobs{end+1}=c;
    end
end

nJ=numel(jobs); RR=cell(nJ,1);
for i=1:nJ
    fprintf('[%d/%d] %s\n',i,nJ,jobs{i}.name);
    RR{i}=run_seeded_config(jobs{i},plant,data);
    r=RR{i};
    fprintf('  cover=%.4f MAE=%.3f/%.3f soft=%.2f unc=%.2f\n', ...
        r.joint_cover,r.MAE(1),r.MAE(2),r.soft_rate,r.uncert_rate);
end

%% CSV
fid=fopen(fullfile(resdir,'copyAF_metrics.csv'),'w');
fprintf(fid,['name,seed,y_max,eps,obs,term,cover,viol,MAE1,MAE2,soft_rate,uncert_rate,' ...
    'primary_fail,soft_ok,n_risk,n_short,n_bound,soft_step_viol\n']);
% op9 cl rows
for r=[r0 r1]
    c=r.cfg;
    fprintf(fid,'%s,%d,%.3f,%s,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,%.6f\n', ...
        r.name,r.seed,c.y_max,c.sigma_eps_mode,c.sigma_obs_mode,c.use_terminal_cost, ...
        r.joint_cover,r.joint_viol,r.MAE(1),r.MAE(2),r.soft_rate,r.uncert_rate, ...
        r.primary_fail,r.soft_ok,r.n_risk,r.n_short,r.n_bound,r.soft_step_viol);
end
for i=1:nJ
    r=RR{i}; c=r.cfg;
    fprintf(fid,'%s,%d,%.3f,%s,%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,%.6f\n', ...
        r.name,r.seed,c.y_max,c.sigma_eps_mode,c.sigma_obs_mode,c.use_terminal_cost, ...
        r.joint_cover,r.joint_viol,r.MAE(1),r.MAE(2),r.soft_rate,r.uncert_rate, ...
        r.primary_fail,r.soft_ok,r.n_risk,r.n_short,r.n_bound,r.soft_step_viol);
end
fclose(fid);

%% Markdown report with aggregates
fmd=fopen(fullfile(resdir,'copyAF_op9_and_p1pp_report.md'),'w');
fprintf(fmd,'# copyAF: Opinion 9 probes + P1++ coverage\n\n');
fprintf(fmd,'**Not theorem closures.**\n\n');

fprintf(fmd,'## Opinion 9 numerical probe\n\n');
fprintf(fmd,'| quantity | value |\n|---|---:|\n');
fprintf(fmd,'| rho(Ahat) | %.6f |\n', op9.rho_A);
fprintf(fmd,'| is_schur | %d |\n', op9.is_schur);
fprintf(fmd,'| dlyap_ok | %d |\n', op9.dlyap_ok);
fprintf(fmd,'| lyap residual F | %.3e |\n', op9.lyap_residual_fro);
fprintf(fmd,'| V-decrease OK ratio | %.3f |\n', op9.V_decrease_ok_ratio);
fprintf(fmd,'| ctrb rank | %d |\n', op9.ctrb_rank);
fprintf(fmd,'| left geometry err | %.3e |\n', op9.left_error);
fprintf(fmd,'\nCL terminal OFF vs ON (y_max=0.70, soft ON):\n\n');
fprintf(fmd,'| | cover | MAE | soft_rate |\n|---|---:|---:|---:|\n');
fprintf(fmd,'| term OFF | %.4f | %.4f/%.4f | %.3f |\n', r0.joint_cover,r0.MAE(1),r0.MAE(2),r0.soft_rate);
fprintf(fmd,'| term ON | %.4f | %.4f/%.4f | %.3f |\n', r1.joint_cover,r1.MAE(1),r1.MAE(2),r1.soft_rate);
fprintf(fmd,'\nNon-claim: Schur + Lyapunov decrease on **unconstrained free dynamics** does **not** prove recursive feasibility or chance-constrained closed-loop stability.\n\n');

% helpers to aggregate
agg = @(pred) local_agg(RR, pred);

fprintf(fmd,'## P1++ Sigma_eps multi-seed tight limits\n\n');
fprintf(fmd,'| y_max | eps | n | mean cover | std cover | mean MAE1 | mean soft_rate |\n');
fprintf(fmd,'|---:|---|---:|---:|---:|---:|---:|\n');
for ym=ym_list
    for mode={'t2','ml','ols'}
        pref=sprintf('P1pp_eps_%s_y%.2f_',mode{1},ym);
        idx=[];
        for i=1:nJ
            if startsWith(RR{i}.name,pref), idx(end+1)=i; end %#ok<AGROW>
        end
        if isempty(idx), continue; end
        covs=zeros(size(idx)); mae1=covs; srate=covs;
        for k=1:numel(idx)
            covs(k)=RR{idx(k)}.joint_cover;
            mae1(k)=RR{idx(k)}.MAE(1);
            srate(k)=RR{idx(k)}.soft_rate;
        end
        fprintf(fmd,'| %.2f | %s | %d | %.4f | %.4f | %.4f | %.3f |\n', ...
            ym, mode{1}, numel(idx), mean(covs), std(covs), mean(mae1), mean(srate));
    end
end

fprintf(fmd,'\n## P1++ Sigma_obs multi-seed y_max=0.50\n\n');
fprintf(fmd,'| obs | n | mean cover | std cover | mean MAE1 | mean soft_step_viol |\n');
fprintf(fmd,'|---|---:|---:|---:|---:|---:|\n');
for mode={'declared_shape','residual_support','additive'}
    pref=sprintf('P1pp_obs_%s_y0.50_',mode{1});
    idx=[];
    for i=1:nJ
        if startsWith(RR{i}.name,pref), idx(end+1)=i; end %#ok<AGROW>
    end
    if isempty(idx), continue; end
    covs=zeros(size(idx)); mae1=covs; sv=covs;
    for k=1:numel(idx)
        covs(k)=RR{idx(k)}.joint_cover;
        mae1(k)=RR{idx(k)}.MAE(1);
        sv(k)=RR{idx(k)}.soft_step_viol;
    end
    fprintf(fmd,'| %s | %d | %.4f | %.4f | %.4f | %.4f |\n', ...
        mode{1}, numel(idx), mean(covs), std(covs), mean(mae1), mean(sv));
end

fprintf(fmd,'\n## Takeaways\n\n');
fprintf(fmd,'1. Opinion 9: if Ahat is Schur and dlyap residual tiny, free dynamics admit a Lyapunov function — still no RF under constraints/chance/noise.\n');
fprintf(fmd,'2. P1++: multi-seed means/std of cover at y_max=0.50/0.45; use only if cover separates.\n');
fprintf(fmd,'3. Soft ON keeps uncert near zero; rankings remain empirical.\n');
fprintf(fmd,'4. Still open: true RF/stability theorems; original-alpha soft cert; DOF optimality; Cov(o) identity; input-conditional PredVARX.\n');
fclose(fmd);

save(fullfile(resdir,'copyAF_op9_and_p1pp_data.mat'),'op9','r0','r1','RR','plant','-v7.3');
fprintf('\nCOPYAF_DONE\n');
