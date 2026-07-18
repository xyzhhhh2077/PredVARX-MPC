%% copyAK_terminal_set_rf
% Opinion 9: terminal set / recursive-feasibility NUMERICAL probes (not theorems).
% - Documents minimal RF condition chain (C1-C6) with proved/unproved flags
% - Optional terminal cost Pterm~dlyap(A',Qf) (default OFF)
% - Optional coarse terminal ellipsoid probe via spectral-box linear rows (default OFF)
% - OFF vs ON closed-loop compare at medium/tight y_max, multi-seed
% Does NOT modify main/ or other copies.
clear; clc; close all;
root = fileparts(mfilename('fullpath'));
addpath(root); addpath(fullfile(root,'lib'));
rng(20260710,'twister');

%% Plant (aligned with copyAF lineage)
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
patchdir=fullfile(root,'patches'); if ~exist(patchdir,'dir'), mkdir(patchdir); end

%% ========== Opinion 9 RF / Lyapunov probe ==========
fprintf('=== Opinion 9 RF condition chain + Lyapunov probe ===\n');
op9 = opinion9_rf_probe(y_off, u_off, plant.ell, plant.tracked, Sn);
fprintf('rho(A)=%.6f Schur=%d dlyap_ok=%d lyap_res=%.3e Vdec=%.3f ctrb=%d\n', ...
    op9.rho_A, op9.is_schur, op9.dlyap_ok, op9.lyap_residual_fro, ...
    op9.V_decrease_ok_ratio, op9.ctrb_rank);
fprintf('C1=%d C2=%d C3=%d C4=%d C5=%d C6=%d (1=probed-pass / 0=unproved-or-fail)\n', ...
    op9.C1_schur, op9.C2_dlyap, op9.C3_Vdec_free, ...
    op9.C4_Xf_invar_constrained, op9.C5_stage_to_Xf_theorem, op9.C6_chance_RF);
if isfield(op9,'alpha_cal') && isfield(op9.alpha_cal,'alpha_recommend')
    fprintf('alpha_recommend (unit-ball V q0.9)=%.4g\n', op9.alpha_cal.alpha_recommend);
end

% Calibrate alpha_term on offline latent energy scale
alpha_term = NaN;
if op9.dlyap_ok
    Zc = op9.R' * (y_off - mean(y_off,2));  % rough; runner uses stats.y_mean
    % re-id mean via probe model: use R from op9
    [~,~,~,~,~,st0] = split_control_free_ivr_varx(y_off,u_off,plant.ell,plant.tracked,Sn);
    Zc = op9.R' * (y_off - st0.y_mean);
    Voff = sum((op9.Pterm * Zc) .* Zc, 1);
    v_med = median(Voff); v_p90 = quantile(Voff, 0.90);
    a_unit = op9.alpha_cal.alpha_recommend;
    % Practical CL level: mix unit direction quantile with offline energy
    alpha_term = max(a_unit, 0.75*v_p90);
    fprintf('alpha_term calib: unit_q=%.4g Voff_med=%.4g Voff_p90=%.4g -> alpha=%.4g\n', ...
        a_unit, v_med, v_p90, alpha_term);
end

%% ========== CL jobs: terminal OFF vs ON ==========
base=struct('sigma_eps_mode','t2','use_cross_cov',false,'sigma_obs_mode','declared_shape', ...
    'enable_soft_recovery',true,'use_terminal_cost',false,'use_terminal_set',false, ...
    'alpha_term',alpha_term,'terminal_box_scale',1.0, ...
    'input_residualize',false,'y_max',0.70,'alpha_joint',0.10, ...
    'T_cl',500,'N',18,'cl_seed',91001,'name','');

seeds_med = [91001 91002];
seeds_tgt = [92001];
ym_med = 0.70;
ym_tgt = 0.50;

jobs = {}; %#ok<*SAGROW>
% Medium y_max: OFF vs COST vs SET, 2 seeds
for s = seeds_med
    c=base; c.name=sprintf('MED_OFF_y%.2f_s%d',ym_med,s); c.y_max=ym_med; c.cl_seed=s;
    c.use_terminal_cost=false; c.use_terminal_set=false; c.T_cl=400; jobs{end+1}=c;
    c=base; c.name=sprintf('MED_COST_y%.2f_s%d',ym_med,s); c.y_max=ym_med; c.cl_seed=s;
    c.use_terminal_cost=true; c.use_terminal_set=false; c.T_cl=400; jobs{end+1}=c;
    c=base; c.name=sprintf('MED_SET_y%.2f_s%d',ym_med,s); c.y_max=ym_med; c.cl_seed=s;
    c.use_terminal_cost=true; c.use_terminal_set=true; c.alpha_term=alpha_term; c.T_cl=400; jobs{end+1}=c;
end
% Tight y_max: OFF vs COST only (1 seed) — SET optional probe
for s = seeds_tgt
    c=base; c.name=sprintf('TGT_OFF_y%.2f_s%d',ym_tgt,s); c.y_max=ym_tgt; c.cl_seed=s;
    c.use_terminal_cost=false; c.use_terminal_set=false; c.T_cl=400; jobs{end+1}=c;
    c=base; c.name=sprintf('TGT_COST_y%.2f_s%d',ym_tgt,s); c.y_max=ym_tgt; c.cl_seed=s;
    c.use_terminal_cost=true; c.use_terminal_set=false; c.T_cl=400; jobs{end+1}=c;
    c=base; c.name=sprintf('TGT_SET_y%.2f_s%d',ym_tgt,s); c.y_max=ym_tgt; c.cl_seed=s;
    c.use_terminal_cost=true; c.use_terminal_set=true; c.alpha_term=alpha_term; c.T_cl=400; jobs{end+1}=c;
end

nJ=numel(jobs); RR=cell(nJ,1);
for i=1:nJ
    fprintf('[%d/%d] %s\n',i,nJ,jobs{i}.name);
    RR{i}=run_seeded_config(jobs{i},plant,data);
    r=RR{i};
    fprintf('  cover=%.4f MAE=%.3f/%.3f qp=%.3f soft=%.3f unc=%.3f Vmed=%.3g\n', ...
        r.joint_cover,r.MAE(1),r.MAE(2),r.qp_rate,r.soft_rate,r.uncert_rate,r.V_term_median);
end

%% CSV
fid=fopen(fullfile(resdir,'copyAK_metrics.csv'),'w');
fprintf(fid,['name,seed,y_max,term_cost,term_set,alpha_term,rho_A,lyap_res,Vdec,' ...
    'cover,viol,MAE1,MAE2,qp_rate,soft_rate,uncert_rate,' ...
    'primary_fail,soft_ok,n_risk,n_short,n_bound,soft_step_viol,' ...
    'V_term_med,V_term_p90,V_term_max\n']);
for i=1:nJ
    r=RR{i}; c=r.cfg;
    at = c.alpha_term; if ~isfinite(at), at=-1; end
    fprintf(fid,['%s,%d,%.3f,%d,%d,%.6g,%.6f,%.3e,%.4f,' ...
        '%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,' ...
        '%d,%d,%d,%d,%d,%.6f,%.6g,%.6g,%.6g\n'], ...
        r.name,r.seed,c.y_max,c.use_terminal_cost,c.use_terminal_set,at, ...
        r.rho_A,r.lyap_residual_fro,r.V_decrease_ok_ratio, ...
        r.joint_cover,r.joint_viol,r.MAE(1),r.MAE(2),r.qp_rate,r.soft_rate,r.uncert_rate, ...
        r.primary_fail,r.soft_ok,r.n_risk,r.n_short,r.n_bound, ...
        local_nan(r.soft_step_viol), local_nan(r.V_term_median), ...
        local_nan(r.V_term_p90), local_nan(r.V_term_max));
end
fclose(fid);

%% Markdown report
fmd=fopen(fullfile(resdir,'copyAK_terminal_set_rf_report.md'),'w');
fprintf(fmd,'# copyAK: Opinion 9 terminal set / RF numerical probes\n\n');
fprintf(fmd,'**Not theorem closures. No chance-constraint RF. No closed-loop stability theorem.**\n\n');

fprintf(fmd,'## Minimal RF condition chain\n\n');
fprintf(fmd,'| ID | Condition | Status in this probe |\n|---|---|---|\n');
fprintf(fmd,'| C1 | rho(Ahat)<1 (Schur latent model) | %s (rho=%.6f) |\n', ...
    tf(op9.C1_schur), op9.rho_A);
fprintf(fmd,'| C2 | Pterm=dlyap(A'',Qf), residual small | %s (res=%.3e) |\n', ...
    tf(op9.C2_dlyap), op9.lyap_residual_fro);
fprintf(fmd,'| C3 | V-decrease on free z+=A z | %s (ratio=%.3f) |\n', ...
    tf(op9.C3_Vdec_free), op9.V_decrease_ok_ratio);
fprintf(fmd,'| C4 | terminal law keeps Xf +invariant under constraints | **UNPROVED** |\n');
fprintf(fmd,'| C5 | stage feasible => z_N in exact Xf | **UNPROVED** (optional box rows only) |\n');
fprintf(fmd,'| C6 | chance-constraint RF under noise | **UNPROVED** |\n\n');

fprintf(fmd,'### What is implemented (engineering probe)\n\n');
fprintf(fmd,'1. Terminal cost (opt, default OFF): `Pterm≈dlyap(A'',Qf)`, add `z_N''Pterm z_N`.\n');
fprintf(fmd,'2. Terminal set (opt, default OFF): spectral-box outer approx of\n');
fprintf(fmd,'   `z_N'' Pterm z_N ≤ α_term` as linear QP rows (not exact ellipsoid SOCP).\n');
fprintf(fmd,'3. `α_term` calibrated from free unit-ball V quantile + offline latent energy.\n');
fprintf(fmd,'4. Soft recovery remains engineering ladder (not original-α cert).\n\n');

fprintf(fmd,'### alpha_term used\n\n');
fprintf(fmd,'`alpha_term = %.6g`\n\n', alpha_term);

fprintf(fmd,'## Offline Lyapunov probe\n\n');
fprintf(fmd,'| quantity | value |\n|---|---:|\n');
fprintf(fmd,'| rho(Ahat) | %.6f |\n', op9.rho_A);
fprintf(fmd,'| is_schur | %d |\n', op9.is_schur);
fprintf(fmd,'| dlyap_ok | %d |\n', op9.dlyap_ok);
fprintf(fmd,'| lyap residual F | %.3e |\n', op9.lyap_residual_fro);
fprintf(fmd,'| V-decrease OK ratio | %.3f |\n', op9.V_decrease_ok_ratio);
fprintf(fmd,'| ctrb rank | %d |\n', op9.ctrb_rank);
fprintf(fmd,'| left geometry err | %.3e |\n\n', op9.left_error);

fprintf(fmd,'## Closed-loop OFF vs ON (multi-seed)\n\n');
fprintf(fmd,'Modes: OFF = no terminal; COST = terminal cost only; SET = cost + box terminal set.\n\n');

groups = { ...
    sprintf('MED_*_y%.2f_',ym_med), ym_med; ...
    sprintf('TGT_*_y%.2f_',ym_tgt), ym_tgt};
% Better: aggregate by prefix family
fams = {'MED_OFF','MED_COST','MED_SET','TGT_OFF','TGT_COST','TGT_SET'};
fprintf(fmd,['| family | n | mean cover | mean MAE1 | mean MAE2 | mean qp | mean soft | mean unc | ' ...
    'mean V_term_med |\n']);
fprintf(fmd,'|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for fi=1:numel(fams)
    pref = fams{fi};
    idx=[];
    for i=1:nJ
        if startsWith(RR{i}.name, pref), idx(end+1)=i; end %#ok<AGROW>
    end
    if isempty(idx), continue; end
    covs=zeros(size(idx)); m1=covs; m2=covs; qp=covs; sf=covs; un=covs; vm=covs;
    for k=1:numel(idx)
        rr=RR{idx(k)};
        covs(k)=rr.joint_cover; m1(k)=rr.MAE(1); m2(k)=rr.MAE(2);
        qp(k)=rr.qp_rate; sf(k)=rr.soft_rate; un(k)=rr.uncert_rate;
        vm(k)=rr.V_term_median;
    end
    fprintf(fmd,'| %s | %d | %.4f | %.4f | %.4f | %.3f | %.3f | %.3f | %.4g |\n', ...
        pref, numel(idx), mean(covs), mean(m1), mean(m2), mean(qp), mean(sf), mean(un), mean(vm,'omitnan'));
end

fprintf(fmd,'\n### Per-run detail\n\n');
fprintf(fmd,['| name | cover | MAE1 | MAE2 | qp | soft | unc | lyap_res | Vdec | V_term_med |\n']);
fprintf(fmd,'|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i=1:nJ
    r=RR{i};
    fprintf(fmd,'| %s | %.4f | %.4f | %.4f | %.3f | %.3f | %.3f | %.2e | %.3f | %.4g |\n', ...
        r.name, r.joint_cover, r.MAE(1), r.MAE(2), r.qp_rate, r.soft_rate, r.uncert_rate, ...
        r.lyap_residual_fro, r.V_decrease_ok_ratio, r.V_term_median);
end

fprintf(fmd,'\n## Non-claims (mandatory)\n\n');
fprintf(fmd,'1. No recursive feasibility theorem under chance constraints / noise.\n');
fprintf(fmd,'2. No closed-loop stability theorem (with or without terminal cost/set).\n');
fprintf(fmd,'3. Free-dynamics Schur + V-decrease ≠ constrained Xf positive invariance.\n');
fprintf(fmd,'4. Spectral-box rows are a **coarse outer approximation**, not exact ellipsoid Xf.\n');
fprintf(fmd,'5. Soft recovery ≠ original joint-risk certificate.\n');
fprintf(fmd,'6. Empirical cover/MAE shifts under terminal ON are probes, not RF proofs.\n\n');

fprintf(fmd,'## Takeaways\n\n');
fprintf(fmd,'1. If C1-C3 pass, free latent dynamics admit a discrete Lyapunov function (numeric).\n');
fprintf(fmd,'2. Terminal cost/set are optional regularizers/constraints default OFF.\n');
fprintf(fmd,'3. Compare OFF/COST/SET at medium and tight y_max for empirical effect size.\n');
fprintf(fmd,'4. Still open: true Xf construction under chance rows; kappa terminal law; CL stability.\n');
fclose(fmd);

save(fullfile(resdir,'copyAK_terminal_set_rf_data.mat'),'op9','RR','plant','alpha_term','-v7.3');
fprintf('\nCOPYAK_DONE\n');

function s = tf(b)
if b, s='PROBED-PASS'; else, s='FAIL/UNPROVED'; end
end

function x = local_nan(v)
if isempty(v) || (isnumeric(v) && ~isfinite(v)), x = -1; else, x = v; end
end
