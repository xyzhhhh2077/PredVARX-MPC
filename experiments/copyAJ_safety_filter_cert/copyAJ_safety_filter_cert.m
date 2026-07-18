%% copyAJ_safety_filter_cert
% Opinion 8 guarantee-layer upgrade: AFTER primary original-alpha fail,
% switch certified OBJECT to one-step deterministic mean hard-bound safety
% filter (L2-prime level = mean_safety_filter).
%
% Compare recovery modes under tight y_max=0.55:
%   none | soft_relaxed | mean_safety_filter
%
% Certificate naming (hard rule):
%   NEVER call recovery qp_original / original_alpha.
%   L2-prime label is mean_safety_filter (explicit non-chance object).
%
% Non-claims:
%   - Not recursive feasibility
%   - Not original-alpha chance constraint recovery success
%   - Not infinite-horizon joint chance guarantee
clear; clc; close all;
root=fileparts(mfilename('fullpath'));
addpath(root); addpath(fullfile(root,'lib'));
rng(20260710,'twister');

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

base=struct('recovery_mode','none','y_max',0.55,'alpha_joint',0.10, ...
    'T_cl',500,'N',18,'cl_seed',1,'name','');
% >=3 seeds as required
seeds=[95001 95002 95003];
modes={'none','soft_relaxed','mean_safety_filter'};

jobs={}; %#ok<*SAGROW>
for im=1:numel(modes)
    for s=seeds
        c=base; c.name=sprintf('AJ_%s_s%d',modes{im},s);
        c.recovery_mode=modes{im}; c.cl_seed=s; jobs{end+1}=c;
    end
end

fprintf('Smoke mean_safety_filter...\n');
cs=base; cs.name='smoke'; cs.T_cl=40; cs.y_max=2.0; cs.recovery_mode='mean_safety_filter';
rs=run_safety_config(cs,plant,data);
fprintf('smoke cover=%.3f msf_rate=%.3f uncert=%.3f left=%.1e\n', ...
    rs.joint_cover, rs.msf_rate, rs.uncert_rate, rs.left_err);

nJ=numel(jobs); R=cell(nJ,1);
for i=1:nJ
    fprintf('[%d/%d] %s\n',i,nJ,jobs{i}.name);
    R{i}=run_safety_config(jobs{i},plant,data);
    r=R{i};
    fprintf(['  cover=%.4f MAE=%.3f/%.3f pf=%d uncert=%.3f msf=%.3f soft=%.3f ' ...
        'msf_given_fail=%.3f stages qp/prev/um=%d/%d/%d msf_also_orig=%.3f\n'], ...
        r.joint_cover,r.MAE(1),r.MAE(2),r.primary_fail,r.uncert_rate,r.msf_rate, ...
        r.soft_relaxed_rate,r.msf_success_given_fail, ...
        r.n_msf_qp,r.n_msf_prev,r.n_msf_umean,r.msf_also_orig_rate);
end

resdir=fullfile(root,'results'); if ~exist(resdir,'dir'), mkdir(resdir); end
fid=fopen(fullfile(resdir,'copyAJ_safety_filter_cert_metrics.csv'),'w');
fprintf(fid,['name,seed,mode,cover,MAE1,MAE2,pf,uncert_rate,msf_rate,soft_rate,' ...
    'msf_succ_given_fail,soft_succ_given_fail,msf_also_orig,soft_also_mean,' ...
    'n_msf,n_soft,n_unc,n_msf_qp,n_msf_prev,n_msf_umean\n']);
for i=1:nJ
    r=R{i};
    fprintf(fid,['%s,%d,%s,%.6f,%.6f,%.6f,%d,%.6f,%.6f,%.6f,' ...
        '%.6f,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,%d\n'], ...
        r.name,r.seed,r.cfg.recovery_mode,r.joint_cover,r.MAE(1),r.MAE(2), ...
        r.primary_fail,r.uncert_rate,r.msf_rate,r.soft_relaxed_rate, ...
        r.msf_success_given_fail,r.soft_success_given_fail, ...
        r.msf_also_orig_rate,r.soft_also_mean_rate, ...
        r.n_msf,r.n_soft_rel,r.n_unc,r.n_msf_qp,r.n_msf_prev,r.n_msf_umean);
end
fclose(fid);

% Build REPORT.md
fmd=fopen(fullfile(resdir,'REPORT.md'),'w');
fprintf(fmd,'# copyAJ safety-filter certificate (Opinion 8 L2-prime)\n\n');
fprintf(fmd,'Guarantee-layer upgrade after primary alpha fail: **change the certified object**.\n\n');
fprintf(fmd,'Stress: `y_max=0.55`, `T_cl=500`, `alpha_joint=0.10`, seeds=%s.\n\n', mat2str(seeds));
fprintf(fmd,'## Certificate levels (naming hard rule)\n\n');
fprintf(fmd,'| level | meaning |\n|---|---|\n');
fprintf(fmd,'| `qp_primary` | primary N-step chance QP success (original alpha by construction) |\n');
fprintf(fmd,'| `soft_relaxed` | old soft recovery (risk inflate / short N / bound-only) |\n');
fprintf(fmd,'| `mean_safety_filter` | **L2-prime** one-step deterministic mean hard-bound filter |\n');
fprintf(fmd,'| `uncertified_fallback` | sat(u_mean), no certificate |\n\n');
fprintf(fmd,'**Never** labeled as `qp_original` / `original_alpha` for recovery.\n\n');

fprintf(fmd,'## Aggregate by mode\n\n');
fprintf(fmd,['| mode | n | mean cover | mean MAE1 | mean pf | mean uncert | ' ...
    'mean msf_rate | mean soft_rate | mean msf_given_fail | mean msf_also_orig |\n']);
fprintf(fmd,'|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for im=1:numel(modes)
    mode=modes{im};
    covs=[]; mae=[]; pfv=[]; ur=[]; mr=[]; sr=[]; sg=[]; ao=[];
    for i=1:nJ
        if strcmp(R{i}.cfg.recovery_mode,mode)
            covs(end+1)=R{i}.joint_cover; %#ok<AGROW>
            mae(end+1)=R{i}.MAE(1);
            pfv(end+1)=R{i}.primary_fail;
            ur(end+1)=R{i}.uncert_rate;
            mr(end+1)=R{i}.msf_rate;
            sr(end+1)=R{i}.soft_relaxed_rate;
            sg(end+1)=R{i}.msf_success_given_fail;
            ao(end+1)=R{i}.msf_also_orig_rate;
        end
    end
    fprintf(fmd,'| %s | %d | %.4f | %.4f | %.1f | %.4f | %.4f | %.4f | %.4f | %.4f |\n', ...
        mode,numel(covs),mean(covs),mean(mae),mean(pfv),mean(ur),mean(mr),mean(sr),mean(sg),mean(ao));
end

fprintf(fmd,'\n## mean_safety_filter stage mix (per seed)\n\n');
fprintf(fmd,['| seed | pf | n_msf | n_unc | qp | prev | umean | cover | MAE1 | ' ...
    'msf_given_fail | msf_also_orig |\n']);
fprintf(fmd,'|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i=1:nJ
    if ~strcmp(R{i}.cfg.recovery_mode,'mean_safety_filter'), continue; end
    r=R{i};
    fprintf(fmd,'| %d | %d | %d | %d | %d | %d | %d | %.4f | %.4f | %.4f | %.4f |\n', ...
        r.seed,r.primary_fail,r.n_msf,r.n_unc,r.n_msf_qp,r.n_msf_prev,r.n_msf_umean, ...
        r.joint_cover,r.MAE(1),r.msf_success_given_fail,r.msf_also_orig_rate);
end

fprintf(fmd,'\n## L2-prime object definition\n\n');
fprintf(fmd,'After primary fail, solve for first-step `uk` in input bounds s.t.\n\n');
fprintf(fmd,'```text\nH * y_mean_pred(uk)  <=  h     (and lower-side mirror)\n');
fprintf(fmd,'y_mean_pred(uk) = y_bar + P A z + P B (uk - u_bar)\n```\n\n');
fprintf(fmd,'Optional backup: accept `u_prev` / `u_mean` only if they pass the same mean rows.\n\n');

fprintf(fmd,'## Honest non-claims\n\n');
fprintf(fmd,'1. **Not recursive feasibility.**\n');
fprintf(fmd,'2. **Not original-alpha chance-constraint recovery.** `msf_also_orig` is diagnostic only.\n');
fprintf(fmd,'3. **Not** a multi-step joint chance guarantee.\n');
fprintf(fmd,'4. Mean hard bound is a **different certified object** (deterministic one-step mean).\n');
fprintf(fmd,'5. Soft remains an engineering relaxed certificate; MSF is L2-prime mean-safety.\n');
fclose(fmd);

save(fullfile(resdir,'copyAJ_safety_filter_cert_data.mat'),'R','plant','-v7.3');
fprintf('COPYAJ_DONE\n');
