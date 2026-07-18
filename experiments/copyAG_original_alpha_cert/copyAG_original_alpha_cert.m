%% copyAG_original_alpha_cert
% Opinion 8 guarantee-layer experiment:
% compare recovery modes under tight y_max:
%   none | soft_relaxed (old, not original alpha) | original_alpha (recheck)
% Does NOT modify earlier copies.
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

base=struct('sigma_eps_mode','t2','recovery_mode','none', ...
    'y_max',0.55,'alpha_joint',0.10,'T_cl',500,'N',18,'cl_seed',1,'name','');
seeds=[93001 93002 93003 93004 93005];
modes={'none','soft_relaxed','original_alpha'};

jobs={}; %#ok<*SAGROW>
for im=1:numel(modes)
    for s=seeds
        c=base; c.name=sprintf('AG_%s_s%d',modes{im},s);
        c.recovery_mode=modes{im}; c.cl_seed=s; jobs{end+1}=c;
    end
end

fprintf('Smoke cert helper...\n');
cfg0=base; cfg0.name='smoke'; cfg0.cl_seed=1; cfg0.T_cl=30; cfg0.y_max=2.0;
cfg0.recovery_mode='original_alpha';
rs=run_cert_config(cfg0,plant,data);
fprintf('smoke cover=%.3f orig_cert_rate=%.3f left=%.1e\n', ...
    rs.joint_cover, rs.orig_cert_rate, rs.left_err);

nJ=numel(jobs); R=cell(nJ,1);
for i=1:nJ
    fprintf('[%d/%d] %s\n',i,nJ,jobs{i}.name);
    R{i}=run_cert_config(jobs{i},plant,data);
    r=R{i};
    fprintf(['  cover=%.4f MAE=%.3f/%.3f pf=%d orig_cert=%.3f ' ...
        'soft_rel=%d orig_rec=%d unc=%d stages sh/1s/prev/um=%d/%d/%d/%d\n'], ...
        r.joint_cover,r.MAE(1),r.MAE(2),r.primary_fail,r.orig_cert_rate, ...
        r.n_soft_rel,r.n_orig_rec,r.n_unc,r.n_short,r.n_onestep,r.n_prev,r.n_umean);
end

resdir=fullfile(root,'results'); if ~exist(resdir,'dir'), mkdir(resdir); end
fid=fopen(fullfile(resdir,'copyAG_original_alpha_cert_metrics.csv'),'w');
fprintf(fid,['name,seed,mode,y_max,cover,viol,MAE1,MAE2,primary_fail,' ...
    'orig_cert_rate,uncert_rate,soft_rel_rate,orig_rec_rate,' ...
    'n_soft_rel,n_orig_rec,n_unc,n_short,n_onestep,n_prev,n_umean,' ...
    'orig_success_given_fail,soft_also_orig_rate\n']);
for i=1:nJ
    r=R{i}; c=r.cfg;
    fprintf(fid,['%s,%d,%s,%.3f,%.6f,%.6f,%.6f,%.6f,%d,' ...
        '%.6f,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,%d,%d,%.6f,%.6f\n'], ...
        r.name,r.seed,c.recovery_mode,c.y_max,r.joint_cover,r.joint_viol, ...
        r.MAE(1),r.MAE(2),r.primary_fail, ...
        r.orig_cert_rate,r.uncert_rate,r.soft_relaxed_rate,r.orig_rec_rate, ...
        r.n_soft_rel,r.n_orig_rec,r.n_unc,r.n_short,r.n_onestep,r.n_prev, ...
        r.n_umean,r.orig_success_given_fail,r.soft_also_orig_rate);
end
fclose(fid);

fmd=fopen(fullfile(resdir,'copyAG_original_alpha_cert_report.md'),'w');
fprintf(fmd,'# copyAG original-alpha certificate experiment\n\n');
fprintf(fmd,'Opinion 8 guarantee layer. **Not recursive feasibility / stability.**\n\n');
fprintf(fmd,'Tight stress: `y_max=0.55`, `T_cl=500`, `alpha_joint=0.10`, 5 seeds.\n\n');
fprintf(fmd,'Recovery modes:\n\n');
fprintf(fmd,'1. `none` — uncertified sat(u_mean)\n');
fprintf(fmd,'2. `soft_relaxed` — old risk-inflate/short-N/bound-only (NOT original alpha)\n');
fprintf(fmd,'3. `original_alpha` — only accept inputs that re-pass original-alpha constant-hold rows\n\n');

fprintf(fmd,'## Aggregate by mode\n\n');
fprintf(fmd,['| mode | n | mean cover | mean MAE1 | mean pf | ' ...
    'mean orig_cert_rate | mean uncert_rate | mean orig_succ|fail | mean soft_also_orig |\n']);
fprintf(fmd,'|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for im=1:numel(modes)
    mode=modes{im}; idx=[];
    for i=1:nJ
        if strcmp(R{i}.cfg.recovery_mode,mode), idx(end+1)=i; end %#ok<AGROW>
    end
    covs=zeros(size(idx)); mae=covs; pfv=covs; oc=covs; ur=covs; sg=covs; sa=covs;
    for k=1:numel(idx)
        covs(k)=R{idx(k)}.joint_cover;
        mae(k)=R{idx(k)}.MAE(1);
        pfv(k)=R{idx(k)}.primary_fail;
        oc(k)=R{idx(k)}.orig_cert_rate;
        ur(k)=R{idx(k)}.uncert_rate;
        sg(k)=R{idx(k)}.orig_success_given_fail;
        sa(k)=R{idx(k)}.soft_also_orig_rate;
    end
    fprintf(fmd,'| %s | %d | %.4f | %.4f | %.1f | %.4f | %.4f | %.4f | %.4f |\n', ...
        mode, numel(idx), mean(covs), mean(mae), mean(pfv), mean(oc), mean(ur), mean(sg), mean(sa));
end

fprintf(fmd,'\n## Per-seed original_alpha stage mix\n\n');
fprintf(fmd,['| seed | pf | orig_rec | unc | shortN | one-step | prev | umean | ' ...
    'cover | orig_cert_rate |\n']);
fprintf(fmd,'|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i=1:nJ
    if ~strcmp(R{i}.cfg.recovery_mode,'original_alpha'), continue; end
    r=R{i};
    fprintf(fmd,'| %d | %d | %d | %d | %d | %d | %d | %d | %.4f | %.4f |\n', ...
        r.seed,r.primary_fail,r.n_orig_rec,r.n_unc,r.n_short,r.n_onestep, ...
        r.n_prev,r.n_umean,r.joint_cover,r.orig_cert_rate);
end

fprintf(fmd,'\n## Certificate definition\n\n');
fprintf(fmd,'Applied first-step `uk` is original-alpha certified iff\n\n');
fprintf(fmd,'```text\nmax(A_ch(alpha) * repmat(uk,N,1) - b_ch(alpha)) <= tol\n```\n\n');
fprintf(fmd,'with the **same** `alpha_joint` used by the primary controller.\n');
fprintf(fmd,'Short-horizon solves are only candidate generators; acceptance requires recheck.\n\n');
fprintf(fmd,'## Non-claims\n\n');
fprintf(fmd,'1. Not recursive feasibility.\n');
fprintf(fmd,'2. Not infinite-horizon joint chance guarantee.\n');
fprintf(fmd,'3. soft_relaxed success does **not** imply original-alpha certification.\n');
fprintf(fmd,'4. If original_alpha still has uncertified steps, the guarantee layer is incomplete.\n');
fclose(fmd);

save(fullfile(resdir,'copyAG_original_alpha_cert_data.mat'),'R','plant','-v7.3');
fprintf('COPYAG_DONE\n');
