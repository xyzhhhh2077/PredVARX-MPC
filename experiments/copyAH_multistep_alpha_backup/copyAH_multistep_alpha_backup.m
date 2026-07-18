%% copyAH_multistep_alpha_backup
% Opinion 8: multi-step original-alpha backup after primary fail.
% Compare none / soft_relaxed / original_alpha_multistep under y_max=0.55.
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
seeds=[94001 94002 94003 94004 94005];
modes={'none','soft_relaxed','original_alpha_multistep'};

jobs={}; %#ok<*SAGROW>
for im=1:numel(modes)
    for s=seeds
        c=base; c.name=sprintf('AH_%s_s%d',modes{im},s);
        c.recovery_mode=modes{im}; c.cl_seed=s; jobs{end+1}=c;
    end
end

fprintf('Smoke...\n');
cs=base; cs.name='smoke'; cs.T_cl=40; cs.y_max=2.0; cs.recovery_mode='original_alpha_multistep';
rs=run_backup_config(cs,plant,data);
fprintf('smoke cover=%.3f orig_cert=%.3f left=%.1e\n',rs.joint_cover,rs.orig_cert_rate,rs.left_err);

nJ=numel(jobs); R=cell(nJ,1);
for i=1:nJ
    fprintf('[%d/%d] %s\n',i,nJ,jobs{i}.name);
    R{i}=run_backup_config(jobs{i},plant,data);
    r=R{i};
    fprintf(['  cover=%.4f MAE=%.3f/%.3f pf=%d orig_cert=%.3f unc=%.3f ' ...
        'orig_rec=%d soft=%d succ|fail=%.3f stages sh/1/feas/Q/prev/um=%d/%d/%d/%d/%d/%d\n'], ...
        r.joint_cover,r.MAE(1),r.MAE(2),r.primary_fail,r.orig_cert_rate,r.uncert_rate, ...
        r.n_orig_rec,r.n_soft_rel,r.orig_success_given_fail, ...
        r.n_short,r.n_one,r.n_feas,r.n_redQ,r.n_prev,r.n_umean);
end

resdir=fullfile(root,'results'); if ~exist(resdir,'dir'), mkdir(resdir); end
fid=fopen(fullfile(resdir,'copyAH_multistep_alpha_backup_metrics.csv'),'w');
fprintf(fid,['name,seed,mode,cover,MAE1,MAE2,pf,orig_cert_rate,uncert_rate,' ...
    'orig_rec,soft_rel,succ_given_fail,soft_also_orig,' ...
    'n_short,n_one,n_feas,n_redQ,n_prev,n_umean\n']);
for i=1:nJ
    r=R{i};
    fprintf(fid,['%s,%d,%s,%.6f,%.6f,%.6f,%d,%.6f,%.6f,' ...
        '%d,%d,%.6f,%.6f,%d,%d,%d,%d,%d,%d\n'], ...
        r.name,r.seed,r.cfg.recovery_mode,r.joint_cover,r.MAE(1),r.MAE(2), ...
        r.primary_fail,r.orig_cert_rate,r.uncert_rate, ...
        r.n_orig_rec,r.n_soft_rel,r.orig_success_given_fail,r.soft_also_orig_rate, ...
        r.n_short,r.n_one,r.n_feas,r.n_redQ,r.n_prev,r.n_umean);
end
fclose(fid);

fmd=fopen(fullfile(resdir,'copyAH_multistep_alpha_backup_report.md'),'w');
fprintf(fmd,'# copyAH multi-step original-alpha backup\n\n');
fprintf(fmd,'Opinion 8 guarantee layer v2. **Not RF/stability.**\n\n');
fprintf(fmd,'Stress: y_max=0.55, T_cl=500, alpha=0.10, 5 seeds.\n\n');
fprintf(fmd,'## Aggregate\n\n');
fprintf(fmd,['| mode | mean cover | mean MAE1 | mean pf | mean orig_cert | ' ...
    'mean uncert | mean succ|fail | mean soft_also_orig |\n']);
fprintf(fmd,'|---|---:|---:|---:|---:|---:|---:|---:|\n');
for im=1:numel(modes)
    mode=modes{im};
    covs=[]; mae=[]; pfv=[]; oc=[]; ur=[]; sg=[]; sa=[];
    for i=1:nJ
        if strcmp(R{i}.cfg.recovery_mode,mode)
            covs(end+1)=R{i}.joint_cover; %#ok<AGROW>
            mae(end+1)=R{i}.MAE(1);
            pfv(end+1)=R{i}.primary_fail;
            oc(end+1)=R{i}.orig_cert_rate;
            ur(end+1)=R{i}.uncert_rate;
            sg(end+1)=R{i}.orig_success_given_fail;
            sa(end+1)=R{i}.soft_also_orig_rate;
        end
    end
    fprintf(fmd,'| %s | %.4f | %.4f | %.1f | %.4f | %.4f | %.4f | %.4f |\n', ...
        mode,mean(covs),mean(mae),mean(pfv),mean(oc),mean(ur),mean(sg),mean(sa));
end
fprintf(fmd,'\n## original_alpha_multistep stage mix\n\n');
fprintf(fmd,'| seed | pf | orig_rec | unc | short | one | feas | redQ | prev | umean | cover | succ|fail |\n');
fprintf(fmd,'|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i=1:nJ
    if ~strcmp(R{i}.cfg.recovery_mode,'original_alpha_multistep'), continue; end
    r=R{i};
    fprintf(fmd,'| %d | %d | %d | %d | %d | %d | %d | %d | %d | %d | %.4f | %.4f |\n', ...
        r.seed,r.primary_fail,r.n_orig_rec,r.n_unc,r.n_short,r.n_one,r.n_feas, ...
        r.n_redQ,r.n_prev,r.n_umean,r.joint_cover,r.orig_success_given_fail);
end
fprintf(fmd,'\n## Takeaways placeholder filled after run\n');
fclose(fmd);

save(fullfile(resdir,'copyAH_multistep_alpha_backup_data.mat'),'R','plant','-v7.3');
fprintf('COPYAH_DONE\n');
