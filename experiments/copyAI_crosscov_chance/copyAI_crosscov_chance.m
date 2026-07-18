%% copyAI_crosscov_chance
% Opinion 6: inject cross-covariance into chance Sigma_y; fair OFF vs ON.
% Does NOT re-prove Boole risk allocation under full cross blocks.
clear; clc; close all;
root = fileparts(mfilename('fullpath'));
addpath(root); addpath(fullfile(root, 'lib'));
rng(20260710, 'twister');

plant = struct('n', 6, 'm', 3, 'p', 30, 'ell', 5, 'tracked', [1 2], ...
    'sw', 0.045, 'se', 0.055, 'noise_cycle', 400, ...
    'sw_min', 0.02, 'sw_max', 0.09, 'se_min', 0.025, 'se_max', 0.10, ...
    'noise_phase_e', pi/3, 'u_min', -3, 'u_max', 3);
A = diag([0.94 0.88 0.78 0.64 0.50 0.35]);
A(1, 2) = 0.10; A(2, 3) = -0.06; A(3, 4) = 0.05; A(4, 5) = 0.04;
B = [0.34 -0.10 0.05; 0.12 0.28 -0.06; 0.05 0.12 0.24; ...
    -0.05 0.06 0.18; 0.02 -0.10 0.14; 0.08 0.02 -0.08];
C = zeros(plant.p, plant.n); C(1, 1) = 1; C(1, 3) = 0.16; C(2, 2) = 1; C(2, 4) = -0.12;
for i = 3:plant.p
    C(i, :) = 0.45 * randn(1, plant.n);
    C(i, :) = C(i, :) / max(norm(C(i, :)), 1e-12);
end
rel = linspace(0.55, 1.65, plant.p)'; Corr = eye(plant.p);
Corr(3, 4) = 0.3; Corr(4, 3) = 0.3; Corr(5, 6) = -0.22; Corr(6, 5) = -0.22;
Corr(8, 9) = 0.18; Corr(9, 8) = 0.18;
Sn = diag(plant.se * rel) * Corr * diag(plant.se * rel);
plant.A = A; plant.B = B; plant.C = C; plant.Sigma_n = Sn; plant.L_n = chol(Sn, 'lower');

T_off = 1500; u_off = 1.2 * randn(plant.m, T_off);
x = zeros(plant.n, 1); y_off = zeros(plant.p, T_off);
for k = 1:T_off
    y_off(:, k) = C * x + plant.L_n * randn(plant.p, 1);
    x = A * x + B * u_off(:, k) + plant.sw * randn(plant.n, 1);
end
data.y_off = y_off; data.u_off = u_off;

base = struct('use_cross_cov', false, 'y_max', 0.55, 'alpha_joint', 0.10, ...
    'T_cl', 500, 'N', 18, 'cl_seed', 1, 'name', '');
seeds = [94001 94002 94003];
modes = [false, true];  % OFF then ON
mode_names = {'cross_OFF', 'cross_ON'};

jobs = {}; %#ok<*SAGROW>
for im = 1:numel(modes)
    for s = seeds
        c = base;
        c.use_cross_cov = modes(im);
        c.cl_seed = s;
        c.name = sprintf('AI_%s_s%d', mode_names{im}, s);
        jobs{end + 1} = c;
    end
end

fprintf('Smoke (short T_cl)...\n');
cs = base; cs.name = 'smoke_OFF'; cs.T_cl = 40; cs.y_max = 2.0; cs.use_cross_cov = false;
rs0 = run_crosscov_config(cs, plant, data);
cs1 = cs; cs1.name = 'smoke_ON'; cs1.use_cross_cov = true;
rs1 = run_crosscov_config(cs1, plant, data);
fprintf('smoke OFF cover=%.3f MAE1=%.3f ||Szo||=%.3e drop=%.3e full=%.3e\n', ...
    rs0.joint_cover, rs0.MAE(1), rs0.Sigma_zo_fro, rs0.drop_cross_rel_err, rs0.full_recon_rel_err);
fprintf('smoke ON  cover=%.3f MAE1=%.3f var_ratio=[%.3f %.3f] qp=%.3f\n', ...
    rs1.joint_cover, rs1.MAE(1), rs1.var_ratio(1), rs1.var_ratio(2), rs1.qp_success);

nJ = numel(jobs); R = cell(nJ, 1);
for i = 1:nJ
    fprintf('[%d/%d] %s\n', i, nJ, jobs{i}.name);
    R{i} = run_crosscov_config(jobs{i}, plant, data);
    r = R{i};
    fprintf(['  cover=%.4f MAE=%.3f/%.3f RMSE=%.3f/%.3f qp=%.3f fail=%d ' ...
        'act=%.3f ||Szo||=%.3e drop=%.3e vr=[%.3f %.3f]\n'], ...
        r.joint_cover, r.MAE(1), r.MAE(2), r.RMSE(1), r.RMSE(2), ...
        r.qp_success, r.qp_fail, r.active_rate, r.Sigma_zo_fro, ...
        r.drop_cross_rel_err, r.var_ratio(1), r.var_ratio(2));
end

resdir = fullfile(root, 'results');
if ~exist(resdir, 'dir'), mkdir(resdir); end

fid = fopen(fullfile(resdir, 'copyAI_crosscov_chance_metrics.csv'), 'w');
fprintf(fid, ['name,seed,use_cross,cover,joint_viol,MAE1,MAE2,RMSE1,RMSE2,' ...
    'qp_success,qp_fail,uncert_rate,active_rate,' ...
    'Sigma_zo_fro,drop_cross_rel_err,full_recon_rel_err,' ...
    'var_full1,var_full2,var_drop1,var_drop2,var_ratio1,var_ratio2,' ...
    'cross_block_fro,left_err,dual_err,cross_flag_rate\n']);
for i = 1:nJ
    r = R{i};
    fprintf(fid, ['%s,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,' ...
        '%.6f,%d,%.6f,%.6f,' ...
        '%.6e,%.6e,%.6e,' ...
        '%.6e,%.6e,%.6e,%.6e,%.6f,%.6f,' ...
        '%.6e,%.6e,%.6e,%.6f\n'], ...
        r.name, r.seed, double(r.use_cross_cov), r.joint_cover, r.joint_viol, ...
        r.MAE(1), r.MAE(2), r.RMSE(1), r.RMSE(2), ...
        r.qp_success, r.qp_fail, r.uncert_rate, r.active_rate, ...
        r.Sigma_zo_fro, r.drop_cross_rel_err, r.full_recon_rel_err, ...
        r.var_full(1), r.var_full(2), r.var_drop(1), r.var_drop(2), ...
        r.var_ratio(1), r.var_ratio(2), ...
        r.cross_block_fro, r.left_err, r.dual_err, r.cross_flag_rate);
end
fclose(fid);

% Aggregate OFF vs ON
fmd = fopen(fullfile(resdir, 'REPORT.md'), 'w');
fprintf(fmd, '# copyAI Opinion 6: cross-covariance in chance Sigma_y\n\n');
fprintf(fmd, '## What was derived\n\n');
fprintf(fmd, 'With `y_c = P z + o`, `z = R'' y_c`, and dual `R'' P = I`:\n\n');
fprintf(fmd, '```text\n');
fprintf(fmd, 'Sigma_y_full = P Sigma_z P'' + P Sigma_zo + Sigma_zo'' P'' + Sigma_o\n');
fprintf(fmd, 'Sigma_y_drop = P Sigma_z P'' + Sigma_o     %% default SMPC path\n');
fprintf(fmd, '```\n\n');
fprintf(fmd, 'Samplewise `R'' o = 0` does **not** imply `Sigma_zo = Cov(z,o) = 0`.\n');
fprintf(fmd, 'Chance direction variance: `hq'' Sigma_y_* hq` compared via `var_ratio = full/drop`.\n\n');
fprintf(fmd, '## What was implemented\n\n');
fprintf(fmd, '- `lib/cross_cov_diagnostics.m`: full/drop recon + optional H variances\n');
fprintf(fmd, '- `lib/sigma_y_chance_blocks.m`: pure block builder for chance rows\n');
fprintf(fmd, '- `centered_smpc_step.m` / `build_chance_rows.m`: `opt.use_cross_cov` (default false)\n');
fprintf(fmd, '- Fair runner: same plant, offline data, seeds; only flag differs\n');
fprintf(fmd, '- Recovery = none (primary QP only) so cover/qp reflect chance law\n\n');
fprintf(fmd, '## What was verified\n\n');
fprintf(fmd, '- Focused test: R''o~0, ||Sigma_zo|| can be >0, full recon rel err ~0\n');
fprintf(fmd, '- Closed-loop OFF vs ON metrics in CSV\n\n');
fprintf(fmd, '## Honest non-claims\n\n');
fprintf(fmd, '- **Does NOT** re-prove Boole / risk allocation under the full cross law.\n');
fprintf(fmd, '- Online still uses offline `Sigma_zo` (not rolling re-estimate of cross).\n');
fprintf(fmd, '- `Sigma_obs` remains a proxy, not necessarily empirical `Sigma_o`.\n');
fprintf(fmd, '- Joint cover is empirical; not a certificate of alpha_joint.\n\n');

fprintf(fmd, '## Aggregate (mean over seeds)\n\n');
fprintf(fmd, ['| mode | cover | MAE1 | MAE2 | RMSE1 | qp_success | active | ' ...
    '||Szo||_F | drop_rel | var_ratio1 | var_ratio2 |\n']);
fprintf(fmd, '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for im = 1:numel(modes)
    covs = []; mae1 = []; mae2 = []; rm1 = []; qp = []; act = [];
    szo = []; dr = []; vr1 = []; vr2 = [];
    for i = 1:nJ
        if logical(R{i}.use_cross_cov) == modes(im)
            covs(end+1) = R{i}.joint_cover; %#ok<AGROW>
            mae1(end+1) = R{i}.MAE(1);
            mae2(end+1) = R{i}.MAE(2);
            rm1(end+1) = R{i}.RMSE(1);
            qp(end+1) = R{i}.qp_success;
            act(end+1) = R{i}.active_rate;
            szo(end+1) = R{i}.Sigma_zo_fro;
            dr(end+1) = R{i}.drop_cross_rel_err;
            vr1(end+1) = R{i}.var_ratio(1);
            vr2(end+1) = R{i}.var_ratio(2);
        end
    end
    fprintf(fmd, '| %s | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f | %.3e | %.3e | %.4f | %.4f |\n', ...
        mode_names{im}, mean(covs), mean(mae1), mean(mae2), mean(rm1), ...
        mean(qp), mean(act), mean(szo), mean(dr), mean(vr1), mean(vr2));
end

fprintf(fmd, '\n## Per-seed rows\n\n');
fprintf(fmd, '| name | cover | MAE1 | qp | active | ||Szo|| | drop | vr1 | vr2 |\n');
fprintf(fmd, '|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i = 1:nJ
    r = R{i};
    fprintf(fmd, '| %s | %.4f | %.4f | %.4f | %.4f | %.3e | %.3e | %.4f | %.4f |\n', ...
        r.name, r.joint_cover, r.MAE(1), r.qp_success, r.active_rate, ...
        r.Sigma_zo_fro, r.drop_cross_rel_err, r.var_ratio(1), r.var_ratio(2));
end

fprintf(fmd, '\n## Setup\n\n');
fprintf(fmd, '- y_max=%.2f, alpha_joint=%.2f, T_cl=%d, N=%d, seeds=%s\n', ...
    base.y_max, base.alpha_joint, base.T_cl, base.N, mat2str(seeds));
fprintf(fmd, '- plant: n=%d p=%d ell=%d tracked=[1 2], split free dual + declared Sigma_n\n', ...
    plant.n, plant.p, plant.ell);
fprintf(fmd, '\nCOPYAI_DONE placeholder filled by real run.\n');
fclose(fmd);

save(fullfile(resdir, 'copyAI_crosscov_chance_data.mat'), 'R', 'plant', '-v7.3');
fprintf('COPYAI_DONE\n');
