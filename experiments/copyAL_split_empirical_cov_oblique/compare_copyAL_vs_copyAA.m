%% compare_copyAL_vs_copyAA
% Same-seed plant compare of copyAL (empirical free dual, no Sigma_n)
% against copyAA (declared Sigma_n free dual). Does not mutate either copy.
% Paths are rebuilt with char codes so runner clear/CJK path is robust.
clear; clc; close all;

repo = ['E:', filesep, 'academic_files', filesep, 'phd-learning', filesep, ...
    char([20195 30721]), filesep, 'PredVARX-MPC'];
al_dir = fullfile(repo,'experiments','copyAL_split_empirical_cov_oblique');
aa_dir = fullfile(repo,'experiments','copyAA_split_control_free_oblique');
assert(exist(al_dir,'dir')==7, ['copyAL missing: ' al_dir]);
assert(exist(aa_dir,'dir')==7, ['copyAA missing: ' aa_dir]);

fprintf('=== Running copyAL (empirical free dual) ===\n');
addpath(al_dir);
run(fullfile(al_dir,'copyAL_split_empirical_cov_oblique.m'));
% rebuild absolute paths after runner clear
repo = ['E:', filesep, 'academic_files', filesep, 'phd-learning', filesep, ...
    char([20195 30721]), filesep, 'PredVARX-MPC'];
al_dir = fullfile(repo,'experiments','copyAL_split_empirical_cov_oblique');
aa_dir = fullfile(repo,'experiments','copyAA_split_control_free_oblique');
al_mat_path = fullfile(al_dir,'results','copyAL_split_empirical_cov_oblique_data.mat');
al_metrics_path = fullfile(al_dir,'results','copyAL_split_empirical_cov_oblique_metrics.txt');
assert(exist(al_mat_path,'file')==2, ['AL mat missing: ' al_mat_path]);

fprintf('\n=== Running copyAA (declared Sigma_n free dual) ===\n');
try, rmpath(al_dir); catch, end
addpath(aa_dir);
run(fullfile(aa_dir,'copyAA_split_control_free_oblique.m'));
repo = ['E:', filesep, 'academic_files', filesep, 'phd-learning', filesep, ...
    char([20195 30721]), filesep, 'PredVARX-MPC'];
al_dir = fullfile(repo,'experiments','copyAL_split_empirical_cov_oblique');
aa_dir = fullfile(repo,'experiments','copyAA_split_control_free_oblique');
aa_mat_path = fullfile(aa_dir,'results','copyAA_split_control_free_oblique_data.mat');
aa_metrics_path = fullfile(aa_dir,'results','copyAA_split_control_free_oblique_metrics.txt');
al_mat_path = fullfile(al_dir,'results','copyAL_split_empirical_cov_oblique_data.mat');
al_metrics_path = fullfile(al_dir,'results','copyAL_split_empirical_cov_oblique_metrics.txt');
assert(exist(aa_mat_path,'file')==2, ['AA mat missing: ' aa_mat_path]);
assert(exist(al_mat_path,'file')==2, ['AL mat missing after AA: ' al_mat_path]);

AL = load(al_mat_path);
AA = load(aa_mat_path);

out_txt = fullfile(al_dir,'results','copyAL_vs_copyAA_compare.txt');
fid = fopen(out_txt,'w');
fprintf(fid,'copyAL vs copyAA same-seed plant comparison\n');
fprintf(fid,'seed_note both runners use rng(20260710,twister) at start\n');
fprintf(fid,'AL_metric %s\n', AL.stats.metric_name);
fprintf(fid,'AL_uses_true_Sigma_n %d\n', AL.stats.uses_true_Sigma_n);
fprintf(fid,'AA_uses_declared_Sigma_n 1\n');
fprintf(fid,'AL_MAE %.12f %.12f\n', AL.MAE(1), AL.MAE(2));
fprintf(fid,'AA_MAE %.12f %.12f\n', AA.MAE(1), AA.MAE(2));
fprintf(fid,'AL_RMSE %.12f %.12f\n', AL.RMSE(1), AL.RMSE(2));
fprintf(fid,'AA_RMSE %.12f %.12f\n', AA.RMSE(1), AA.RMSE(2));
fprintf(fid,'AL_qp_success_rate %.12f\n', AL.qp_success_rate);
fprintf(fid,'AA_qp_success_rate %.12f\n', AA.qp_success_rate);
fprintf(fid,'AL_infeasible_count %d\n', AL.infeasible_count);
fprintf(fid,'AA_infeasible_count %d\n', AA.infeasible_count);
fprintf(fid,'AL_constraint_active_rate %.12f\n', AL.constraint_active_rate);
fprintf(fid,'AA_constraint_active_rate %.12f\n', AA.constraint_active_rate);
fprintf(fid,'AL_cover %.12f\n', AL.cover);
fprintf(fid,'AA_cover %.12f\n', mean(AA.exitflag(151:end)>0));
fprintf(fid,'AL_free_emp_cov_objective %.12e\n', AL.stats.free_emp_cov_objective);
fprintf(fid,'AL_free_emp_cov_improvement %.12e\n', AL.stats.free_emp_cov_improvement);
fprintf(fid,'AA_free_noise_objective %.12e\n', AA.stats.free_noise_objective);
fprintf(fid,'AA_free_noise_improvement %.12e\n', AA.stats.free_noise_improvement);
fprintf(fid,'AL_free_oblique_norm %.12e\n', AL.stats.free_oblique_norm);
fprintf(fid,'AA_free_oblique_norm %.12e\n', AA.stats.free_oblique_norm);
fprintf(fid,'AL_dual_error %.12e\n', AL.stats.dual_error);
fprintf(fid,'AA_dual_error %.12e\n', AA.stats.dual_error);
fprintf(fid,'note cover is warm-horizon QP-success time-fraction proxy only\n');
fprintf(fid,'note AL free metric is empirical total free-cov; AA is declared Sigma_n\n');
fprintf(fid,'note neither result is PredVAR min-innovation covariance theorem\n');
fclose(fid);

fprintf('\n=== AL vs AA summary ===\n');
fprintf('MAE  AL=[%.4f %.4f]  AA=[%.4f %.4f]\n', AL.MAE(1),AL.MAE(2),AA.MAE(1),AA.MAE(2));
fprintf('RMSE AL=[%.4f %.4f]  AA=[%.4f %.4f]\n', AL.RMSE(1),AL.RMSE(2),AA.RMSE(1),AA.RMSE(2));
fprintf('QP success AL=%.4f AA=%.4f  fallbacks AL=%d AA=%d\n', ...
    AL.qp_success_rate, AA.qp_success_rate, AL.infeasible_count, AA.infeasible_count);
fprintf('active AL=%.4f AA=%.4f  cover AL=%.4f AA=%.4f\n', ...
    AL.constraint_active_rate, AA.constraint_active_rate, AL.cover, mean(AA.exitflag(151:end)>0));
fprintf('free dual ||R-P|| AL=%.3e AA=%.3e\n', AL.stats.free_oblique_norm, AA.stats.free_oblique_norm);
fprintf('AL free emp-cov obj=%.6e gain=%.6e\n', AL.stats.free_emp_cov_objective, AL.stats.free_emp_cov_improvement);
fprintf('AA free noise  obj=%.6e gain=%.6e\n', AA.stats.free_noise_objective, AA.stats.free_noise_improvement);
fprintf('wrote %s\n', out_txt);
fprintf('AL metrics: %s\n', al_metrics_path);
fprintf('AA metrics: %s\n', aa_metrics_path);
