clc; clear; close all;
% Runner for compare: build ASCII path so CJK folders work in -batch.
repo = ['E:', filesep, 'academic_files', filesep, 'phd-learning', filesep, ...
    char([20195 30721]), filesep, 'PredVARX-MPC'];
an_dir = fullfile(repo,'experiments','copyAN_crte_fixed_surrogate');
am_dir = fullfile(repo,'experiments','copyAM_tracked_cov_only');
assert(exist(an_dir,'dir')==7, ['copyAN missing: ' an_dir]);
assert(exist(am_dir,'dir')==7, ['copyAM missing: ' am_dir]);

fprintf('=== running copyAN (CRTE fixed surrogate) ===\n');
addpath(an_dir); run(fullfile(an_dir,'copyAN_crte_fixed_surrogate.m'));

repo = ['E:', filesep, 'academic_files', filesep, 'phd-learning', filesep, ...
    char([20195 30721]), filesep, 'PredVARX-MPC'];
an_dir = fullfile(repo,'experiments','copyAN_crte_fixed_surrogate');
am_dir = fullfile(repo,'experiments','copyAM_tracked_cov_only');
an_mat = fullfile(an_dir,'results','copyAN_crte_fixed_surrogate_data.mat');
am_mat = fullfile(am_dir,'results','copyAM_tracked_cov_only_data.mat');
assert(exist(an_mat,'file')==2, ['AN mat missing: ' an_mat]);

fprintf('\n=== running copyAM (empirical free dual, no Sigma_n) ===\n');
try, rmpath(an_dir); catch, end
addpath(am_dir); run(fullfile(am_dir,'copyAM_tracked_cov_only.m'));

repo = ['E:', filesep, 'academic_files', filesep, 'phd-learning', filesep, ...
    char([20195 30721]), filesep, 'PredVARX-MPC'];
an_dir = fullfile(repo,'experiments','copyAN_crte_fixed_surrogate');
am_dir = fullfile(repo,'experiments','copyAM_tracked_cov_only');
an_mat = fullfile(an_dir,'results','copyAN_crte_fixed_surrogate_data.mat');
am_mat = fullfile(am_dir,'results','copyAM_tracked_cov_only_data.mat');
assert(exist(am_mat,'file')==2, ['AM mat missing: ' am_mat]);

AN = load(an_mat); AM = load(am_mat);

out_txt = fullfile(an_dir,'results','copyAN_vs_copyAM_compare.txt');
fid = fopen(out_txt,'w');
fprintf(fid,'copyAN vs copyAM same-seed plant comparison\n');
fprintf(fid,'seed_note both runners use rng(20260710,twister)\n');
fprintf(fid,'AN_objective %s\n','fixed CRTE spectral surrogate (free complement)');
fprintf(fid,'AM_objective %s\n','empirical free-coordinate total covariance');
fprintf(fid,'AN_selected_mu %.4f\n', AN.stats.selected_mu);
fprintf(fid,'AN_selected_alpha %.4f\n', AN.stats.selected_alpha);
fprintf(fid,'AN_selected_beta %.4f\n', AN.stats.selected_beta);
fprintf(fid,'AN_selected_validation_nrmse %.6f\n', AN.stats.selected_validation_nrmse);
fprintf(fid,'AN_MAE %.6f %.6f\n', AN.MAE(1), AN.MAE(2));
fprintf(fid,'AM_MAE %.6f %.6f\n', AM.MAE(1), AM.MAE(2));
fprintf(fid,'AN_RMSE %.6f %.6f\n', AN.RMSE(1), AN.RMSE(2));
fprintf(fid,'AM_RMSE %.6f %.6f\n', AM.RMSE(1), AM.RMSE(2));
fprintf(fid,'AN_avg_J %.6f\n', mean(AN.costJ(151:end),'omitnan'));
fprintf(fid,'AM_avg_J %.6f\n', mean(AM.costJ(151:end),'omitnan'));
fprintf(fid,'AN_qp_success %.6f fallbacks %d\n', AN.qp_success_rate, AN.infeasible_count);
fprintf(fid,'AM_qp_success %.6f fallbacks %d\n', AM.qp_success_rate, AM.infeasible_count);
fprintf(fid,'AN_constraint_active_rate %.6f\n', AN.constraint_active_rate);
fprintf(fid,'AM_constraint_active_rate %.6f\n', AM.constraint_active_rate);
fprintf(fid,'AN_cover %.6f\n', AN.cover);
fprintf(fid,'AM_cover %.6f\n', AM.cover);
fprintf(fid,'AN_dual_error %.6e\n', AN.stats.dual_error);
fprintf(fid,'AM_dual_error %.6e\n', AM.stats.dual_error);
fprintf(fid,'note neither result is recursive-feasibility or chance-certificate\n');
fclose(fid);

fprintf('\n=== AN vs AM summary (same plant seed 20260710) ===\n');
fprintf('MAE  AN=[%.4f %.4f]  AM=[%.4f %.4f]\n', AN.MAE(1),AN.MAE(2),AM.MAE(1),AM.MAE(2));
fprintf('RMSE AN=[%.4f %.4f]  AM=[%.4f %.4f]\n', AN.RMSE(1),AN.RMSE(2),AM.RMSE(1),AM.RMSE(2));
fprintf('avg J AN=%.3f  AM=%.3f\n', mean(AN.costJ(151:end),'omitnan'), mean(AM.costJ(151:end),'omitnan'));
fprintf('QP success AN=%.4f AM=%.4f  fallbacks AN=%d AM=%d\n', AN.qp_success_rate, AM.qp_success_rate, AN.infeasible_count, AM.infeasible_count);
fprintf('active AN=%.4f AM=%.4f  cover AN=%.4f AM=%.4f\n', AN.constraint_active_rate, AM.constraint_active_rate, AN.cover, AM.cover);
fprintf('wrote %s\n', out_txt);