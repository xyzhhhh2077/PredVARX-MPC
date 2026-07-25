function aggregate = aggregate_copyAP_results(results_dir, output_file, require_complete)
%AGGREGATE_COPYAP_RESULTS Validate and aggregate per-run copyAP MAT files.
% The aggregate contains one compact record per run plus dense 4-by-20 metric
% arrays. It is always saved with -v7.3.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(results_dir), results_dir = fullfile(here,'results','runs'); end
if nargin < 2 || isempty(output_file)
    output_file = fullfile(here,'results','copyAP_crte_multistep_task_20x4_aggregate.mat');
end
if nargin < 3 || isempty(require_complete), require_complete = true; end
Hset = [1 3 6 18]; seedset = 1:20;
files = dir(fullfile(results_dir,'copyAP_run_H*_seed*.mat'));

shape = [numel(Hset),numel(seedset)];
summary = struct();
summary.present = false(shape);
summary.runtime_seconds = nan(shape);
summary.teacher_objective = nan(shape);
summary.teacher_prediction_term = nan(shape);
summary.teacher_task_term = nan(shape);
summary.teacher_noise_term = nan(shape);
summary.teacher_reach_min = nan(shape);
summary.validation_nrmse = nan(shape);
summary.qp_success_rate = nan(shape);
summary.infeasible_count = nan(shape);
summary.constraint_active_rate = nan(shape);
summary.cost_mean = nan(shape);
summary.cost_sum = nan(shape);
summary.MAE = nan([2 shape]); summary.RMSE = nan([2 shape]);
summary.Bias = nan([2 shape]); summary.upper_violation_rate = nan([2 shape]);
summary.abs_violation_rate = nan([2 shape]);
records = struct([]); seen = strings(0,1);
for k = 1:numel(files)
    f = fullfile(files(k).folder,files(k).name);
    s = load(f,'run');
    assert(isfield(s,'run') && s.run.completed,'Incomplete copyAP result: %s',f);
    r = s.run; H = r.config.task_horizon; seed_id = r.seed.seed_id;
    ih = find(Hset==H,1); is = find(seedset==seed_id,1);
    assert(~isempty(ih) && ~isempty(is),'Unexpected H/seed in %s',f);
    key = sprintf('H%02d_seed%02d',H,seed_id);
    assert(~any(seen==key),'Duplicate run key %s',key); seen(end+1,1)=key; %#ok<AGROW>
    assert(~r.algorithm_contract.uses_true_Sigma_n,'Unknown-noise contract failed in %s',f);
    assert(r.teacher.task_horizon==H,'Teacher/config task horizon mismatch in %s',f);
    summary.present(ih,is)=true;
    summary.runtime_seconds(ih,is)=r.runtime_seconds;
    summary.teacher_objective(ih,is)=r.teacher.objective;
    summary.teacher_prediction_term(ih,is)=r.teacher.prediction_term;
    summary.teacher_task_term(ih,is)=r.teacher.task_term;
    summary.teacher_noise_term(ih,is)=r.teacher.noise_term;
    summary.teacher_reach_min(ih,is)=r.teacher.reach_min;
    summary.validation_nrmse(ih,is)=r.teacher.validation_nrmse;
    summary.qp_success_rate(ih,is)=r.metrics.qp_success_rate;
    summary.infeasible_count(ih,is)=r.metrics.infeasible_count;
    summary.constraint_active_rate(ih,is)=r.metrics.constraint_active_rate;
    summary.cost_mean(ih,is)=r.metrics.cost_mean;
    summary.cost_sum(ih,is)=r.metrics.cost_sum;
    summary.MAE(:,ih,is)=r.metrics.MAE;
    summary.RMSE(:,ih,is)=r.metrics.RMSE;
    summary.Bias(:,ih,is)=r.metrics.Bias;
    summary.upper_violation_rate(:,ih,is)=r.metrics.upper_violation_rate;
    summary.abs_violation_rate(:,ih,is)=r.metrics.abs_violation_rate;
    rec = struct('task_horizon',H,'seed_id',seed_id,'seed',r.seed, ...
        'config',r.config,'teacher',r.teacher,'metrics',r.metrics, ...
        'runtime_seconds',r.runtime_seconds,'result_file',f);
    if isempty(records), records=rec; else, records(end+1)=rec; end %#ok<AGROW>
end

if require_complete
    assert(numel(files)==80,'Expected 80 per-run MAT files, found %d.',numel(files));
    assert(all(summary.present,'all'),'The 4-by-20 H/seed grid is incomplete.');
end
aggregate = struct();
aggregate.schema_version = 'copyAP_crte_multistep_task_20x4_aggregate_v1';
aggregate.task_horizons = Hset; aggregate.seed_ids = seedset;
aggregate.expected_run_count = 80; aggregate.actual_run_count = numel(files);
aggregate.complete = all(summary.present,'all');
aggregate.algorithm_contract = struct('uses_true_Sigma_n',false, ...
    'task_definition','weighted future tracked-output stack t+1:t+H after exact OLS-FWL');
aggregate.summary = summary; aggregate.records = records;
aggregate.source_results_dir = results_dir;
aggregate.created_at = char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss'));
outdir = fileparts(output_file); if ~exist(outdir,'dir'), mkdir(outdir); end
save(output_file,'aggregate','-v7.3');
fprintf('copyAP aggregate saved: %s (runs=%d complete=%d)\n', ...
    output_file,aggregate.actual_run_count,aggregate.complete);
end
