function aggregate = copyAP_run_all_10way(varargin)
%COPYAP_RUN_ALL_10WAY Run the full 4 horizons x 20 seeds using 10 workers.
% Default invocation executes exactly 80 full-scale jobs and then writes the
% validated -v7.3 aggregate MAT. Existing valid per-run files are reused
% unless Overwrite=true.

here = fileparts(mfilename('fullpath')); addpath(here);
p = inputParser;
p.addParameter('Overwrite',false,@(x)islogical(x)&&isscalar(x));
p.addParameter('ResultsDir',fullfile(here,'results','runs'),@(x)ischar(x)||isstring(x));
p.addParameter('AggregateFile',fullfile(here,'results', ...
    'copyAP_crte_multistep_task_20x4_aggregate.mat'),@(x)ischar(x)||isstring(x));
p.parse(varargin{:}); args=p.Results;
results_dir=char(args.ResultsDir); if ~exist(results_dir,'dir'), mkdir(results_dir); end
% Select MathWorks quadprog once before starting the shared-memory pool.
mathworks_qp_dir=fullfile(matlabroot,'toolbox','optim','optim');
qp_before=which('quadprog'); qp_shadow_dir='';
if ~startsWith(qp_before,mathworks_qp_dir,'IgnoreCase',true)
    qp_shadow_dir=fileparts(qp_before); rmpath(qp_shadow_dir); rehash;
    qp_path_guard=onCleanup(@()addpath(qp_shadow_dir,'-begin')); %#ok<NASGU>
end
assert(startsWith(which('quadprog'),mathworks_qp_dir,'IgnoreCase',true));
Hset=[1 3 6 18]; seedset=1:20;
[HH,SS]=ndgrid(Hset,seedset); jobs=[HH(:),SS(:)];
assert(size(jobs,1)==80 && numel(unique(string(jobs(:,1))+"_"+string(jobs(:,2))))==80);

pool=gcp('nocreate');
if ~isempty(pool) && (~isa(pool,'parallel.ThreadPool') || pool.NumWorkers~=10)
    delete(pool); pool=[];
end
if isempty(pool), pool=parpool('Threads',10); end
assert(pool.NumWorkers==10,'copyAP dispatcher requires exactly 10 workers.');

fprintf('copyAP dispatch: %d jobs, task_horizons=[1 3 6 18], seeds=1:20, workers=%d\n', ...
    size(jobs,1),pool.NumWorkers);
t0=tic; status=cell(size(jobs,1),1);
parfor j=1:size(jobs,1)
    H=jobs(j,1); seed_id=jobs(j,2);
    try
        r=copyAP_run_seed(H,seed_id,'ResultsDir',results_dir,'Overwrite',args.Overwrite);
        status{j}=struct('ok',true,'H',H,'seed_id',seed_id, ...
            'runtime_seconds',r.runtime_seconds,'file',r.result_file,'error','');
    catch ME
        status{j}=struct('ok',false,'H',H,'seed_id',seed_id, ...
            'runtime_seconds',NaN,'file','','error',getReport(ME,'extended','hyperlinks','off'));
    end
end
elapsed=toc(t0); ok=cellfun(@(s)s.ok,status);
dispatch=struct('schema_version','copyAP_dispatch_v1','jobs',jobs,'status',{status}, ...
    'num_workers',pool.NumWorkers,'pool_type',class(pool), ...
    'elapsed_seconds',elapsed,'num_ok',sum(ok), ...
    'num_failed',sum(~ok),'created_at',char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss')));
save(fullfile(results_dir,'copyAP_dispatch_status.mat'),'dispatch','-v7.3');
fprintf('copyAP dispatch finished: ok=%d failed=%d elapsed=%.1fs\n',sum(ok),sum(~ok),elapsed);
if any(~ok)
    bad=find(~ok);
    for k=1:numel(bad)
        s=status{bad(k)}; fprintf(2,'FAILED H=%d seed=%d\n%s\n',s.H,s.seed_id,s.error);
    end
    aggregate_copyAP_results(results_dir,char(args.AggregateFile),false);
    error('copyAP_run_all_10way:WorkerFailures','%d of 80 jobs failed.',sum(~ok));
end
aggregate=aggregate_copyAP_results(results_dir,char(args.AggregateFile),true);
end
