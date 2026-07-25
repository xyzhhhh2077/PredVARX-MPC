function status = copyAP_run_group(group_id, num_groups, overwrite)
%COPYAP_RUN_GROUP Run one deterministic shard of the 4-by-20 copyAP grid.
if nargin < 2 || isempty(num_groups), num_groups = 10; end
if nargin < 3 || isempty(overwrite), overwrite = false; end
assert(group_id >= 1 && group_id <= num_groups && group_id == round(group_id));
here = fileparts(mfilename('fullpath')); addpath(here);
results_dir = fullfile(here,'results','runs');
if ~exist(results_dir,'dir'), mkdir(results_dir); end
Hset = [1 3 6 18]; seedset = 1:20;
[HH,SS] = ndgrid(Hset,seedset); jobs = [HH(:),SS(:)];
idx = group_id:num_groups:size(jobs,1);
status = repmat(struct('ok',false,'H',NaN,'seed_id',NaN,'file','','error',''),numel(idx),1);
for ii = 1:numel(idx)
    j = idx(ii); H = jobs(j,1); seed_id = jobs(j,2);
    status(ii).H = H; status(ii).seed_id = seed_id;
    try
        r = copyAP_run_seed(H,seed_id,'ResultsDir',results_dir,'Overwrite',overwrite,'MuGrid',1,'PredictionHorizon',18);
        status(ii).ok = true; status(ii).file = r.result_file;
    catch ME
        status(ii).error = getReport(ME,'extended','hyperlinks','off');
        fprintf(2,'GROUP %d FAILED H=%d seed=%d\n%s\n',group_id,H,seed_id,status(ii).error);
    end
end
save(fullfile(results_dir,sprintf('copyAP_group_%02d_status.mat',group_id)),'status','group_id','num_groups','-v7.3');
assert(all([status.ok]),'copyAP_run_group:Failures','Group %d has %d failed jobs.',group_id,sum(~[status.ok]));
fprintf('copyAP group %02d complete: %d/%d jobs OK\n',group_id,sum([status.ok]),numel(status));
end
