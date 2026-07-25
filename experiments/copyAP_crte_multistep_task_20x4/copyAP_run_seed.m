function run = copyAP_run_seed(task_horizon, seed_id, varargin)
%COPYAP_RUN_SEED Run one copyAP task-horizon/seed experiment and save one MAT.
%   RUN = COPYAP_RUN_SEED(H,SEED_ID) runs one of H=[1 3 6 18] and one of
%   SEED_ID=1:20. Plant/offline and closed-loop RNG streams are recorded
%   separately. True plant Sigma_n is never passed to identification or MPC.

here = fileparts(mfilename('fullpath'));
% This MATLAB installation has a stale MOSEK quadprog wrapper ahead of the
% MathWorks solver. Remove only that shadow directory for this standalone
% call; a thread dispatcher preselects the MathWorks solver before parfor.
quadprog_before = which('quadprog');
mathworks_qp_dir = fullfile(matlabroot,'toolbox','optim','optim');
qp_shadow_dir = '';
if ~startsWith(quadprog_before,mathworks_qp_dir,'IgnoreCase',true)
    qp_shadow_dir = fileparts(quadprog_before);
    rmpath(qp_shadow_dir); rehash;
    path_guard = onCleanup(@() addpath(qp_shadow_dir,'-begin')); %#ok<NASGU>
end
quadprog_used = which('quadprog');
assert(startsWith(quadprog_used,mathworks_qp_dir,'IgnoreCase',true), ...
    'copyAP requires the MathWorks Optimization Toolbox quadprog.');
p = inputParser;
p.addRequired('task_horizon', @(x) isscalar(x) && any(x == [1 3 6 18]));
p.addRequired('seed_id', @(x) isscalar(x) && x == round(x) && x >= 1 && x <= 20);
p.addParameter('Smoke', false, @(x) islogical(x) && isscalar(x));
p.addParameter('Overwrite', false, @(x) islogical(x) && isscalar(x));
p.addParameter('ResultsDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('TOff', [], @(x) isempty(x) || (isscalar(x) && x >= 300));
p.addParameter('TCl', [], @(x) isempty(x) || (isscalar(x) && x >= 50));
p.addParameter('NumRandomSubspaces', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
p.addParameter('MuGrid', [], @(x) isempty(x) || isvector(x));
p.addParameter('PredictionHorizon', [], @(x) isempty(x) || (isscalar(x) && x >= 1));
p.parse(task_horizon, seed_id, varargin{:});
args = p.Results;

if args.Smoke
    if isempty(args.TOff), T_off = 500; else, T_off = args.TOff; end
    if isempty(args.TCl), T_cl = 120; else, T_cl = args.TCl; end
    if isempty(args.NumRandomSubspaces), nrandom = 2; else, nrandom = args.NumRandomSubspaces; end
    if isempty(args.MuGrid), mu_grid = 1; else, mu_grid = args.MuGrid; end
    if isempty(args.PredictionHorizon), N = 6; else, N = args.PredictionHorizon; end
    default_results = fullfile(here, 'results', 'smoke');
    prefix = 'smoke';
else
    if isempty(args.TOff), T_off = 1500; else, T_off = args.TOff; end
    if isempty(args.TCl), T_cl = 1200; else, T_cl = args.TCl; end
    if isempty(args.NumRandomSubspaces), nrandom = 30; else, nrandom = args.NumRandomSubspaces; end
    if isempty(args.MuGrid), mu_grid = 1; else, mu_grid = args.MuGrid; end
    if isempty(args.PredictionHorizon), N = 18; else, N = args.PredictionHorizon; end
    default_results = fullfile(here, 'results', 'runs');
    prefix = 'run';
end
assert(T_off > max(N, task_horizon) + 100, 'Offline record is too short for requested horizons.');
results_dir = char(args.ResultsDir);
if isempty(results_dir), results_dir = default_results; end
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
result_file = fullfile(results_dir, sprintf('copyAP_%s_H%02d_seed%02d.mat', prefix, task_horizon, seed_id));
if exist(result_file, 'file') && ~args.Overwrite
    s = load(result_file, 'run');
    run = s.run;
    fprintf('copyAP reuse H=%d seed=%d: %s\n', task_horizon, seed_id, result_file);
    return;
end

tic_run = tic;
seed = struct();
seed.seed_id = seed_id;
seed.plant_offline = 2026071000 + seed_id;
seed.teacher_search = 2026072500 + seed_id;
seed.closed_loop = 2026073000 + seed_id;

%% Process-like plant and controller parameters
rng(seed.plant_offline, 'twister');
n = 6; m = 3; ny = 30; ell = 5; tracked = [1 2];
sw = 0.045; se = 0.055;
sensor_noise_mode = 'heteroscedastic_correlated_plant_only';
noise_cycle = 400;
sw_min = 0.020; sw_max = 0.090;
se_min = 0.025; se_max = 0.100;
noise_phase_e = pi/3;
u_min = -3.0; u_max = 3.0; y_max = 2.00; alpha_joint = 0.10;

A = diag([0.94 0.88 0.78 0.64 0.50 0.35]);
A(1,2) = 0.10; A(2,3) = -0.06; A(3,4) = 0.05; A(4,5) = 0.04;
B = [0.34 -0.10  0.05; 0.12  0.28 -0.06; 0.05  0.12  0.24; ...
    -0.05  0.06  0.18; 0.02 -0.10  0.14; 0.08  0.02 -0.08];
C = zeros(ny,n);
C(1,1) = 1.00; C(1,3) = 0.16;
C(2,2) = 1.00; C(2,4) = -0.12;
for i = 3:ny
    C(i,:) = 0.45*randn(1,n);
    C(i,:) = C(i,:) / max(norm(C(i,:)), 1e-12);
end

%% Plant-only sensor-noise truth (for sampling and post-hoc audit only)
sensor_rel = linspace(0.55,1.65,ny)';
Corr_n = eye(ny);
Corr_n(3,4) = 0.30; Corr_n(4,3) = 0.30;
Corr_n(5,6) = -0.22; Corr_n(6,5) = -0.22;
Corr_n(8,9) = 0.18; Corr_n(9,8) = 0.18;
D_n = diag(se*sensor_rel);
Sigma_n_plant = D_n*Corr_n*D_n;
L_n_plant = chol(Sigma_n_plant, 'lower');

%% Offline excitation and data
u_off = 1.20*randn(m,T_off);
x_off = zeros(n,T_off+1);
y_off = zeros(ny,T_off);
for k = 1:T_off
    y_off(:,k) = C*x_off(:,k) + L_n_plant*randn(ny,1);
    x_off(:,k+1) = A*x_off(:,k) + B*u_off(:,k) + sw*randn(n,1);
end

%% CRTE unknown-noise profiled teacher with an actual future-task stack
teacher_opt = struct();
teacher_opt.mu_grid = mu_grid;
teacher_opt.alpha = 1.0;
teacher_opt.beta = 1.0;
teacher_opt.prediction_horizon = N;
teacher_opt.omega = ones(1,N)/N;
teacher_opt.task_horizon = task_horizon;
teacher_opt.task_omega = ones(1,task_horizon)/task_horizon;
teacher_opt.Ru = 0.18*eye(m);
teacher_opt.val_fraction = 0.25;
teacher_opt.ridge = 1e-8;
teacher_opt.rank_tol = 1e-9;
teacher_opt.reach_tau = 1e-10;
teacher_opt.num_random_subspaces = nrandom;
teacher_opt.seed = seed.teacher_search;
[Ahat,Bhat,Phat,Rhat,Sigma_eps,stats] = ...
    crte_profiled_teacher_unknown_noise_multitask(y_off,u_off,ell,tracked,teacher_opt);
assert(stats.task_horizon == task_horizon && stats.task_future_rows == numel(tracked)*task_horizon);
assert(stats.selected_mu == 1, 'copyAP fair ablation requires mu=1.');
assert(~stats.uses_true_Sigma_n, 'Unknown-noise contract was violated.');

model = struct('A',Ahat,'B',Bhat,'P',Phat,'R',Rhat, ...
    'y_mean',stats.y_mean,'u_mean',stats.u_mean, ...
    'Sigma_eps',(Sigma_eps+Sigma_eps')/2);
IminusPR = eye(ny)-Phat*Rhat';
obs_off = IminusPR*(y_off-stats.y_mean);
obs_off = obs_off-mean(obs_off,2);
Sigma_obs_proxy = (obs_off*obs_off')/max(size(obs_off,2)-1,1);
Sigma_obs_proxy = (Sigma_obs_proxy+Sigma_obs_proxy')/2;
ridge_obs = 1e-8*max(trace(Sigma_obs_proxy)/max(ny,1),1);
model.Sigma_obs = Sigma_obs_proxy+ridge_obs*eye(ny);
assert(~isfield(model,'Sigma_n'), 'True Sigma_n must not enter the controller model.');

opt = struct();
opt.N = N; opt.Q = zeros(ny); opt.Q(1,1) = 80; opt.Q(2,2) = 80;
opt.Ru = 0.18*eye(m); opt.u_min = u_min; opt.u_max = u_max;
opt.H = zeros(numel(tracked),ny); opt.H(:,tracked) = eye(numel(tracked));
opt.h = y_max*ones(numel(tracked),1); opt.alpha_joint = alpha_joint;

%% Closed-loop reference
Rf = zeros(ny,T_cl);
levels = [0.25 1.50 0.65 1.85 0.45; 0.35 1.25 1.75 0.80 1.55];
seg_len = floor(T_cl/5);
for s = 1:5
    ix = (s-1)*seg_len+1:min(s*seg_len,T_cl);
    Rf(1,ix) = levels(1,s); Rf(2,ix) = levels(2,s);
end
if seg_len*5 < T_cl
    Rf(1,seg_len*5+1:end) = levels(1,end);
    Rf(2,seg_len*5+1:end) = levels(2,end);
end

%% Closed-loop SMPC
rng(seed.closed_loop, 'twister');
[sigma_w_profile,sigma_e_profile] = smooth_noise_profile( ...
    T_cl,noise_cycle,sw_min,sw_max,se_min,se_max,noise_phase_e);
x = zeros(n,1); y = zeros(ny,T_cl); yhat = zeros(ny,T_cl); u = zeros(m,T_cl);
exitflag = zeros(1,T_cl); max_cc_violation = nan(1,T_cl); costJ = nan(1,T_cl);
estimated_sigma_eps = nan(1,T_cl); estimated_sigma_obs = nan(1,T_cl);
true_sigma_w = nan(1,T_cl); true_sigma_e = nan(1,T_cl);
noise_window = 40; eps_buffer = zeros(ell,noise_window); eps_count = 0;
obs_buffer = zeros(ny,noise_window); obs_count = 0; infeasible_count = 0;
for k = 1:T_cl
    vk = (sigma_e_profile(k)/se)*L_n_plant*randn(ny,1);
    yk = C*x + vk; true_sigma_e(k) = norm(vk)/sqrt(ny); y(:,k) = yk;
    zk_now = model.R'*(yk-model.y_mean);
    obs_res = IminusPR*(yk-model.y_mean);
    io = mod(k-1,noise_window)+1; obs_buffer(:,io) = obs_res;
    obs_count = min(obs_count+1,noise_window);
    if k >= 2
        eps_res = zk_now-model.A*z_prev-model.B*(u(:,k-1)-model.u_mean);
        ie = mod(k-2,noise_window)+1; eps_buffer(:,ie) = eps_res;
        eps_count = min(eps_count+1,noise_window);
    end
    if eps_count >= 5
        Ewin = eps_buffer(:,1:eps_count);
        model.Sigma_eps = (Ewin*Ewin')/max(eps_count,1)+1e-8*eye(ell);
        model.Sigma_eps = (model.Sigma_eps+model.Sigma_eps')/2;
    end
    if obs_count >= 5
        Owin = obs_buffer(:,1:obs_count)-mean(obs_buffer(:,1:obs_count),2);
        sigma_obs_k = norm(Owin,'fro')/sqrt(max((ny-ell)*(obs_count-1),1));
        model.Sigma_obs = max(sigma_obs_k^2,1e-8)*mean(diag(Sigma_obs_proxy))*eye(ny)/ny+1e-8*eye(ny);
    end
    estimated_sigma_eps(k) = sqrt(trace(model.Sigma_eps)/ell);
    estimated_sigma_obs(k) = sqrt(trace(model.Sigma_obs)/ny);
    rk = Rf(:,min(k+1,T_cl));
    try
        [~,ypred,U,info] = centered_smpc_step(yk,rk,model,opt);
        uk = U(1:m); yhat(:,k) = ypred; exitflag(k) = info.exitflag;
        max_cc_violation(k) = max(info.A_ch*U-info.b_ch); costJ(k) = info.cost;
    catch ME
        infeasible_count = infeasible_count+1; exitflag(k) = -1;
        uk = min(max(model.u_mean,u_min),u_max);
        if k>1, yhat(:,k)=yhat(:,k-1); else, yhat(:,k)=model.y_mean; end
        if infeasible_count<=3
            fprintf('copyAP H=%d seed=%d QP fallback k=%d: %s\n',task_horizon,seed_id,k,ME.message);
        end
    end
    u(:,k) = uk;
    wk = sigma_w_profile(k)*randn(n,1); true_sigma_w(k) = norm(wk)/sqrt(n);
    x = A*x+B*uk+wk; z_prev = zk_now;
end

%% Metrics and versioned result
warm_start = min(151, max(1,floor(T_cl/4)+1)); warm = warm_start:T_cl;
err = y(tracked,warm)-Rf(tracked,warm);
metrics = struct();
metrics.MAE = mean(abs(err),2); metrics.RMSE = sqrt(mean(err.^2,2));
metrics.Bias = mean(err,2);
metrics.upper_violation_count = sum(y(tracked,:)>y_max,2);
metrics.upper_violation_rate = metrics.upper_violation_count/T_cl;
metrics.abs_violation_count = sum(abs(y(tracked,:))>y_max,2);
metrics.abs_violation_rate = metrics.abs_violation_count/T_cl;
metrics.qp_success_rate = mean(exitflag>0); metrics.infeasible_count = infeasible_count;
metrics.constraint_active_rate = mean(max_cc_violation(warm)>-1e-3,'omitnan');
metrics.cover_warm_qp_success_time_fraction = mean(exitflag(warm)>0);
metrics.cost_mean = mean(costJ(warm),'omitnan'); metrics.cost_sum = sum(costJ(warm),'omitnan');
metrics.u_rms = sqrt(mean(u.^2,2));
metrics.u_sat_rate = mean(abs(u-u_min)<1e-8 | abs(u-u_max)<1e-8,2);
finite_qp = max_cc_violation(isfinite(max_cc_violation));
if isempty(finite_qp), metrics.max_recorded_qp_constraint = NaN;
else, metrics.max_recorded_qp_constraint = max(finite_qp); end

teacher = struct();
teacher.objective = stats.selected_teacher_objective;
teacher.prediction_term = stats.selected_prediction_term;
teacher.task_term = stats.selected_task_term;
teacher.noise_term = stats.selected_noise_term;
teacher.reach_min = stats.selected_reach_min;
teacher.validation_nrmse = stats.selected_validation_nrmse;
teacher.selected_mu = stats.selected_mu; teacher.selected_source = stats.selected_source;
teacher.best_index = stats.best_index; teacher.num_candidates = stats.num_candidates;
teacher.num_feasible = stats.num_feasible; teacher.task_horizon = stats.task_horizon;
teacher.task_omega = stats.task_omega;

config = struct();
config.task_horizon = task_horizon; config.task_horizons_ablation = [1 3 6 18];
config.seed_ids_ablation = 1:20; config.prediction_horizon = N;
config.T_off = T_off; config.T_cl = T_cl; config.n = n; config.m = m;
config.p = ny; config.ell = ell; config.tracked = tracked;
config.mu_grid = mu_grid; config.num_random_subspaces = nrandom;
config.smoke = args.Smoke; config.noise_window = noise_window;
config.sensor_noise_mode = sensor_noise_mode; config.u_min = u_min;
config.u_max = u_max; config.y_max = y_max; config.alpha_joint = alpha_joint;
config.quadprog_before_local_guard = quadprog_before;
config.quadprog_used = quadprog_used;

run = struct();
run.schema_version = 'copyAP_crte_multistep_task_20x4_v1';
run.completed = true; run.config = config; run.seed = seed;
run.algorithm_contract = struct('uses_true_Sigma_n',false, ...
    'noise_object',stats.noise_object, ...
    'task_definition','weighted stacked tracked outputs t+1:t+H after exact OLS-FWL', ...
    'search_claim',stats.search_claim);
run.teacher = teacher; run.metrics = metrics; run.teacher_stats = stats;
run.plant_truth = struct('A',A,'B',B,'C',C,'Sigma_n_plant',Sigma_n_plant, ...
    'sensor_rel',sensor_rel,'Corr_n',Corr_n,'sw',sw,'se',se);
run.identification = struct('Ahat',Ahat,'Bhat',Bhat,'Phat',Phat,'Rhat',Rhat, ...
    'Sigma_eps_offline',Sigma_eps,'Sigma_obs_proxy_offline',Sigma_obs_proxy);
run.offline = struct('x',x_off,'y',y_off,'u',u_off);
run.closed_loop = struct('Rf',Rf,'y',y,'yhat',yhat,'u',u,'exitflag',exitflag, ...
    'max_cc_violation',max_cc_violation,'cost',costJ, ...
    'estimated_sigma_eps',estimated_sigma_eps,'estimated_sigma_obs',estimated_sigma_obs, ...
    'true_sigma_w',true_sigma_w,'true_sigma_e',true_sigma_e, ...
    'sigma_w_profile',sigma_w_profile,'sigma_e_profile',sigma_e_profile);
run.runtime_seconds = toc(tic_run);
run.created_at = char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss'));
run.result_file = result_file;
save(result_file, 'run', '-v7.3');

fprintf(['copyAP H=%d seed=%02d done %.2fs teacher=[J %.6g pred %.6g task %.6g noise %.6g] ' ...
    'MAE=[%.4f %.4f] RMSE=[%.4f %.4f] Bias=[%.4f %.4f] QP=%.4f upper=[%.4f %.4f] cost=%.4f\n'], ...
    task_horizon,seed_id,run.runtime_seconds,teacher.objective,teacher.prediction_term, ...
    teacher.task_term,teacher.noise_term,metrics.MAE(1),metrics.MAE(2), ...
    metrics.RMSE(1),metrics.RMSE(2),metrics.Bias(1),metrics.Bias(2), ...
    metrics.qp_success_rate,metrics.upper_violation_rate(1), ...
    metrics.upper_violation_rate(2),metrics.cost_mean);
fprintf('saved %s\n', result_file);
end
