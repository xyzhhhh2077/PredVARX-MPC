function cfg = build_cdr_closed_loop_config()
%BUILD_CDR_CLOSED_LOOP_CONFIG AU-compatible control setup for copyAX CDR.
here = fileparts(mfilename('fullpath'));
identified = load(fullfile(here,'results','copyAX_cdr_soft_preference_data.mat'));
plant = load(fullfile(here,'data','copyAX_cdr_plant.mat'));

cfg = struct();
cfg.here = here;
cfg.p = 30; cfg.m = 8; cfg.q = 2; cfg.ell = 5;
cfg.N = 18; cfg.T = 1200; cfg.warmup = 150;
cfg.alpha_joint = 0.10;
cfg.u_min = -1; cfg.u_max = 1;
cfg.sensor_noise_std = 0.01;
cfg.seed = 20260802;
cfg.Etask = identified.Etask;
cfg.model = identified.model;
cfg.y_scale = identified.y_scale;
cfg.y_offset = identified.y_offset;
cfg.u_scale = identified.u_scale;
cfg.u_offset = identified.u_offset;
cfg.A_plant = plant.A_plant;
cfg.B_plant = plant.B_plant;
cfg.C_plant = plant.C_plant;
cfg.x0 = plant.x0;

% Preserve AU exactly in the standardized coordinates used for identification.
% Only the plant bridge converts the QP action back to physical CDR units.
model = cfg.model;
model.Sigma_obs = diag((cfg.sensor_noise_std./cfg.y_scale).^2);
cfg.model_control = model;
cfg.Q = 80*(cfg.Etask*cfg.Etask');
cfg.Ru = 0.18*eye(cfg.m);
cfg.H = cfg.Etask';

% Choose an attainable, non-trivial two-axis reference using only the
% identified AU model.  The true CDR matrices remain simulation-only and do
% not participate in controller or reference design.
Gss=cfg.H*model.P*((eye(cfg.ell)-model.A)\model.B);
[Ug,Sg,~] = svd(Gss,'econ');
gain = max(diag(Sg));
if isempty(gain) || gain < 1e-8, error('CDR task steady-state gain is singular.'); end
cfg.reference_amplitude = min(0.20,0.35*gain);
cfg.reference_basis = Ug(:,1:cfg.q);

% A chance bound must contain both the requested task level and the finite-
% horizon Gaussian tightening.  Derive it from the identified AU covariance;
% choosing only a multiple of the reference makes the QP infeasible before
% any control action is considered.
Sigma_z=zeros(cfg.ell); risk_each=cfg.alpha_joint/(2*cfg.q*cfg.N);
z_quantile=norminv(1-risk_each); max_tight=0;
for j=1:cfg.N
    Sigma_z=model.A*Sigma_z*model.A'+model.Sigma_eps;
    Sigma_task=cfg.H*model.P*Sigma_z*model.P'*cfg.H';
    max_tight=max(max_tight,max(z_quantile*sqrt(max(diag(Sigma_task),0))));
end
cfg.max_chance_tightening=max_tight;
cfg.task_limit=1.25*(cfg.reference_amplitude+max_tight);
end
