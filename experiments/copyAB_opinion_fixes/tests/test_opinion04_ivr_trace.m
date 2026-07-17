function test_opinion04_ivr_trace
% Opinion 4: export IVR iteration trace diagnostics in stats and assert
% fields always exist (including tracked-only r=0 and free-IVR r>0).
% Stop criterion is numerical stationarity, not proven monotone convergence.

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
rng(1704,'twister');

p = 9; m = 2; tracked = [1 2]; T = 400;
F = diag([0.91 0.77 0.58]);
G = [0.40 -0.15; 0.18 0.32; -0.08 0.24];
C = randn(p,3);
C(1,:) = [1 0 0];
C(2,:) = [0 1 0];
u = randn(m,T);
x = zeros(3,T+1);
y_clean = zeros(p,T);
for k = 1:T
    y_clean(:,k) = C*x(:,k);
    x(:,k+1) = F*x(:,k) + G*u(:,k) + 0.035*randn(3,1);
end
sigma = 0.06;
Sigma_n = sigma^2*eye(p);
y = y_clean + sigma*randn(p,T);

req = {'ivr_iter','ivr_trace','ivr_subspace_delta'};

% --- Case A: free IVR runs (ell > q) ---
ell_free = 5;
[~,~,~,~,~,st_free] = split_control_free_ivr_varx(y,u,ell_free,tracked,Sigma_n);
for i = 1:numel(req)
    assert(isfield(st_free,req{i}),'free case missing stats.%s',req{i});
end
assert(isscalar(st_free.ivr_iter) && st_free.ivr_iter >= 1, ...
    'free case ivr_iter should be >= 1, got %g', st_free.ivr_iter);
assert(isnumeric(st_free.ivr_trace) && isvector(st_free.ivr_trace), ...
    'free case ivr_trace must be a numeric vector');
assert(numel(st_free.ivr_trace) == st_free.ivr_iter, ...
    'ivr_trace length must equal ivr_iter');
assert(all(isfinite(st_free.ivr_trace)),'ivr_trace has non-finite entries');
assert(isscalar(st_free.ivr_subspace_delta) && isfinite(st_free.ivr_subspace_delta), ...
    'ivr_subspace_delta must be a finite scalar');
assert(st_free.ivr_subspace_delta >= 0, 'ivr_subspace_delta must be >= 0');

% --- Case B: tracked-only (r=0) still exports the three fields ---
ell_track = numel(tracked);
[~,~,~,~,~,st0] = split_control_free_ivr_varx(y,u,ell_track,tracked,Sigma_n);
for i = 1:numel(req)
    assert(isfield(st0,req{i}),'tracked-only missing stats.%s',req{i});
end
assert(st0.ivr_iter == 0,'tracked-only ivr_iter must be 0');
assert(isempty(st0.ivr_trace),'tracked-only ivr_trace must be empty');
assert(st0.ivr_subspace_delta == 0,'tracked-only ivr_subspace_delta must be 0');

fprintf(['PASS opinion04 ivr_trace: free iter=%d len(trace)=%d delta=%.3e; ' ...
    'tracked-only iter=%d empty_trace=%d delta=%g\n'], ...
    st_free.ivr_iter, numel(st_free.ivr_trace), st_free.ivr_subspace_delta, ...
    st0.ivr_iter, isempty(st0.ivr_trace), st0.ivr_subspace_delta);
end
