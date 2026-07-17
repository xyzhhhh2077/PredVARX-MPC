function test_opinion08_fallback_flag
% TEST_OPINION08_FALLBACK_FLAG  Logic unit test (no real QP).
% Opinion 8: sat(u_mean) fallback is always uncertified for chance constraints.
% A one-step deterministic mean check is only a diagnostic hard-bound residual.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(fullfile(root, 'lib'));

%% Synthetic model / opt (tiny, no quadprog)
p = 4; ell = 2; nu = 2; nq = 2;
model.A = 0.8 * eye(ell);
model.B = [0.3 0; 0 0.25];
model.P = [eye(ell); zeros(p-ell, ell)];
model.R = model.P;
model.y_mean = zeros(p,1);
model.u_mean = zeros(nu,1);
model.Sigma_eps = 1e-3*eye(ell);
model.Sigma_obs = 1e-3*eye(p);

opt.H = zeros(nq, p);
opt.H(1,1) = 1;
opt.H(2,2) = 1;
opt.h = [1.0; 1.0];   % y_max = 1 on tracked axes
opt.u_min = -2; opt.u_max = 2;

%% Case A: mean prediction safely under y_max — STILL uncertified
y_safe = [0.1; 0.1; 0; 0];
uk = min(max(model.u_mean, opt.u_min), opt.u_max);  % sat(u_mean)
infoA = fallback_certify_step(y_safe, model, opt, uk);

assert(strcmp(infoA.cc_cert_level, 'uncertified_fallback'), ...
    'Case A: cert level must be uncertified_fallback');
assert(infoA.exitflag == -1, 'Case A: exitflag must be -1');
assert(isnan(infoA.max_cc_violation), ...
    'Case A: max_cc_violation must be NaN (no QP residual)');
assert(infoA.uncertified_fallback == true, 'Case A: uncertified flag');
assert(~isempty(infoA.mean_pred), 'Case A: mean prediction should exist');
assert(infoA.det_mean_safe == true, 'Case A: mean should be under y_max');
assert(infoA.det_mean_violation <= 0, 'Case A: det residual <= 0');
assert(contains(lower(infoA.note), 'not chance-constraint'), ...
    'Case A: note must deny chance-constraint certification');

%% Case B: mean prediction violates hard y_max — still uncertified, det unsafe
% Drive z large via y so one-step mean exceeds h.
y_hot = [3.0; 3.0; 0; 0];
infoB = fallback_certify_step(y_hot, model, opt, uk);

assert(strcmp(infoB.cc_cert_level, 'uncertified_fallback'), ...
    'Case B: cert level must remain uncertified_fallback');
assert(infoB.exitflag == -1, 'Case B: exitflag must be -1');
assert(isnan(infoB.max_cc_violation), 'Case B: max_cc_violation NaN');
assert(infoB.det_mean_safe == false, 'Case B: mean should violate y_max');
assert(infoB.det_mean_violation > 0, 'Case B: det residual > 0');
assert(contains(lower(infoB.note), 'not chance-constraint'), ...
    'Case B: note must deny chance-constraint certification');

%% Case C: bookkeeping mimics runner split (qp vs uncertified_fallback)
cc_cert_level = {'qp','qp','uncertified_fallback','qp','uncertified_fallback'};
exitflag = [1, 1, -1, 1, -1];
max_cc_violation = [-0.1, -0.05, NaN, 0.0, NaN];

qp_certified_count = sum(strcmp(cc_cert_level, 'qp'));
uncertified_fallback_count = sum(strcmp(cc_cert_level, 'uncertified_fallback'));
assert(qp_certified_count == 3, 'Case C: qp count');
assert(uncertified_fallback_count == 2, 'Case C: fallback count');
assert(mean(exitflag > 0) == 0.6, 'Case C: qp success rate');

qp_mask = strcmp(cc_cert_level, 'qp');
max_qp = max(max_cc_violation(qp_mask & ~isnan(max_cc_violation)));
assert(abs(max_qp - 0.0) < 1e-15, 'Case C: only QP residuals enter max');

% Fallbacks must not be counted as chance-constraint-guaranteed successes.
assert(qp_certified_count + uncertified_fallback_count == numel(cc_cert_level));
assert(uncertified_fallback_count > 0, 'Case C: has fallbacks');
% Explicit policy: do NOT claim fallback inherits chance guarantees.
assert(~any(strcmp(cc_cert_level, 'qp') & exitflag <= 0), ...
    'Case C: exitflag<=0 must never be labeled qp');

fprintf(['PASS opinion08: cert levels ok; safe-mean still uncertified; ' ...
    'violating-mean uncertified; bookkeeping splits qp=%d fallback=%d\n'], ...
    qp_certified_count, uncertified_fallback_count);
end
