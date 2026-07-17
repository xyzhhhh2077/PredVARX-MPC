function test_opinion01_geometry
% Opinion 1: force geometry residuals (left/right/column/dual) via
% split_control_free_ivr_varx asserts, and explicitly verify left=0, col=0.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(here);
rng(20260718,'twister');

p = 9; m = 2; ell = 5; tracked = [1 2]; T = 650;
q = numel(tracked);
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

% Heteroscedastic + correlated sensor noise (non-orthogonal free dual).
scales = [0.025 0.030 0.045 0.080 0.130 0.210 0.330 0.480 0.700]';
Corr = eye(p);
Corr(3,4) = 0.35; Corr(4,3) = 0.35;
Corr(5,6) = -0.25; Corr(6,5) = -0.25;
D = diag(scales);
Sigma_n = D*Corr*D;
L = chol(Sigma_n,'lower');
y = y_clean + L*randn(p,T);

[Ahat,Bhat,P,R,Sigma_eps,stats] = split_control_free_ivr_varx( ...
    y,u,ell,tracked,Sigma_n);

% Internal assert path is on by default.
assert(isfield(stats,'enforce_geometry') && isequal(stats.enforce_geometry,true), ...
    'enforce_geometry should default to true');

% Explicit external checks requested by Opinion 1: left=0, column=0
% (and the full geometry quartet for completeness).
tol = 1e-8;
assert(stats.tracked_left_error <= tol, ...
    'left residual not ~0: %.3e', stats.tracked_left_error);
assert(stats.tracked_column_error <= tol, ...
    'column residual not ~0: %.3e', stats.tracked_column_error);
assert(stats.tracked_right_error <= tol, ...
    'right residual not ~0: %.3e', stats.tracked_right_error);
assert(stats.dual_error <= tol, ...
    'dual residual not ~0: %.3e', stats.dual_error);

% Sanity: tracked columns of R/P are exact basis vectors E.
E = zeros(p,q); E(tracked,:) = eye(q);
assert(norm(R(:,1:q)-E,'fro') <= tol, 'R tracked columns != E');
assert(norm(P(:,1:q)-E,'fro') <= tol, 'P tracked columns != E');
assert(isequal(size(Ahat),[ell ell]) && isequal(size(Bhat),[ell m]));
assert(min(eig((Sigma_eps+Sigma_eps')/2)) > -1e-8);

% Optional switch can be turned off without throwing (geometry still reported).
[~,~,~,~,~,stats_off] = split_control_free_ivr_varx( ...
    y,u,ell,tracked,Sigma_n,'enforce_geometry',false);
assert(isequal(stats_off.enforce_geometry,false), ...
    'enforce_geometry=false not honored');
assert(stats_off.tracked_left_error <= tol && stats_off.tracked_column_error <= tol, ...
    'residuals should still be numerically zero even when assert is off');

fprintf(['PASS opinion01 geometry: left=%.3e right=%.3e col=%.3e dual=%.3e ' ...
    '(enforce=%d, free_oblique=%.3e)\n'], ...
    stats.tracked_left_error, stats.tracked_right_error, ...
    stats.tracked_column_error, stats.dual_error, ...
    stats.enforce_geometry, stats.free_oblique_norm);
end
