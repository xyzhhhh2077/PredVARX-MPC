function test_opinion10_input_residualize
% Opinion 10: optional input residualization IVR mode (default OFF).
% 1) Default path keeps original geometry assertions and flag false.
% 2) residualize-on-U mode runs without crash; geometry still holds;
%    stats.ivr_input_conditional is true (candidate only, not proved optimum).

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
rng(1710,'twister');

p = 9; m = 2; ell = 5; tracked = [1 2]; T = 650;
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

scales = [0.025 0.030 0.045 0.080 0.130 0.210 0.330 0.480 0.700]';
Corr = eye(p);
Corr(3,4)=0.35; Corr(4,3)=0.35;
Corr(5,6)=-0.25; Corr(6,5)=-0.25;
D = diag(scales);
Sigma_het = D*Corr*D;
L = chol(Sigma_het,'lower');
y = y_clean + L*randn(p,T);
E = zeros(p,numel(tracked)); E(tracked,:) = eye(numel(tracked));

% ---- Case A: default (input_residualize off / omitted) ----
[Ad,Bd,Pd,Rd,Sd,std] = split_control_free_ivr_varx(y,u,ell,tracked,Sigma_het);
assert(isfield(std,'ivr_input_conditional'),'missing stats.ivr_input_conditional');
assert(std.ivr_input_conditional == false, ...
    'default path must set ivr_input_conditional=false');
assert(norm(Rd'*Pd-eye(ell),'fro') < 1e-8,'default dual identity failed');
assert(norm(Pd*Rd'*E-E,'fro') < 1e-8,'default right coverage failed');
assert(norm(E'*Pd*Rd'-E','fro') < 1e-8,'default left value preservation failed');
assert(norm(Rd(:,1:numel(tracked))-E,'fro') < 1e-12,'default tracked columns changed');
assert(std.free_noise_objective <= std.free_noise_baseline + 1e-10, ...
    'default free-noise objective regressed');
assert(isequal(size(Ad),[ell ell]) && isequal(size(Bd),[ell m]));
assert(min(eig((Sd+Sd')/2)) > -1e-8,'default Sigma not PSD');

% Explicit false name-value must match omitted default
[~,~,P0,R0,~,st0] = split_control_free_ivr_varx( ...
    y,u,ell,tracked,Sigma_het,'input_residualize',false);
assert(st0.ivr_input_conditional == false);
assert(norm(P0-Pd,'fro') < 1e-12 && norm(R0-Rd,'fro') < 1e-12, ...
    'explicit false must match default path');

% ---- Case B: residualize-on-U candidate mode ----
[Ar,Br,Pr,Rr,Sr,str] = split_control_free_ivr_varx( ...
    y,u,ell,tracked,Sigma_het,'input_residualize',true);
assert(str.ivr_input_conditional == true, ...
    'residualize mode must set ivr_input_conditional=true');
assert(norm(Rr'*Pr-eye(ell),'fro') < 1e-8,'residualize dual identity failed');
assert(norm(Pr*Rr'*E-E,'fro') < 1e-8,'residualize right coverage failed');
assert(norm(E'*Pr*Rr'-E','fro') < 1e-8,'residualize left value preservation failed');
assert(norm(Rr(:,1:numel(tracked))-E,'fro') < 1e-12, ...
    'residualize tracked extractor columns changed');
assert(str.free_noise_objective <= str.free_noise_baseline + 1e-10, ...
    'residualize free-noise objective regressed vs R=P free block');
assert(isequal(size(Ar),[ell ell]) && isequal(size(Br),[ell m]));
assert(min(eig((Sr+Sr')/2)) > -1e-8,'residualize Sigma not PSD');
assert(isfinite(str.ivr_iter) && str.ivr_iter >= 1, ...
    'residualize IVR did not iterate');
assert(all(isfinite(str.ivr_trace(:))),'residualize ivr_trace not finite');

% Free subspace may differ from default (candidate conditionalization);
% do not assert equality — only that both paths are well-posed.
sub_delta = norm(Pr*Pr'-Pd*Pd','fro');

fprintf(['PASS opinion10 input residualize: default flag=%d; ' ...
    'residualize flag=%d; free-subspace delta=%.3e; ' ...
    'default dual=%.3e residualize dual=%.3e\\n'], ...
    std.ivr_input_conditional, str.ivr_input_conditional, sub_delta, ...
    norm(Rd'*Pd-eye(ell),'fro'), norm(Rr'*Pr-eye(ell),'fro'));
end
