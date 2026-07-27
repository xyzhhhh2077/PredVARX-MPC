function test_learn_output_directions
% Regression gates for learning output directions from one fixed dataset.
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
rng(20260728,'twister');

p = 8; m = 2; T = 900; H = 12;
A = diag([0.88 0.72 0.50]);
B = [0.55 0.00; 0.05 0.12; 0.00 0.06];
C = randn(p,3);
C(1,:) = [1 0 0];
C(2,:) = [0 1 0];
u = randn(m,T);
x = zeros(3,T+1); y = zeros(p,T);
for k = 1:T
    y(:,k) = C*x(:,k) + 0.02*randn(p,1);
    x(:,k+1) = A*x(:,k) + B*u(:,k) + 0.015*randn(3,1);
end

Ec = zeros(p,2); Ec(1,1) = 1; Ec(2,2) = 1;
[E_sup,st_sup] = learn_output_directions(y,u,2, ...
    struct('mode','supervised','task_outputs',Ec,'reach_horizon',H,'Ru',eye(m)));
[E_auth,st_auth] = learn_output_directions(y,u,2, ...
    struct('mode','authority','reach_horizon',H,'Ru',eye(m)));

ang_sup = max(subspace_angles(E_sup,Ec));
assert(ang_sup < 3.0,'Supervised learner did not recover selected output span.');
assert(norm(st_sup.R_anchor'*st_sup.P_anchor-eye(2),'fro') < 1e-7, ...
    'Supervised anchor pair is not dual.');
assert(norm(st_auth.R_anchor'*st_auth.P_anchor-eye(2),'fro') < 1e-7, ...
    'Authority anchor pair is not dual.');
assert(st_auth.authority_fraction > 0.55, ...
    'Learned authority subspace captures too little input authority.');
assert(all(isfinite(E_auth),'all') && rank(E_auth)==2, ...
    'Authority directions must be finite and independent.');

[Ah,Bh,Ph,Rh,Sh,st_fit] = fit_anchored_varx(y,u,E_auth,5,struct('ridge',1e-8));
assert(isequal(size(Ah),[5 5]) && isequal(size(Bh),[5 m]));
assert(norm(Rh'*Ph-eye(5),'fro') < 1e-7,'Full anchored VARX basis is not dual.');
assert(norm(Ph*Rh'*E_auth-E_auth,'fro') < 1e-7, ...
    'Full latent subspace does not preserve the learned anchor.');
assert(min(eig((Sh+Sh')/2)) > -1e-8,'Innovation covariance must be PSD.');
assert(st_fit.spectral_radius < 1.10,'Synthetic anchored VARX is unstable.');

fprintf('PASS output directions: supervised max angle=%.3f deg authority fraction=%.3f\n', ...
    ang_sup,st_auth.authority_fraction);
end

function a = subspace_angles(E,F)
[Qe,~] = qr(E,0); [Qf,~] = qr(F,0);
s = svd(Qe'*Qf); s = min(max(s,-1),1);
a = acosd(s);
end
