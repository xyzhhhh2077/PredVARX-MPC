function test_predvarx_identify_oblique
% Regression target: the oblique realization must retain R'P = I,
% rather than replacing the IVR dual basis by R=P after QR.

rng(7, 'twister');
p = 6; ell = 2; m = 2; T = 120;
A = diag([0.82, 0.61, 0.37]);
B = randn(3,m);
C = randn(p,3); [C,~] = qr(C,0);
u = randn(m,T);
x = zeros(3,T+1); y = zeros(p,T);
for k = 1:T
    y(:,k) = C*x(:,k) + 0.05*randn(p,1);
    x(:,k+1) = A*x(:,k) + B*u(:,k) + 0.03*randn(3,1);
end

[~,~,P,Pbar,R,Rbar,~,Sigma_eps,Sigma_ebar] = ...
    predvarx_identify_oblique(y,u,ell,0.5,2,A,B,C,3,m,p);

assert(isequal(size(P), [p,ell]));
assert(isequal(size(R), [p,ell]));
assert(norm(R' * P - eye(ell), 'fro') < 1e-8, ...
    'Oblique dual bases must satisfy R''P = I.');
assert(norm(R' * Pbar, 'fro') < 1e-8, ...
    'The static complement must be annihilated by R''.');
assert(norm(Rbar' * P, 'fro') < 1e-8, ...
    'The dynamic subspace must be annihilated by Rbar''.');
assert(all(eig((Sigma_eps + Sigma_eps')/2) >= -1e-10));
assert(all(eig((Sigma_ebar + Sigma_ebar')/2) >= -1e-10));

fprintf('PASS: oblique PredVARX dual-basis invariants hold.\n');
end
