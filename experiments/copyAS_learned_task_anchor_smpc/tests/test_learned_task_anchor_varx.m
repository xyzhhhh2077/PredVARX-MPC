function test_learned_task_anchor_varx
% RED/GREEN gates for learned task anchor without tracked-output indices.
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
rng(20260728,'twister');

p = 10; m = 2; ell = 5; q = 2; T = 700;
F = diag([0.91 0.76 0.58]);
G = [0.35 -0.12; 0.16 0.28; -0.09 0.22];
C = randn(p,3);
u = randn(m,T);
x = zeros(3,T+1); y = zeros(p,T);
for k = 1:T
    y(:,k) = C*x(:,k) + 0.04*randn(p,1);
    x(:,k+1) = F*x(:,k) + G*u(:,k) + 0.03*randn(3,1);
end

% Full-output task reference: no tracked channel index is supplied.
t = linspace(0,8*pi,T);
Rtask = C*[0.8*sin(t); 0.6*cos(0.7*t); 0.45*sin(0.35*t+0.4)];
opt = struct('task_reference',Rtask,'Ru',0.2*eye(m), ...
    'anchor_weights',[1 0.6 0.5 0.4],'mu_grid',[0 0.5 1]);

[A,B,P,R,S,st] = learned_task_anchor_varx(y,u,ell,q,opt);

assert(isequal(size(A),[ell ell]));
assert(isequal(size(B),[ell m]));
assert(isequal(size(P),[p ell]) && isequal(size(R),[p ell]));
assert(isequal(size(st.E_task_anchor),[p q]));
assert(isempty(st.tracked_indices),'method must not require tracked indices');
assert(norm(R'*P-eye(ell),'fro') < 1e-7,'dual identity failed');
assert(norm(st.R_task_anchor'*st.P_task_anchor-eye(q),'fro') < 1e-7, ...
    'task anchor pair is not dual');
assert(norm(P*R'*st.P_task_anchor-st.P_task_anchor,'fro') < 1e-7, ...
    'learned task subspace is not preserved');
assert(min(eig((S+S')/2)) > -1e-8,'Sigma_eps must be PSD');
assert(st.uses_true_Sigma_n==false);
assert(st.anchor_eigengap >= 0);
assert(st.spectral_radius < 1.10);

fprintf(['PASS learned task anchor: dual=%.2e preserve=%.2e gap=%.3e ' ...
    'spectral=%.4f mu=%.2f val=%.4f\n'], ...
    st.dual_error,st.task_preservation_error,st.anchor_eigengap, ...
    st.spectral_radius,st.selected_mu,st.selected_validation_nrmse);
end
