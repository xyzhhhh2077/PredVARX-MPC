function test_control_aware_oblique_ivr_varx
here=fileparts(mfilename('fullpath')); addpath(fileparts(here));
rng(17,'twister'); p=8; m=2; ell=5; tracked=[1 2]; T=500;
u=randn(m,T); F=diag([.90 .75 .55]); G=randn(3,m); x=zeros(3,T+1);
C=randn(p,3); C(1,:)=[1 0 0]; C(2,:)=[0 1 0];
noise_scale=[.02 .03 .05 .10 .20 .35 .50 .80]';
y=zeros(p,T);
for k=1:T
    y(:,k)=C*x(:,k)+noise_scale.*randn(p,1);
    x(:,k+1)=F*x(:,k)+G*u(:,k)+.04*randn(3,1);
end
[A,B,P,R,Sigma,stats]=control_aware_oblique_ivr_varx(y,u,ell,tracked);
E=zeros(p,2); E(tracked,:)=eye(2);
assert(norm(P'*P-eye(ell),'fro')<1e-10,'P is not orthonormal');
assert(norm(R'*P-eye(ell),'fro')<1e-8,'dual identity failed');
assert(norm(P-R,'fro')>1e-4,'oblique extractor collapsed to R=P');
assert(norm(P*R'*E-E,'fro')<1e-8,'tracked oblique coverage failed');
assert(norm(P*R'-(P*R')','fro')>1e-4,'projector is not oblique');
assert(isequal(size(A),[ell ell]) && isequal(size(B),[ell m]));
assert(min(eig((Sigma+Sigma')/2))>-1e-8,'Sigma is not PSD');

% Boundary regression: ell equal to the number of tracked axes means there
% is no free IVR complement. The identifier must return a valid tracked-only
% model instead of referencing an undefined loop variable.
ell_tracked=numel(tracked);
[A0,B0,P0,R0,Sigma0,stats0]=control_aware_oblique_ivr_varx(y,u,ell_tracked,tracked,0.02);
assert(stats0.ivr_iter==0,'tracked-only model must report zero IVR iterations');
assert(stats0.ivr_subspace_delta==0,'tracked-only model must report zero subspace change');
assert(isequal(size(A0),[ell_tracked ell_tracked]) && isequal(size(B0),[ell_tracked m]));
assert(norm(R0'*P0-eye(ell_tracked),'fro')<1e-8,'tracked-only dual identity failed');
assert(min(eig((Sigma0+Sigma0')/2))>-1e-8,'tracked-only Sigma is not PSD');

fprintf('PASS oblique test: dual=%.3e cover=%.3e ||P-R||=%.3e asym=%.3e; tracked-only ivr_iter=%d\n', ...
    stats.dual_error,stats.tracked_oblique_error,norm(P-R,'fro'),stats.pr_asymmetry,stats0.ivr_iter);
end
