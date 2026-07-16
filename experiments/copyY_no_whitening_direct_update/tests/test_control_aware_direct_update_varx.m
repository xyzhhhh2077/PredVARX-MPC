function test_control_aware_direct_update_varx
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
[A,B,P,R,Sigma,stats]=control_aware_direct_update_varx(y,u,ell,tracked);
E=zeros(p,2); E(tracked,:)=eye(2);
assert(norm(P'*P-eye(ell),'fro')<1e-10,'P is not orthonormal');
assert(norm(R'*P-eye(ell),'fro')<1e-8,'dual identity failed');
assert(norm(P-R,'fro')>1e-4,'oblique extractor collapsed to R=P');
assert(norm(P*R'*E-E,'fro')<1e-8,'tracked oblique coverage failed');
assert(norm(P*R'-(P*R')','fro')>1e-4,'projector is not oblique');
assert(isequal(size(A),[ell ell]) && isequal(size(B),[ell m]));
assert(min(eig((Sigma+Sigma')/2))>-1e-8,'Sigma is not PSD');
assert(~stats.normalization_applied,'normalization must remain disabled');
assert(strcmp(stats.update_rule,'direct_eigenspace_replacement'),'wrong update rule');
fprintf('PASS no-whitening direct-update test: dual=%.3e cover=%.3e ||P-R||=%.3e asym=%.3e\n', ...
    stats.dual_error,stats.tracked_oblique_error,norm(P-R,'fro'),stats.pr_asymmetry);
end
