function test_predvarx_identify_strict_whitened
here=fileparts(mfilename('fullpath')); addpath(fileparts(here));
% Regression test for the rank-full, s=1 Mo--Qin-style whitening ablation:
% whitened IVR eigenvectors are de-normalized directly, with no QR or
% post-hoc SVD realignment. This does not certify complete Algorithm 1.

rng(11, 'twister');
p=6; n=3; m=2; ell=2; T=500;
A=diag([.80,.62,.43]); B=0.4*randn(n,m);
C=randn(p,n); [C,~]=qr(C,0);
u=randn(m,T); x=zeros(n,T+1); y=zeros(p,T);
for k=1:T
    y(:,k)=C*x(:,k)+0.04*randn(p,1);
    x(:,k+1)=A*x(:,k)+B*u(:,k)+0.03*randn(n,1);
end

[~,~,P,Pbar,R,Rbar,~,Sigma_eps,Sigma_ebar,~,~,~,info] = ...
    predvarx_identify_strict_whitened(y,u,ell);

assert(norm(R'*P-eye(ell),'fro') < 1e-8, 'Eq. (28)/(29) de-normalization must yield R''P=I.');
assert(norm(R'*Pbar,'fro') < 1e-8, 'Eq. (34) complement must satisfy R''Pbar=0.');
assert(norm(Rbar'*P,'fro') < 1e-8, 'Eq. (34) complement must satisfy Rbar''P=0.');
assert(norm(Rbar'*Pbar-eye(p-ell),'fro') < 1e-8, 'Full dual basis must satisfy Rbar''Pbar=I.');
assert(info.whitening_applied,'strict whitening flag must be true');
assert(~info.complete_algorithm1,'this ablation must not claim complete Algorithm 1');
assert(info.normalized_orthogonality < 1e-10,'Pstar must be orthonormal in whitened coordinates');
assert(norm(info.R_raw'*info.P_raw-eye(ell),'fro') < 1e-8, ...
    'Raw Mo--Qin de-normalization must already satisfy R''P=I; no SVD alignment is allowed.');
assert(all(eig((Sigma_eps+Sigma_eps')/2) >= -1e-10), 'Dynamic covariance must be PSD.');
assert(all(eig((Sigma_ebar+Sigma_ebar')/2) >= -1e-10), 'Static covariance must be PSD.');
fprintf('PASS: whitening and direct rank-full s=1 de-normalization hold; complete Algorithm 1 is not claimed.\n');
end
