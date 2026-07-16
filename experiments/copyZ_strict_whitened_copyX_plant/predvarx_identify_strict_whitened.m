function [Ahat,Bhat,P,Pbar,R,Rbar,G,Sigma_eps,Sigma_ebar,F,H,lambda,info] = predvarx_identify_strict_whitened(y,u,ell)
% PREDVARX_IDENTIFY_STRICT_WHITENED
% Rank-full, first-order (s=1) whitening/IVR/de-normalization ablation based
% on the Mo--Qin normalized-space construction, followed by a separately
% fitted centered VARX regression.  This is NOT a complete Algorithm-1
% implementation: the paper supports general VAR order s, uses an
% increase-only IVR stopping rule, constructs the full complement including
% the covariance nullspace, and computes its prescribed covariance outputs.
% Here p=r for this noisy synthetic data, and the implemented de-normalization
% is P=UD^(1/2)P*, R=UD^(-1/2)P*. No QR or post-hoc SVD alignment is applied.

[p,T]=size(y); m=size(u,1);
y_mean=mean(y,2); u_mean=mean(u,2); yc=y-y_mean; uc=u-u_mean;
Sigma_y=(yc*yc')/T; [U,D]=eig((Sigma_y+Sigma_y')/2);
[d,idx]=sort(real(diag(D)),'descend'); U=real(U(:,idx));
d=max(d,1e-10); Dsqrt=diag(sqrt(d)); Dinv=diag(1./sqrt(d));
Ystar=Dinv*U'*yc; N=T-1; Y0=Ystar(:,1:N); Ys=Ystar(:,2:T);
Pi0=Y0'/(Y0*Y0'+1e-9*eye(p))*Y0;
Pstar=top_eigs(Ys*Pi0*Ys'/N,ell);
trace_prev=-inf;
for iter=1:30
    z=Pstar'*Ystar; V=z(:,1:N)'; Vs=z(:,2:T)';
    Bvar=(V'*V+1e-9*eye(ell))\(V'*Vs); Vpred=V*Bvar;
    trace_now=real(trace((Vpred'*Vpred)/N));
    Pi_pred=Vpred/(Vpred'*Vpred+1e-9*eye(ell))*Vpred';
    Pnext=top_eigs(Ys*Pi_pred*Ys'/N,ell);
    if iter>1 && abs(trace_now-trace_prev)<1e-6*max(1,abs(trace_now))
        Pstar=Pnext; break;
    end
    Pstar=Pnext; trace_prev=trace_now;
end
Pbarstar=orth_complement(Pstar);
P=U*Dsqrt*Pstar; R=U*Dinv*Pstar;
Pbar=U*Dsqrt*Pbarstar; Rbar=U*Dinv*Pbarstar;

% PredVARX extension: centered latent coordinates and explicit u channel.
% Theorem 1 concerns the innovations of the original VAR realization.  This
% experiment does not claim that a separately fitted VARX residual inherits
% its sample cross-covariance cancellation, so it retains the direct
% whitening/de-normalization realization rather than applying an invalid
% post-fit N update.
z=R'*yc; zn=z(:,2:end); zc=z(:,1:end-1); ur=uc(:,1:end-1);
Theta=([zc;ur]*[zc;ur]'+1e-9*eye(ell+m))\([zc;ur]*zn');
Ahat=Theta(1:ell,:)'; Bhat=Theta(ell+1:end,:)';
E=zn-Ahat*zc-Bhat*ur; Sigma_eps=(E*E')/max(N-1,1); Sigma_eps=(Sigma_eps+Sigma_eps')/2;
ebar=Rbar'*yc(:,2:end); Sigma_ebar=(ebar*ebar')/max(N-1,1); Sigma_ebar=(Sigma_ebar+Sigma_ebar')/2;
Sigma_ebar_eps_final=(ebar*E')/max(N-1,1);
G=eye(ell); F=Ahat; H=Bhat; lambda=0.95;
info.whitening_applied=true; info.realization='rank_full_s1_moqin_style_whitening_ablation';
info.complete_algorithm1=false;
info.algorithm1_limitations=['s=1 only; convergence does not require strict trace increase; ', ...
    'rank-null complement and prescribed Algorithm-1 covariance outputs are not implemented; ', ...
    'final dynamics/covariance are separately fitted VARX quantities'];
info.normalized_orthogonality=norm(Pstar'*Pstar-eye(ell),'fro');
info.y_mean=y_mean; info.u_mean=u_mean; info.U=U; info.D=sqrt(d); info.Pstar=Pstar;
info.P_raw=P; info.R_raw=R; info.cross_cov_after=Sigma_ebar_eps_final; info.ivr_iterations=iter; info.predicted_trace=trace_now;
info.dual_errors=[norm(R'*P-eye(ell),'fro'),norm(R'*Pbar,'fro'),norm(Rbar'*P,'fro'),norm(Rbar'*Pbar-eye(p-ell),'fro')];
end

function V=top_eigs(M,ell)
M=(M+M')/2; [U,D]=eig(M); [~,idx]=sort(real(diag(D)),'descend'); V=real(U(:,idx(1:ell)));
end
function V=orth_complement(P)
[U,~,~]=svd(P','econ'); %#ok<ASGLU>
V=null(P');
end
