function [Ahat,Bhat,P,R,Sigma_eps,stats] = control_aware_subspace_varx(y,u,ell,tracked)
% CONTROL_AWARE_SUBSPACE_VARX Reduced orthogonal model retaining tracked axes.
% P=[E_tracked, Q_perp] makes each tracked output exactly representable.
p=size(y,1); m=size(u,1); q=numel(tracked);
if ell<q, error('ell must be at least the number of tracked outputs.'); end
y_mean=mean(y,2); u_mean=mean(u,2); yc=y-y_mean; uc=u-u_mean;
E=zeros(p,q); E(tracked,:)=eye(q);
Nperp=null(E');
[U,~,~]=svd(Nperp'*yc,'econ');
P=[E, Nperp*U(:,1:ell-q)]; R=P;
z=R'*yc; zn=z(:,2:end); zc=z(:,1:end-1); ur=uc(:,1:end-1);
Phi=[zc;ur]; Theta=(Phi*Phi'+1e-8*eye(ell+m))\(Phi*zn');
Ahat=Theta(1:ell,:)'; Bhat=Theta(ell+1:end,:)';
Eps=zn-Ahat*zc-Bhat*ur; Sigma_eps=(Eps*Eps')/max(size(Eps,2)-1,1); Sigma_eps=(Sigma_eps+Sigma_eps')/2;
stats.y_mean=y_mean; stats.u_mean=u_mean; stats.tracked_projection_error=norm(P*P'*E-E,'fro');
stats.reconstruction_residual=norm(yc-P*(P'*yc),'fro')/norm(yc,'fro');
end
