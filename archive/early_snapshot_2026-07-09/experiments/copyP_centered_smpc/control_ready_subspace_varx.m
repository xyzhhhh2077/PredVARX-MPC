function [Ahat,Bhat,P,R,Sigma_eps,stats] = control_ready_subspace_varx(y,u,ell)
% CONTROL_READY_SUBSPACE_VARX Centered full-order subspace VARX for MPC.
% Uses a numerically orthonormal output subspace; this is the orthogonal
% special case of the dual-basis realization (R=P, R'P=I).
y_mean=mean(y,2); u_mean=mean(u,2); yc=y-y_mean; uc=u-u_mean;
[U,~,~]=svd(yc,'econ'); P=U(:,1:ell); R=P;
z=R'*yc; zn=z(:,2:end); zc=z(:,1:end-1); ur=uc(:,1:end-1);
Phi=[zc;ur]; Theta=(Phi*Phi'+1e-8*eye(ell+size(u,1)))\(Phi*zn');
Ahat=Theta(1:ell,:)'; Bhat=Theta(ell+1:end,:)';
E=zn-Ahat*zc-Bhat*ur; Sigma_eps=(E*E')/max(size(E,2)-1,1); Sigma_eps=(Sigma_eps+Sigma_eps')/2;
stats.y_mean=y_mean; stats.u_mean=u_mean; stats.subspace_residual=norm(yc-P*(P'*yc),'fro')/norm(yc,'fro');
end
