function [Ahat,Bhat,P,R,Sigma_eps,stats] = fit_anchored_varx(y,u,Eanchor,ell,opt)
%FIT_ANCHORED_VARX Preserve a learned output anchor and fit a dual VARX.
if nargin < 5 || isempty(opt), opt=struct(); end
if ~isfield(opt,'ridge'), opt.ridge=1e-8; end
[p,T]=size(y); m=size(u,1); q=size(Eanchor,2); r=ell-q;
assert(size(u,2)==T,'y and u sample counts must match.');
assert(r>=0 && ell<=p && rank(Eanchor)==q,'Invalid anchor or latent dimension.');

ymean=mean(y,2); umean=mean(u,2);
yc=y-ymean; uc=u-umean;
scale=max(trace(yc*yc'/T)/p,1e-12);
ridge=opt.ridge*max(scale,1)+1e-12;
Cy=yc*yc'/T+ridge*eye(p); Cy=(Cy+Cy')/2;

% Normalize the anchor in the output covariance metric.
[Qe,~]=qr(Eanchor,0);
Ptask=Qe/real_sqrtm(Qe'*Cy*Qe,ridge);
Rtask=Cy*Ptask;
assert(norm(Rtask'*Ptask-eye(q),'fro')<1e-7,'Anchor normalization failed.');

if r>0
    Ntask=null(Rtask');
    % Predictable residual directions from the same old data only.
    Ylag=yc(:,1:end-1); Ycur=yc(:,2:end); Ulag=uc(:,1:end-1);
    Phi=[Ylag;Ulag];
    Theta=(Phi*Phi'+ridge*eye(p+m))\(Phi*Ycur');
    Ypred=Theta'*Phi;
    Sp=(Ypred*Ypred')/size(Ypred,2); Sp=(Sp+Sp')/2;
    Sf=Ntask'*Sp*Ntask; Cf=Ntask'*Cy*Ntask;
    [Uc,Dc]=eig((Cf+Cf')/2); dc=max(real(diag(Dc)),ridge);
    Cih=Uc*diag(1./sqrt(dc))*Uc';
    M=Cih*((Sf+Sf')/2)*Cih; M=(M+M')/2;
    [V,D]=eig(M); [~,ord]=sort(real(diag(D)),'descend');
    V=deterministic_sign(V(:,ord(1:r)));
    Pfree=Ntask*(Cih*V);
    Zfree=null(Ptask'); K=Zfree'*Pfree;
    assert(rank(K)==r,'Free loading is dependent on anchor.');
    Rfree=Zfree*K/(K'*K);
    P=[Ptask,Pfree]; R=[Rtask,Rfree];
else
    P=Ptask; R=Rtask;
end
assert(norm(R'*P-eye(ell),'fro')<1e-7,'Full dual identity failed.');

z=R'*yc; zn=z(:,2:end); zc=z(:,1:end-1); ur=uc(:,1:end-1);
Phi=[zc;ur]; Theta=(Phi*Phi'+ridge*eye(ell+m))\(Phi*zn');
Ahat=Theta(1:ell,:)'; Bhat=Theta(ell+1:end,:)';
Eps=zn-Ahat*zc-Bhat*ur;
Sigma_eps=Eps*Eps'/size(Eps,2); Sigma_eps=(Sigma_eps+Sigma_eps')/2;

stats=struct();
stats.y_mean=ymean; stats.u_mean=umean; stats.C_y=Cy;
stats.P_anchor=Ptask; stats.R_anchor=Rtask;
stats.dual_error=norm(R'*P-eye(ell),'fro');
stats.anchor_preservation_error=norm(P*R'*Ptask-Ptask,'fro');
stats.spectral_radius=max(abs(eig(Ahat)));
stats.uses_new_training_data=false;
end

function S=real_sqrtm(A,ridge)
A=(A+A')/2; [U,D]=eig(A); d=max(real(diag(D)),ridge);
S=U*diag(sqrt(d))*U';
end

function X=deterministic_sign(X)
for j=1:size(X,2)
    [~,i]=max(abs(X(:,j)));
    if X(i,j)<0, X(:,j)=-X(:,j); end
end
end
