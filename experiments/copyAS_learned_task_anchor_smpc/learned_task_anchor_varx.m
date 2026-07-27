function [Ahat,Bhat,P,R,Sigma_eps,stats] = learned_task_anchor_varx(y,u,ell,q,opt)
%LEARNED_TASK_ANCHOR_VARX Learn a task-output anchor, then fit dual VARX.
% Extension beyond CRTE draft: no tracked-output indices are accepted.
% Task information must enter through task_reference or Qy.
% True sensor-noise covariance is never accepted; residual covariance is a proxy.

if nargin < 5 || isempty(opt), opt = struct(); end
defaults = struct('task_reference',[],'Qy',[],'Ru',[], ...
    'anchor_weights',[1 0.6 0.5 0.4], ... % task,predict,reach,noise
    'mu_grid',[0 0.25 0.5 0.75 1], ...
    'val_fraction',0.25,'ridge',1e-8,'reach_horizon',12);
f = fieldnames(defaults);
for i=1:numel(f)
    if ~isfield(opt,f{i}) || isempty(opt.(f{i})), opt.(f{i})=defaults.(f{i}); end
end

[p,T] = size(y); m = size(u,1); r = ell-q;
assert(size(u,2)==T,'y and u sample counts must match.');
assert(q>=1 && r>=1 && ell<=p,'Require 1 <= q < ell <= p.');
assert(numel(opt.anchor_weights)==4,'anchor_weights must have four entries.');
assert(~isempty(opt.task_reference) || ~isempty(opt.Qy), ...
    'A task_reference or Qy is required; task-free anchor learning is undefined.');
if ~isempty(opt.task_reference)
    assert(isequal(size(opt.task_reference),[p T]),'task_reference must be p-by-T.');
end

nval=max(50,round(opt.val_fraction*T)); ntrain=T-nval;
assert(ntrain>p+m+20,'Training segment too short.');
train=1:ntrain; val=ntrain+1:T;
ym=mean(y(:,train),2); um=mean(u(:,train),2);
yc=y-ym; uc=u-um;
scale=max(trace((yc(:,train)*yc(:,train)')/ntrain)/p,1e-12);
ridge=opt.ridge*max(scale,1)+1e-12;
Cy=(yc(:,train)*yc(:,train)')/ntrain + ridge*eye(p);
Cy=(Cy+Cy')/2;

% Full-output one-step prediction content and unknown-noise residual proxy.
Ylag=yc(:,train(1:end-1)); Ycur=yc(:,train(2:end));
Ulag=uc(:,train(1:end-1)); Phi=[Ylag;Ulag];
Theta=(Phi*Phi'+ridge*eye(p+m))\(Phi*Ycur');
Ypred=Theta'*Phi;
Syu=(Ypred*Ypred')/size(Ypred,2); Syu=(Syu+Syu')/2;
Eres=Ycur-Ypred;
Cn=(Eres*Eres')/size(Eres,2); Cn=(Cn+Cn')/2+ridge*eye(p);
Ay=Theta(1:p,:)'; By=Theta(p+1:end,:)';

% Task matrix from full-output reference or declared output cost.
Mtask=zeros(p);
if ~isempty(opt.task_reference)
    Rt=opt.task_reference(:,train); Rt=Rt-mean(Rt,2);
    Mtask=(Rt*Rt')/size(Rt,2);
end
if ~isempty(opt.Qy)
    Qy=(opt.Qy+opt.Qy')/2;
    assert(isequal(size(Qy),[p p]) && min(eig(Qy))>=-1e-10,'Qy must be p-by-p PSD.');
    Mtask=Mtask+Qy;
end
Mtask=(Mtask+Mtask')/2;
assert(norm(Mtask,'fro')>1e-12,'Task matrix is numerically zero.');

% Full-output finite-horizon input authority from the fitted output VARX.
if isempty(opt.Ru), Ru=eye(m); else, Ru=(opt.Ru+opt.Ru')/2; end
assert(isequal(size(Ru),[m m]) && min(eig(Ru))>0,'Ru must be SPD.');
Rui=Ru\eye(m); Wy=zeros(p);
for h=0:opt.reach_horizon-1
    Gh=(Ay^h)*By; Wy=Wy+Gh*Rui*Gh';
end
Wy=(Wy+Wy')/2;

w=opt.anchor_weights(:);
ME=w(1)*fro_norm(Mtask)+w(2)*fro_norm(Syu)+w(3)*fro_norm(Wy)-w(4)*fro_norm(Cn);
ME=(ME+ME')/2;

% Generalized eigensystem for learned task anchor Etask' Cy Etask = I.
[Uc,Dc]=eig(Cy); dc=max(real(diag(Dc)),ridge);
Cih=Uc*diag(1./sqrt(dc))*Uc';
M=Cih*ME*Cih; M=(M+M')/2;
[Z,Dz]=eig(M); [anchor_evals,ord]=sort(real(diag(Dz)),'descend');
Z=deterministic_sign(Z(:,ord));
Etask=Cih*Z(:,1:q);
% Re-normalize defensively.
Etask=Etask/real_sqrtm(Etask'*Cy*Etask,ridge);
Ptask=Etask; Rtask=Cy*Etask;
assert(norm(Rtask'*Ptask-eye(q),'fro')<1e-7,'Learned task anchor pair is not dual.');
Ntask=null(Rtask');
assert(size(Ntask,2)==p-q,'Task-anchor complement has wrong dimension.');

% Fixed free-space score; candidate mu only changes free metric.
Sfree=Ntask'*Syu*Ntask; Tfree=Ntask'*Mtask*Ntask;
Nfree=Ntask'*Cn*Ntask; Wfree=Ntask'*Wy*Ntask;
Afree=fro_norm(Sfree)+w(1)*fro_norm(Tfree)+w(3)*fro_norm(Wfree)-w(4)*fro_norm(Nfree);
Afree=(Afree+Afree')/2;
Sigmafree=Ntask'*Cy*Ntask; Sigmafree=(Sigmafree+Sigmafree')/2;
Gfree=eye(p-q); tau=trace(Gfree\Sigmafree)/(p-q);

rows=struct([]);
for im=1:numel(opt.mu_grid)
    mu=opt.mu_grid(im);
    Cmu=(1-mu)*Sigmafree+mu*tau*Gfree+ridge*eye(p-q); Cmu=(Cmu+Cmu')/2;
    [U,D]=eig(Cmu); dd=max(real(diag(D)),ridge); Cfi=U*diag(1./sqrt(dd))*U';
    H=Cfi*Afree*Cfi; H=(H+H')/2;
    [V,Dv]=eig(H); [ev,o]=sort(real(diag(Dv)),'descend');
    V=deterministic_sign(V(:,o(1:r)));
    Pfree=Ntask*(Cfi\V); % Cmu^(1/2)*V
    Zfree=null(Ptask'); K=Zfree'*Pfree;
    assert(rank(K)>=r,'Free loading is not independent of task anchor.');
    Rfree=Zfree*K/(K'*K);
    Pc=[Ptask,Pfree]; Rc=[Rtask,Rfree];
    [Ac,Bc,Sc]=refit(yc,uc,Pc,Rc,train,ridge);
    rows(im).mu=mu; rows(im).eigenvalues=ev(1:r);
    rows(im).A=Ac; rows(im).B=Bc; rows(im).P=Pc; rows(im).R=Rc; rows(im).Sigma_eps=Sc;
    rows(im).dual_error=norm(Rc'*Pc-eye(ell),'fro');
    rows(im).spectral_radius=max(abs(eig(Ac)));
    rows(im).val_nrmse=validation_nrmse(yc,uc,Ac,Bc,Pc,Rc,val);
end
valid=[rows.dual_error]<1e-7 & [rows.spectral_radius]<1.10;
assert(any(valid),'No stable dual-consistent learned-anchor candidate.');
vi=find(valid); [~,j]=min([rows(vi).val_nrmse]); best=rows(vi(j));
P=best.P; R=best.R;

stats.y_mean=mean(y,2); stats.u_mean=mean(u,2);
% Refit final model on all centered data with selected P/R.
ycf=y-stats.y_mean; ucf=u-stats.u_mean;
[Ahat,Bhat,Sigma_eps]=refit(ycf,ucf,P,R,1:T,ridge);
stats.E_task_anchor=Etask; stats.P_task_anchor=Ptask; stats.R_task_anchor=Rtask;
stats.tracked_indices=[]; stats.anchor_eigenvalues=anchor_evals;
if q<p, stats.anchor_eigengap=max(anchor_evals(q)-anchor_evals(q+1),0); else, stats.anchor_eigengap=NaN; end
stats.anchor_score_matrix=ME; stats.C_y=Cy; stats.M_task=Mtask;
stats.S_yu_output=Syu; stats.W_y=Wy; stats.Sigma_noise_proxy=Cn;
stats.anchor_weights=w; stats.selected_mu=best.mu;
stats.selected_validation_nrmse=best.val_nrmse; stats.candidates=rows;
stats.dual_error=norm(R'*P-eye(ell),'fro');
stats.task_preservation_error=norm(P*R'*Ptask-Ptask,'fro');
stats.spectral_radius=max(abs(eig(Ahat)));
stats.uses_true_Sigma_n=false;
stats.method='learned task anchor extension; not CRTE draft algorithm';
stats.selection='full-output validation NRMSE after stable dual gate';
end

function A=fro_norm(A)
A=(A+A')/2; n=norm(A,'fro');
if n<=1e-12, A=zeros(size(A)); else, A=A/n; end
end

function X=deterministic_sign(X)
for j=1:size(X,2)
    [~,i]=max(abs(X(:,j))); if X(i,j)<0, X(:,j)=-X(:,j); end
end
end

function S=real_sqrtm(A,ridge)
A=(A+A')/2; [U,D]=eig(A); d=max(real(diag(D)),ridge); S=U*diag(sqrt(d))*U';
end

function [A,B,S]=refit(yc,uc,P,R,idx,ridge)
z=R'*yc(:,idx); zn=z(:,2:end); zc=z(:,1:end-1); ur=uc(:,idx(1:end-1));
Phi=[zc;ur]; Theta=(Phi*Phi'+ridge*eye(size(Phi,1)))\(Phi*zn');
ell=size(P,2); A=Theta(1:ell,:)'; B=Theta(ell+1:end,:)';
E=zn-A*zc-B*ur; S=(E*E')/size(E,2); S=(S+S')/2;
end

function score=validation_nrmse(yc,uc,A,B,P,R,idx)
ytrue=yc(:,idx(2:end)); yhat=zeros(size(ytrue));
for k=1:numel(idx)-1
    z=R'*yc(:,idx(k)); yhat(:,k)=P*(A*z+B*uc(:,idx(k)));
end
rmse=sqrt(mean((ytrue-yhat).^2,2)); scale=max(std(ytrue,0,2),1e-8); score=mean(rmse./scale);
end
