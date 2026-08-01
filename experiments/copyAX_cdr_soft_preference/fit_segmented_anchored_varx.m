function [Ahat,Bhat,P,R,Sigma_eps,stats] = fit_segmented_anchored_varx(y,u,run_id,Eanchor,ell,opt)
%FIT_SEGMENTED_ANCHORED_VARX copyAU anchored CRTE without cross-run lags.
% CDR alignment contract: u(:,k) drives y(:,k+1) (standard causality).
% Transitions pair (z(:,k), u(:,k)) -> z(:,k+1), so ur=uc(:,valid).
if nargin < 6 || isempty(opt), opt=struct(); end
if ~isfield(opt,'ridge'), opt.ridge=1e-8; end
if ~isfield(opt,'mu'), opt.mu=0.10; end
if ~isfield(opt,'ntr_epsilon'), opt.ntr_epsilon=10e-6; end
[p,T]=size(y); m=size(u,1); q=size(Eanchor,2); r=ell-q;
assert(size(u,2)==T && numel(run_id)==T,'Sample counts must match.');
assert(r>=0 && ell<=p && rank(Eanchor)==q,'Invalid anchor or latent dimension.');
assert(opt.mu>0 && opt.mu<1,'mu must be an interior metric blend.');
assert(opt.ntr_epsilon>0,'ntr_epsilon must be positive.');
run_id=run_id(:)'; valid=find(run_id(1:end-1)==run_id(2:end));
assert(~isempty(valid),'No within-run transitions are available.');

ymean=mean(y,2); umean=mean(u,2);
yc=y-ymean; uc=u-umean;
scale=max(trace(yc*yc'/T)/p,1e-12);
ridge=opt.ridge*max(scale,1)+1e-12;
Cy=yc*yc'/T+ridge*eye(p); Cy=(Cy+Cy')/2;
[Qe,~]=qr(Eanchor,0); Ptask=deterministic_sign(Qe); Rtask=Ptask;

if r>0
    Ntask=null(Ptask'); d=size(Ntask,2); Yp=Ntask'*yc;
    Sigma_perp=Yp*Yp'/T; Sigma_perp=(Sigma_perp+Sigma_perp')/2+ridge*eye(d);
    tau_G=trace(Sigma_perp)/d;
    Cmu=(1-opt.mu)*Sigma_perp+opt.mu*tau_G*eye(d);
    Cmu=(Cmu+Cmu')/2+ridge*eye(d);

    Ylag=Yp(:,valid); Ycur=Yp(:,valid+1); Ulag=uc(:,valid);
    Phi_free=[Ylag;Ulag];
    Theta_free=(Phi_free*Phi_free'+ridge*eye(d+m))\(Phi_free*Ycur');
    Ypred=Theta_free'*Phi_free;
    S_yu=Ypred*Ypred'/size(Ypred,2); S_yu=(S_yu+S_yu')/2;
    Efree=Ycur-Ypred;
    C_n=Efree*Efree'/size(Efree,2); C_n=(C_n+C_n')/2+ridge*eye(d);

    base=[Ptask'*yc(:,valid);uc(:,valid)]; W=Yp(:,valid);
    Tfuture=Ptask'*yc(:,valid+1);
    H0=base'/(base*base'+ridge*eye(size(base,1)))*base;
    M0=eye(size(H0))-H0;
    task_gram=Tfuture'*Tfuture;
    A_T=W*M0*task_gram*M0*W'; A_T=(A_T+A_T')/2;

    NtrS=paper_ntr(S_yu,opt.ntr_epsilon);
    NtrA=paper_ntr(A_T,opt.ntr_epsilon);
    NtrN=paper_ntr(C_n,opt.ntr_epsilon);
    Acrte=NtrS+NtrA-NtrN; Acrte=(Acrte+Acrte')/2;
    [Uc,Dc]=eig(Cmu); dc=max(real(diag(Dc)),ridge);
    Csqrt=Uc*diag(sqrt(dc))*Uc'; Cinvhalf=Uc*diag(1./sqrt(dc))*Uc';
    Mgev=Cinvhalf*Acrte*Cinvhalf; Mgev=(Mgev+Mgev')/2;
    [V,D]=eig(Mgev); [evals,ord]=sort(real(diag(D)),'descend');
    V=deterministic_sign(V(:,ord(1:r)));
    Pfree=Ntask*(Csqrt*V); Rfree=Ntask*(Cinvhalf*V);
    P=[Ptask,Pfree]; R=[Rtask,Rfree];
else
    P=Ptask; R=Rtask; d=0; Sigma_perp=[]; tau_G=[]; Cmu=[];
    S_yu=[]; A_T=[]; C_n=[]; NtrS=[]; NtrA=[]; NtrN=[]; evals=[];
end
assert(norm(R'*P-eye(ell),'fro')<1e-7,'Full dual identity failed.');

z=R'*yc; zn=z(:,valid+1); zc=z(:,valid); ur=uc(:,valid);
Phi=[zc;ur]; Theta=(Phi*Phi'+ridge*eye(ell+m))\(Phi*zn');
Ahat=Theta(1:ell,:)'; Bhat=Theta(ell+1:end,:)';
Eps=zn-Ahat*zc-Bhat*ur;
Sigma_eps=Eps*Eps'/size(Eps,2); Sigma_eps=(Sigma_eps+Sigma_eps')/2;

stats=struct('y_mean',ymean,'u_mean',umean,'C_y',Cy, ...
    'P_anchor',Ptask,'R_anchor',Rtask, ...
    'dual_error',norm(R'*P-eye(ell),'fro'), ...
    'anchor_preservation_error',max(norm(P*R'*Ptask-Ptask,'fro'), ...
    norm(Ptask'*P*R'-Ptask','fro')), ...
    'spectral_radius',max(abs(eig(Ahat))),'uses_new_training_data',true, ...
    'selected_mu',opt.mu,'ntr_epsilon',opt.ntr_epsilon, ...
    'ntr_mode','paper_trace_normalize_by_mean_trace', ...
    'ntr_formula','A/max(abs(trace(A))/d,10e-6)', ...
    'free_dimension',d,'Sigma_perp',Sigma_perp,'tau_G',tau_G,'C_mu',Cmu, ...
    'S_yu',S_yu,'A_T',A_T,'C_n',C_n,'Ntr_S_yu',NtrS, ...
    'Ntr_A_T',NtrA,'Ntr_C_n',NtrN, ...
    'crte_eigenvalues',evals(1:min(r,numel(evals))), ...
    'transition_count',numel(valid),'segment_count',numel(unique(run_id)));
end

function A=paper_ntr(A,epsilon_ntr)
A=(A+A')/2; d=size(A,1);
A=A/max(abs(trace(A))/d,epsilon_ntr);
end

function X=deterministic_sign(X)
for j=1:size(X,2)
    [~,i]=max(abs(X(:,j)));
    if X(i,j)<0, X(:,j)=-X(:,j); end
end
end
