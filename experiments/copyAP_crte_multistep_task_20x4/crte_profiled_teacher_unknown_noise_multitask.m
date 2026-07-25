function [Ahat,Bhat,P,R,Sigma_eps,stats] = crte_profiled_teacher_unknown_noise_multitask(y,u,ell,tracked,opt)
%CRTE_PROFILED_TEACHER_UNKNOWN_NOISE_MULTITASK Profiled CRTE teacher.
%
% Decision object: a complete free subspace X in St(d,ell-q), together with
% metric parameter mu. Every candidate rebuilds (P,R), re-extracts z, and
% refits VARX before evaluating all teacher terms.
%
% The FWL task uses an actual future stack with task_horizon blocks:
%   [sqrt(w1)y_T(t+1); ...; sqrt(wH)y_T(t+H)].
% This changes the task Gram itself; it is independent of prediction_horizon.
%
% Unknown-noise variant: true Sigma_n is forbidden. The readout-noise term
% uses a two-fold cross-fitted one-step free-output residual covariance.
% This preserves the teacher structure but is not the known-Sigma_n formula.
%
% Search is a deterministic finite candidate study, not a claim of global
% optimization on the Stiefel/Grassmann manifold.

if nargin < 5 || isempty(opt), opt = struct(); end
def = struct('mu_grid',[0 0.25 0.5 0.75 1], ...
    'alpha',1,'beta',1,'omega',[], 'prediction_horizon',18, ...
    'task_horizon',1,'task_omega',[], ...
    'Ru',eye(size(u,1)),'val_fraction',0.25,'ridge',1e-8, ...
    'rank_tol',1e-9,'reach_tau',1e-10,'num_random_subspaces',30, ...
    'G',[],'seed',20260725);
fn=fieldnames(def);
for i=1:numel(fn)
    if ~isfield(opt,fn{i}) || isempty(opt.(fn{i})), opt.(fn{i})=def.(fn{i}); end
end

p=size(y,1); m=size(u,1); T=size(y,2); q=numel(tracked); r=ell-q;
assert(size(u,2)==T,'y and u sample counts differ.');
assert(r>0 && ell<=p,'Need 0 < ell-q and ell <= p.');
assert(all(opt.mu_grid>=0 & opt.mu_grid<=1),'mu must lie in [0,1].');
assert(isequal(size(opt.Ru),[m m]) && min(eig((opt.Ru+opt.Ru')/2))>0,'Ru must be SPD.');
Np=opt.prediction_horizon;
if isempty(opt.omega), omega=ones(1,Np)/Np; else, omega=opt.omega(:)'/sum(opt.omega); end
assert(numel(omega)==Np && all(omega>=0),'omega must be nonnegative and match horizon.');
Ht=opt.task_horizon;
assert(isscalar(Ht) && Ht==round(Ht) && Ht>=1,'task_horizon must be a positive integer.');
if isempty(opt.task_omega)
    task_omega=ones(1,Ht)/Ht;
else
    task_omega=opt.task_omega(:)'/sum(opt.task_omega);
end
assert(numel(task_omega)==Ht && all(task_omega>=0), ...
    'task_omega must be nonnegative and match task_horizon.');

nval=max(50,round(opt.val_fraction*T)); ntr=T-nval;
tr=1:ntr; va=ntr+1:T;
assert(ntr>max(Np,Ht)+ell+m+20,'Training segment too short.');

ybar=mean(y(:,tr),2); ubar=mean(u(:,tr),2);
yc=y-ybar; uc=u-ubar;
E=zeros(p,q); E(tracked,:)=eye(q); Nperp=null(E'); d=size(Nperp,2);
Yp=Nperp'*yc;
Sigma_perp=cov(Yp(:,tr)',1); Sigma_perp=(Sigma_perp+Sigma_perp')/2;
scale=max(trace(Sigma_perp)/d,1e-12); ridge=opt.ridge*max(scale,1)+1e-12;
Sigma_perp=Sigma_perp+ridge*eye(d);
if isempty(opt.G), G=eye(d); else, G=(opt.G+opt.G')/2; end
assert(isequal(size(G),[d d]) && min(eig(G))>0,'G must be d-by-d SPD.');
tauG=trace(G\Sigma_perp)/d;

% Exact OLS FWL residual maker on the common Ht-step task window.
% Ht=1 is exactly the copyAO one-step construction. For Ht>1, Tfuture is
% a genuine q*Ht-by-nobs future-task stack, not a prediction-horizon alias.
task_idx=tr(1):(tr(end)-Ht);
base=[yc(tracked,task_idx);uc(:,task_idx)];
W=Yp(:,task_idx);
Tfuture=zeros(q*Ht,numel(task_idx));
for h=1:Ht
    rr=(h-1)*q+(1:q);
    Tfuture(rr,:)=sqrt(task_omega(h))*yc(tracked,task_idx+h);
end
H0=base'*pinv(base*base')*base; H0=(H0+H0')/2;
M0=eye(size(H0))-H0; M0=(M0+M0')/2;
Z0=W*M0;
B_T=Z0*Z0';
TaskGram=Tfuture'*Tfuture;
A_T=Z0*((TaskGram+TaskGram')/2)*Z0';
B_T=(B_T+B_T')/2; A_T=(A_T+A_T')/2;
proj_idempotency=norm(H0*H0-H0,'fro');

% Strict Sec. 5.3 support parameterization: Z0=Us*Ss*Vs'.
% All candidates are generated as V0=Us*Ss^{-1}*Zeta, so they are born in
% range(B_T) and have V0'*B_T*V0=I before metric renormalization.
[Us,Ss,~]=svd(Z0,'econ');
svals=diag(Ss);
svd_tol=opt.rank_tol*max(max(svals),1);
fwl_rank=sum(svals>svd_tol);
if fwl_rank<r
    error('crte_profiled_teacher_unknown_noise:InsufficientFWLRank', ...
        'FWL support rank %d is smaller than requested ell_f=%d.',fwl_rank,r);
end
Us=Us(:,1:fwl_rank); svals=svals(1:fwl_rank);
Qsupport=Us*diag(1./svals);
Psupport=Us*Us';
support_B_identity_error=norm(Qsupport'*B_T*Qsupport-eye(fwl_rank),'fro');

% Unknown-noise object, estimated without using a candidate or true Sigma_n.
[Sigma_noise_proxy,noise_crossfit_info]=crossfit_free_residual_cov(Yp,uc,tr,ridge);

% Candidate pool generated only in the compact-SVD FWL support.
pool={}; ncan=0;
Syu=fixed_prediction_content(Yp,uc,tr,ridge);
NtrS=fro_norm(Syu); NtrA=fro_norm(A_T); NtrN=fro_norm(Sigma_noise_proxy);
for im=1:numel(opt.mu_grid)
    mu=opt.mu_grid(im); [Csqrt,~,Cmu]=metric_roots(Sigma_perp,G,tauG,mu,ridge);
    Ceta=Qsupport'*Cmu*Qsupport; Ceta=(Ceta+Ceta')/2;
    [~,Cei]=spd_roots(Ceta,ridge);
    for aa=[0 0.5 1]
        for bb=[0 0.5 1]
            Aeta=Qsupport'*(NtrS+aa*NtrA-bb*NtrN)*Qsupport;
            M=Cei*((Aeta+Aeta')/2)*Cei; M=(M+M')/2;
            [U,D]=eig(M); [~,ord]=sort(real(diag(D)),'descend'); Zeta=fix_sign(U(:,ord(1:r)));
            Vread=Qsupport*(Cei*Zeta);
            Vread=metric_normalize(Vread,Cmu,ridge);
            Vload=Cmu*Vread;
            X=Csqrt*Vread;
            ncan=ncan+1; pool{ncan}=make_candidate(mu,X,Vread,Vload,'spectral',aa,bb,Cmu); %#ok<AGROW>
        end
    end
end
rng(opt.seed,'twister');
for im=1:numel(opt.mu_grid)
    mu=opt.mu_grid(im); [Csqrt,~,Cmu]=metric_roots(Sigma_perp,G,tauG,mu,ridge);
    for k=1:opt.num_random_subspaces
        [Zeta,~]=qr(randn(fwl_rank,r),0); Zeta=fix_sign(Zeta);
        V0=Qsupport*Zeta;
        Vread=metric_normalize(V0,Cmu,ridge);
        Vload=Cmu*Vread;
        X=Csqrt*Vread;
        ncan=ncan+1; pool{ncan}=make_candidate(mu,X,Vread,Vload,'random',NaN,NaN,Cmu); %#ok<AGROW>
    end
end

rows=struct([]);
for i=1:ncan
    cand=pool{i};
    Vread=cand.Vread;
    Vload=cand.Vload;
    [Ac,Bc,Pc,Rc,~,det]=refit(yc,uc,E,Nperp,Vload,Vread,tr,ridge);

    % Complete profiled teacher term 1: multi-step residual in selected free DLVs.
    pred_term=multistep_free_residual(yc,uc,Rc,Ac,Bc,q,tr,Np,omega);

    % Term 2: exact multidimensional FWL trace on its declared domain.
    Bt=Vread'*B_T*Vread; Bt=(Bt+Bt')/2;
    At=Vread'*A_T*Vread; At=(At+At')/2;
    b_eigs=eig(Bt); fwl_valid=min(b_eigs)>opt.rank_tol*max(max(b_eigs),1);
    if fwl_valid
        task_term=real(trace(Bt\At))/size(W,2);
    else
        task_term=NaN;
    end

    % Term 3: unknown-noise readout proxy in the actual readout coordinates.
    noise_term=real(trace(Vread'*Sigma_noise_proxy*Vread));

    % Candidate-specific Ru-weighted finite-horizon authority.
    Wc=zeros(ell); Ru=(opt.Ru+opt.Ru')/2;
    for h=0:Np-1
        AhB=(Ac^h)*Bc;
        Wc=Wc+(AhB/Ru)*AhB';
    end
    reach_axes=diag(Wc(q+1:ell,q+1:ell));
    reach_min=min(reach_axes);

    teacher=pred_term-opt.alpha*task_term+opt.beta*noise_term;
    val=validation_multistep_nrmse(yc,uc,tracked,Ac,Bc,Pc,Rc,va,Np,omega);
    feasible=fwl_valid && reach_min>=opt.reach_tau && det.dual_error<1e-8 && det.rho<1.05 && isfinite(teacher);

    rows(i).mu=cand.mu; rows(i).source=cand.source;
    rows(i).initializer_alpha=cand.aa; rows(i).initializer_beta=cand.bb;
    rows(i).prediction_term=pred_term; rows(i).task_term=task_term;
    rows(i).noise_term=noise_term; rows(i).teacher_objective=teacher;
    rows(i).reach_min=reach_min; rows(i).fwl_min_eig=min(b_eigs);
    rows(i).support_residual=norm((eye(d)-Psupport)*Vread,'fro');
    rows(i).fwl_valid=fwl_valid; rows(i).validation_nrmse=val;
    rows(i).spectral_radius=det.rho; rows(i).dual_error=det.dual_error;
    rows(i).feasible=feasible; rows(i).X=cand.X; rows(i).Vread=Vread;
end
feas=[rows.feasible]; assert(any(feas),'No feasible profiled-teacher candidate.');
idx=find(feas); [~,k]=min([rows(idx).teacher_objective]); bestidx=idx(k); best=rows(bestidx);

% Refit the selected subspace on all offline data. X and mu remain frozen.
ybarF=mean(y,2); ubarF=mean(u,2); ycF=y-ybarF; ucF=u-ubarF;
YpF=Nperp'*ycF; SigF=cov(YpF',1); SigF=(SigF+SigF')/2+ridge*eye(d);
tauF=trace(G\SigF)/d; [~,~,CmuF]=metric_roots(SigF,G,tauF,best.mu,ridge);
Vread=metric_normalize(best.Vread,CmuF,ridge);
Vload=CmuF*Vread;
[Ahat,Bhat,P,R,Sigma_eps,det]=refit(ycF,ucF,E,Nperp,Vload,Vread,1:T,ridge);
Pbar=null(R'); Qfull=[P,Pbar]; assert(rcond(Qfull)>1e-12,'Complete dual basis ill-conditioned.');
Dmat=inv(Qfull'); Rbar=Dmat(:,ell+1:end);

stats.y_mean=ybarF; stats.u_mean=ubarF; stats.rows=rows; stats.best_index=bestidx;
stats.selected_mu=best.mu; stats.selected_source=best.source;
stats.selected_teacher_objective=best.teacher_objective;
stats.selected_prediction_term=best.prediction_term; stats.selected_task_term=best.task_term;
stats.selected_noise_term=best.noise_term; stats.selected_reach_min=best.reach_min;
stats.selected_validation_nrmse=best.validation_nrmse;
stats.num_candidates=ncan; stats.num_feasible=sum(feas);
stats.alpha=opt.alpha; stats.beta=opt.beta; stats.omega=omega; stats.horizon=Np;
stats.prediction_horizon=Np; stats.task_horizon=Ht; stats.task_omega=task_omega;
stats.task_future_rows=size(Tfuture,1); stats.task_num_observations=size(Tfuture,2);
stats.task_future_stack=Tfuture; stats.task_index=task_idx;
stats.uses_true_Sigma_n=false; stats.noise_object='two-fold forward-chaining blocked cross-fitted free-output residual covariance proxy';
stats.noise_crossfit_scheme=noise_crossfit_info.scheme;
stats.noise_crossfit_num_residuals=noise_crossfit_info.num_residuals;
stats.search_claim='finite candidate verification; no Stiefel global-optimum claim';
stats.H0_idempotency_error=proj_idempotency; stats.B_T=B_T; stats.A_T=A_T;
stats.Z0=Z0; stats.fwl_support_rank=fwl_rank; stats.fwl_support_tol=svd_tol;
stats.fwl_support_basis=Us; stats.fwl_support_singular_values=svals;
stats.support_B_identity_error=support_B_identity_error;
stats.max_candidate_support_residual=max([rows.support_residual]);
stats.Sigma_noise_proxy=Sigma_noise_proxy; stats.C_mu=CmuF;
stats.dual_error=norm(R'*P-eye(ell),'fro');
stats.tracked_right_error=norm(P*R'*E-E,'fro'); stats.tracked_left_error=norm(E'*P*R'-E','fro');
stats.Pi_idempotency_error=norm(P*R'*P*R'-P*R','fro');
stats.dual_errors_4piece=[stats.dual_error,norm(R'*Pbar,'fro'),norm(Rbar'*P,'fro'),norm(Rbar'*Pbar-eye(p-ell),'fro')];
stats.spectral_radius=max(abs(eig(Ahat))); stats.detail=det;
end

function c=make_candidate(mu,X,Vread,Vload,source,aa,bb,Cmu)
c=struct('mu',mu,'X',X,'Vread',Vread,'Vload',Vload, ...
    'source',source,'aa',aa,'bb',bb,'Cmu',Cmu);
end

function V=metric_normalize(V0,C,ridge)
Gram=V0'*C*V0; Gram=(Gram+Gram')/2;
[~,Gih]=spd_roots(Gram,ridge);
V=V0*Gih;
end

function [Gh,Gih]=spd_roots(G,ridge)
[U,D]=eig((G+G')/2); d=max(real(diag(D)),ridge);
Gh=U*diag(sqrt(d))*U'; Gih=U*diag(1./sqrt(d))*U';
end

function [Cs,Cih,C]=metric_roots(S,G,tau,mu,ridge)
C=(1-mu)*S+mu*tau*G; C=(C+C')/2+ridge*eye(size(C));
[U,D]=eig(C); d=max(real(diag(D)),ridge); Cs=U*diag(sqrt(d))*U'; Cih=U*diag(1./sqrt(d))*U';
end

function S=fixed_prediction_content(Yp,uc,tr,ridge)
Y0=Yp(:,tr(1:end-1)); Y1=Yp(:,tr(2:end)); U0=uc(:,tr(1:end-1)); Phi=[Y0;U0];
Th=(Phi*Phi'+ridge*eye(size(Phi,1)))\(Phi*Y1'); Yh=Th'*Phi; S=Yh*Yh'/size(Yh,2); S=(S+S')/2;
end

function [C,info]=crossfit_free_residual_cov(Yp,uc,tr,ridge)
% Two forward-chaining blocked folds. Each test block lies strictly after
% its fit block; no adjacent future samples leak back into the fitted model.
idx=tr(1:end-1); n=numel(idx);
cut1=max(2,floor(n/3)); cut2=max(cut1+1,floor(2*n/3));
fitBlocks={1:cut1,1:cut2};
testBlocks={cut1+1:cut2,cut2+1:n};
E=[];
for f=1:2
    fit=fitBlocks{f}; tst=testBlocks{f};
    assert(~isempty(fit) && ~isempty(tst) && max(fit)<min(tst), ...
        'Forward cross-fit blocks must be nonempty and ordered.');
    Phi=[Yp(:,idx(fit));uc(:,idx(fit))]; Tar=Yp(:,idx(fit)+1);
    Th=(Phi*Phi'+ridge*eye(size(Phi,1)))\(Phi*Tar');
    Phit=[Yp(:,idx(tst));uc(:,idx(tst))];
    Et=Yp(:,idx(tst)+1)-Th'*Phit;
    E=[E,Et]; %#ok<AGROW>
end
C=E*E'/size(E,2); C=(C+C')/2+ridge*eye(size(C));
info.scheme='two-fold forward-chaining blocked split';
info.num_residuals=size(E,2);
end

function [A,B,P,R,S,det]=refit(yc,uc,E,N,Vload,Vread,idx,ridge)
q=size(E,2); ell=q+size(Vload,2); m=size(uc,1); P=[E,N*Vload]; R=[E,N*Vread];
assert(norm(R'*P-eye(ell),'fro')<1e-7,'Metric dual failed.'); z=R'*yc(:,idx);
Phi=[z(:,1:end-1);uc(:,idx(1:end-1))]; Tar=z(:,2:end);
Th=(Phi*Phi'+ridge*eye(ell+m))\(Phi*Tar'); A=Th(1:ell,:)'; B=Th(ell+1:end,:)';
Ep=Tar-A*z(:,1:end-1)-B*uc(:,idx(1:end-1)); S=Ep*Ep'/size(Ep,2); S=(S+S')/2;
det.dual_error=norm(R'*P-eye(ell),'fro'); det.rho=max(abs(eig(A)));
end

function J=multistep_free_residual(yc,uc,R,A,B,q,tr,Np,omega)
z=R'*yc; ell=size(z,1); total=0; count=0; last=tr(end)-Np;
for t=tr(1):last
    zh=z(:,t);
    for h=1:Np
        zh=A*zh+B*uc(:,t+h-1);
        e=z(q+1:ell,t+h)-zh(q+1:ell); total=total+omega(h)*(e'*e); count=count+1;
    end
end
J=total/max(count/Np,1);
end

function v=validation_multistep_nrmse(yc,uc,tracked,A,B,P,R,va,Np,omega)
z=R'*yc; ss=zeros(numel(tracked),1); nn=0; last=va(end)-Np;
for t=va(1):last
    zh=z(:,t);
    for h=1:Np
        zh=A*zh+B*uc(:,t+h-1); e=yc(tracked,t+h)-P(tracked,:)*zh;
        ss=ss+omega(h)*(e.^2); nn=nn+1;
    end
end
rmse=sqrt(ss/max(nn/Np,1)); sc=std(yc(tracked,va),0,2); v=mean(rmse./max(sc,1e-8));
end

function A=fro_norm(A)
A=(A+A')/2; s=norm(A,'fro'); if s>1e-12, A=A/s; else, A=zeros(size(A)); end
end
function X=fix_sign(X)
for j=1:size(X,2), [~,i]=max(abs(X(:,j))); if X(i,j)<0, X(:,j)=-X(:,j); end, end
end
