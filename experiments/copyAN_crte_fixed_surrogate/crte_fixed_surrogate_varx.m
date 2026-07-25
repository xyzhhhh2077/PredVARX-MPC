function [Ahat,Bhat,P,R,Sigma_eps,stats] = crte_fixed_surrogate_varx(y,u,ell,tracked,opt)
%CRTE_FIXED_SURROGATE_VARX Fixed CRTE spectral surrogate in the free complement.
% This function implements the paper's executable fixed surrogate, not the
% profiled teacher objective. Each candidate builds a complete metric dual,
% refits VARX, and is checked by task/noise/reachability gates.
%
% No true sensor-noise covariance is accepted. The noise term is the
% one-step free-output residual covariance estimated from training data.
%
% Optional struct fields (with defaults):
%   mu_grid            (1,:)  [0 0.25 0.5 0.75 1]
%   alpha_grid         (1,:)  [0 0.5 1]
%   beta_grid          (1,:)  [0 0.25 0.5]
%   val_fraction       (1,1)  0.25
%   ridge              (1,1)  1e-8
%   task_gate_fraction (1,1)  0.10
%   noise_gate_factor  (1,1)  2.0
%   reach_gate_fraction(1,1)  0.05
%   reach_horizon      (1,1)  18
%   G                  SPD or empty for eye(d)

if nargin < 5 || isempty(opt)
    opt = struct();
end
defaults = struct( ...
    'mu_grid',[0 0.25 0.5 0.75 1], ...
    'alpha_grid',[0 0.5 1], ...
    'beta_grid',[0 0.25 0.5], ...
    'val_fraction',0.25, ...
    'ridge',1e-8, ...
    'task_gate_fraction',0.10, ...
    'noise_gate_factor',2.0, ...
    'reach_gate_fraction',0.05, ...
    'reach_horizon',18, ...
    'G',[]);
fld = fieldnames(defaults);
for i=1:numel(fld)
    if ~isfield(opt,fld{i}) || isempty(opt.(fld{i}))
        opt.(fld{i}) = defaults.(fld{i});
    end
end
% q and r are defined further down using tracked and ell. The legacy early
% assert is intentionally deferred.

p = size(y,1);
m = size(u,1);
T = size(y,2);
q = numel(tracked);
r = ell-q;
assert(size(u,2)==T,'y and u must have the same sample count.');
assert(r>0,'CRTE copy requires at least one free latent coordinate.');
assert(all(opt.mu_grid>=0 & opt.mu_grid<=1),'mu_grid must be in [0,1].');
assert(opt.val_fraction>0 && opt.val_fraction<0.5,'val_fraction must be in (0,0.5).');

n_val = max(50,round(opt.val_fraction*T));
n_train = T-n_val;
assert(n_train>ell+m+20,'Training segment is too short.');
train_idx = 1:n_train;
val_idx = n_train+1:T;

E = zeros(p,q); E(tracked,:) = eye(q);
Nperp = null(E');
d = size(Nperp,2);
if isempty(opt.G)
    G = eye(d);
else
    G = (opt.G+opt.G')/2;
    assert(isequal(size(G),[d d]) && min(eig(G))>0,'G must be d-by-d SPD.');
end

y_mean = mean(y(:,train_idx),2);
u_mean = mean(u(:,train_idx),2);
yc = y-y_mean;
uc = u-u_mean;
Yp = Nperp'*yc;
Ytrain = Yp(:,train_idx);
Sigma_perp = cov(Ytrain',1);
Sigma_perp = (Sigma_perp+Sigma_perp')/2;
scale_perp = max(trace(Sigma_perp)/d,1e-12);
ridge = opt.ridge*max(scale_perp,1)+1e-12;
Sigma_perp = Sigma_perp+ridge*eye(d);
tau_G = trace(G\Sigma_perp)/d;

% Fixed prediction content S_yu from one-step free-output OLS.
Ylag = Yp(:,train_idx(1:end-1));
Ycur = Yp(:,train_idx(2:end));
Ulag = uc(:,train_idx(1:end-1));
Phi_free = [Ylag;Ulag];
Theta_free = (Phi_free*Phi_free'+ridge*eye(d+m))\(Phi_free*Ycur');
Ypred = Theta_free'*Phi_free;
S_yu = (Ypred*Ypred')/max(size(Ypred,2),1);
S_yu = (S_yu+S_yu')/2;
Efree = Ycur-Ypred;
Sigma_noise_proxy = (Efree*Efree')/max(size(Efree,2),1);
Sigma_noise_proxy = (Sigma_noise_proxy+Sigma_noise_proxy')/2+ridge*eye(d);

% Exact fixed-space FWL task matrices. Base regressors are past tracked
% outputs and centered inputs. Candidate regressors are current free-source
% coordinates. Tfuture contains one-step future tracked outputs.
base = [yc(tracked,train_idx(1:end-1)); uc(:,train_idx(1:end-1))];
W = Yp(:,train_idx(1:end-1));
Tfuture = yc(tracked,train_idx(2:end));
H0 = base'/(base*base'+ridge*eye(size(base,1)))*base;
M0 = eye(size(H0))-H0;
Qt = eye(q);
B_T = W*M0*W';
A_T = W*M0*Tfuture'*Qt*Tfuture*M0*W';
B_T = (B_T+B_T')/2;
A_T = (A_T+A_T')/2;

NtrS = safe_fro_normalize(S_yu);
NtrA = safe_fro_normalize(A_T);
NtrN = safe_fro_normalize(Sigma_noise_proxy);

rows = struct([]);
krow = 0;
for imu = 1:numel(opt.mu_grid)
    mu = opt.mu_grid(imu);
    Cmu = (1-mu)*Sigma_perp + mu*tau_G*G;
    Cmu = (Cmu+Cmu')/2+ridge*eye(d);
    [Uc,Dc] = eig(Cmu);
    dc = max(real(diag(Dc)),ridge);
    Csqrt = Uc*diag(sqrt(dc))*Uc';
    Cinvhalf = Uc*diag(1./sqrt(dc))*Uc';
    for ia = 1:numel(opt.alpha_grid)
        alpha = opt.alpha_grid(ia);
        for ib = 1:numel(opt.beta_grid)
            beta = opt.beta_grid(ib);
            Acrte = NtrS + alpha*NtrA - beta*NtrN;
            Acrte = (Acrte+Acrte')/2;
            % Symmetric whitening avoids a non-invariant gamma*I shift.
            Mgev = Cinvhalf*Acrte*Cinvhalf;
            Mgev = (Mgev+Mgev')/2;
            [Ux,Dx] = eig(Mgev);
            [evals,ord] = sort(real(diag(Dx)),'descend');
            X = Ux(:,ord(1:r));
            X = deterministic_sign(X);
            Vread = Cinvhalf*X;
            Vload = Csqrt*X;
            [Ac,Bc,Pc,Rc,Sc,detail] = refit_candidate( ...
                yc,uc,y_mean,u_mean,E,Nperp,Vload,Vread,train_idx,ridge);

            % Candidate-specific scalar FWL scores and finite-horizon
            % latent authority. c_x is the unit basis for each free DLV in
            % the s=1 companion state; the gate uses the worst selected axis.
            task_axis = zeros(r,1);
            noise_axis = zeros(r,1);
            reach_axis = zeros(r,1);
            Wc = zeros(ell);
            for h = 0:opt.reach_horizon-1
                AhB = (Ac^h)*Bc;
                Wc = Wc+AhB*AhB';
            end
            for j = 1:r
                vj = Vread(:,j);
                den = max(vj'*B_T*vj,ridge);
                task_axis(j) = max(real(vj'*A_T*vj/den)/size(W,2),0);
                noise_axis(j) = max(real(vj'*Sigma_noise_proxy*vj),0);
                cj = zeros(ell,1); cj(q+j)=1;
                reach_axis(j) = max(real(cj'*Wc*cj),0);
            end

            val_nrmse = validation_nrmse(yc,uc,y_mean,tracked,Ac,Bc,Pc,Rc,val_idx);
            krow = krow+1;
            rows(krow).mu = mu;
            rows(krow).alpha = alpha;
            rows(krow).beta = beta;
            rows(krow).eigenvalues = evals(1:r);
            rows(krow).task_score = mean(task_axis);
            rows(krow).task_min = min(task_axis);
            rows(krow).noise_score = mean(noise_axis);
            rows(krow).noise_max = max(noise_axis);
            rows(krow).reach_score = mean(reach_axis);
            rows(krow).reach_min = min(reach_axis);
            rows(krow).val_nrmse = val_nrmse;
            rows(krow).spectral_radius = max(abs(eig(Ac)));
            rows(krow).dual_error = detail.dual_error;
            rows(krow).P = Pc;
            rows(krow).R = Rc;
            rows(krow).A = Ac;
            rows(krow).B = Bc;
            rows(krow).Sigma_eps = Sc;
        end
    end
end

% Outer percentile-style gates scaled to this fixed candidate grid.
task_all = [rows.task_min];
noise_all = [rows.noise_max];
reach_all = [rows.reach_min];
task_gate = opt.task_gate_fraction*max(task_all);
noise_gate = opt.noise_gate_factor*median(noise_all);
reach_gate = opt.reach_gate_fraction*max(reach_all);
valid = task_all>=task_gate & noise_all<=noise_gate & reach_all>=reach_gate & ...
    [rows.dual_error]<1e-8 & [rows.spectral_radius]<1.05;
if ~any(valid)
    warning('crte_fixed_surrogate_varx:NoCandidatePassed', ...
        'No candidate passed all gates; using the minimum validation NRMSE among stable dual candidates.');
    valid = [rows.dual_error]<1e-8 & [rows.spectral_radius]<1.05;
end
assert(any(valid),'No stable dual-consistent CRTE candidate.');
valid_idx = find(valid);
[~,loc] = min([rows(valid_idx).val_nrmse]);
best_idx = valid_idx(loc);
best = rows(best_idx);

% Rebuild the selected hyperparameters on all offline data, as promised by
% the paper pipeline, then fully refit the dual basis and VARX quantities.
y_mean_full = mean(y,2);
u_mean_full = mean(u,2);
yc_full = y-y_mean_full;
uc_full = u-u_mean_full;
Yp_full = Nperp'*yc_full;
Sigma_full = cov(Yp_full',1);
Sigma_full = (Sigma_full+Sigma_full')/2+ridge*eye(d);
tau_full = trace(G\Sigma_full)/d;
Ylag_f = Yp_full(:,1:end-1); Ycur_f = Yp_full(:,2:end); Ulag_f = uc_full(:,1:end-1);
Phi_f = [Ylag_f;Ulag_f];
Theta_f = (Phi_f*Phi_f'+ridge*eye(d+m))\(Phi_f*Ycur_f');
Ypred_f = Theta_f'*Phi_f;
S_f = (Ypred_f*Ypred_f')/size(Ypred_f,2); S_f=(S_f+S_f')/2;
E_f = Ycur_f-Ypred_f;
N_f = (E_f*E_f')/size(E_f,2); N_f=(N_f+N_f')/2+ridge*eye(d);
base_f = [yc_full(tracked,1:end-1); uc_full(:,1:end-1)];
W_f = Yp_full(:,1:end-1); T_f = yc_full(tracked,2:end);
H_f = base_f'/(base_f*base_f'+ridge*eye(size(base_f,1)))*base_f;
M_f = eye(size(H_f))-H_f;
A_f = W_f*M_f*T_f'*Qt*T_f*M_f*W_f'; A_f=(A_f+A_f')/2;
C_f = (1-best.mu)*Sigma_full+best.mu*tau_full*G; C_f=(C_f+C_f')/2+ridge*eye(d);
[Uc,Dc]=eig(C_f); dc=max(real(diag(Dc)),ridge);
Csqrt=Uc*diag(sqrt(dc))*Uc'; Cinvhalf=Uc*diag(1./sqrt(dc))*Uc';
Acrte_f = safe_fro_normalize(S_f)+best.alpha*safe_fro_normalize(A_f)-best.beta*safe_fro_normalize(N_f);
Mgev=Cinvhalf*((Acrte_f+Acrte_f')/2)*Cinvhalf; Mgev=(Mgev+Mgev')/2;
[Ux,Dx]=eig(Mgev); [evals_f,ord]=sort(real(diag(Dx)),'descend');
X=deterministic_sign(Ux(:,ord(1:r)));
Vread=Cinvhalf*X; Vload=Csqrt*X;
[Ahat,Bhat,P,R,Sigma_eps,detail]=refit_candidate( ...
    yc_full,uc_full,y_mean_full,u_mean_full,E,Nperp,Vload,Vread,1:T,ridge);

% Complete static dual basis.
Pbar = null(R');
Qfull = [P,Pbar];
assert(rcond(Qfull)>=1e-12,'Complete dual basis is ill-conditioned.');
DualMat = inv(Qfull');
Rbar = DualMat(:,ell+1:end);

stats.y_mean = y_mean_full;
stats.u_mean = u_mean_full;
stats.selected_mu = best.mu;
stats.selected_alpha = best.alpha;
stats.selected_beta = best.beta;
stats.selected_index = best_idx;
stats.selected_validation_nrmse = best.val_nrmse;
stats.selected_eigenvalues = evals_f(1:r);
stats.candidates = rows;
stats.valid_candidates = valid;
stats.task_gate = task_gate;
stats.noise_gate = noise_gate;
stats.reach_gate = reach_gate;
stats.uses_true_Sigma_n = false;
stats.noise_object = 'one-step free-output residual covariance proxy';
stats.objective = 'fixed CRTE spectral surrogate, not profiled teacher objective';
stats.C_mu = C_f;
stats.tau_G = tau_full;
stats.S_yu = S_f;
stats.A_T = A_f;
stats.Sigma_noise_proxy = N_f;
stats.dual_error = norm(R'*P-eye(ell),'fro');
stats.tracked_right_error = norm(P*R'*E-E,'fro');
stats.tracked_left_error = norm(E'*P*R'-E','fro');
stats.tracked_column_error = norm(R(:,1:q)-E,'fro');
stats.Pi_idempotency_err = norm(P*R'*P*R'-P*R','fro');
stats.dual_errors_4piece = [stats.dual_error,norm(R'*Pbar,'fro'), ...
    norm(Rbar'*P,'fro'),norm(Rbar'*Pbar-eye(p-ell),'fro')];
stats.dual_basis_completion = max(stats.dual_errors_4piece)<1e-9;
stats.pr_asymmetry = norm(P*R'-(P*R')','fro');
stats.reconstruction_residual = norm(yc_full-P*(R'*yc_full),'fro')/max(norm(yc_full,'fro'),eps);
stats.spectral_radius = max(abs(eig(Ahat)));
stats.ridge = ridge;
stats.train_count = n_train;
stats.validation_count = n_val;
stats.detail = detail;
end

function A = safe_fro_normalize(A)
A = (A+A')/2;
scale = norm(A,'fro');
if scale <= 1e-12
    A = zeros(size(A));
else
    A = A/scale;
end
end

function X = deterministic_sign(X)
for j=1:size(X,2)
    [~,i]=max(abs(X(:,j)));
    if X(i,j)<0, X(:,j)=-X(:,j); end
end
end

function [Ahat,Bhat,P,R,Sigma_eps,detail] = refit_candidate(yc,uc,y_mean,u_mean,E,Nperp,Vload,Vread,idx,ridge) %#ok<INUSD>
q=size(E,2); ell=q+size(Vload,2); m=size(uc,1);
P=[E,Nperp*Vload]; R=[E,Nperp*Vread];
assert(norm(R'*P-eye(ell),'fro')<1e-7,'Metric dual construction failed.');
z=R'*yc(:,idx);
uidx=idx;
zn=z(:,2:end); zc=z(:,1:end-1); ur=uc(:,uidx(1:end-1));
Phi=[zc;ur];
Theta=(Phi*Phi'+ridge*eye(ell+m))\(Phi*zn');
Ahat=Theta(1:ell,:)'; Bhat=Theta(ell+1:end,:)';
Eps=zn-Ahat*zc-Bhat*ur;
Sigma_eps=(Eps*Eps')/max(size(Eps,2),1); Sigma_eps=(Sigma_eps+Sigma_eps')/2;
detail.dual_error=norm(R'*P-eye(ell),'fro');
detail.spectral_radius=max(abs(eig(Ahat)));
detail.y_mean=y_mean; detail.u_mean=u_mean;
end

function score = validation_nrmse(yc,uc,y_mean,tracked,A,B,P,R,val_idx) %#ok<INUSD>
% One-step validation NRMSE on controlled outputs using frozen train fit.
yhat=zeros(numel(tracked),numel(val_idx)-1);
ytrue=yc(tracked,val_idx(2:end));
for k=1:numel(val_idx)-1
    z=R'*yc(:,val_idx(k));
    zn=A*z+B*uc(:,val_idx(k));
    yhat(:,k)=P(tracked,:)*zn;
end
rmse=sqrt(mean((ytrue-yhat).^2,2));
scale=std(ytrue,0,2); scale=max(scale,1e-8);
score=mean(rmse./scale);
end
