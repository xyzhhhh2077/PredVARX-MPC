function [Ahat,Bhat,P,R,Sigma_eps,stats] = split_empirical_cov_ivr_varx(y,u,ell,tracked,oblique_alpha)
% SPLIT_EMPIRICAL_COV_IVR_VARX
% Control-aware split construction WITHOUT oracle sensor-noise Sigma_n:
%   P = [E, Nperp*V]   fixed tracked axes + IVR free loading axes
%   R = [E, Nperp*W]   fixed tracked extractors + free dual under the
%                      empirical total free-coordinate covariance metric
%
% Free dual (alpha=1, default):
%   Yp = Nperp' * yc
%   C_emp = (Yp*Yp') / (T-1) + ridge
%   W = C_emp^{-1} V (V' C_emp^{-1} V)^{-1}
% which uniquely solves
%   min_W  tr(W' C_emp W)  s.t.  W' V = I
% under the empirical-total-covariance metric on free coordinates.
% This is NOT a sensor-noise optimum (true Sigma_n is unknown and unused).
% It is also NOT the PredVAR minimum-innovation-covariance theorem.
%
% oblique_alpha in [0,1]: W = V + alpha*(Wfull - V).
%   alpha=1 is the conditional optimum under the empirical free metric;
%   alpha<1 is only a numerical interpolation toward the orthogonal dual.
%
% Signature deliberately has no Sigma_n argument.

p = size(y,1);
m = size(u,1);
q = numel(tracked);
if nargin < 5 || isempty(oblique_alpha)
    oblique_alpha = 1;
end
assert(oblique_alpha >= 0 && oblique_alpha <= 1, ...
    'oblique_alpha must be in [0,1]; alpha=1 is empirical-metric conditional optimum.');
assert(ell >= q,'ell must be at least the number of tracked outputs.');
assert(size(y,2) == size(u,2),'y and u must have the same sample count.');

y_mean = mean(y,2);
u_mean = mean(u,2);
yc = y-y_mean;
uc = u-u_mean;

E = zeros(p,q);
E(tracked,:) = eye(q);
Nperp = null(E');
r = ell-q;
T = size(yc,2);

if r == 0
    V = zeros(size(Nperp,2),0);
    Wfree = V;
    C_emp = zeros(size(Nperp,2));
    ridge = 0;
    iter = 0;
    trace_hist = [];
    delta = 0;
else
    Yp = Nperp'*yc;
    Ylag = Yp(:,1:end-1);
    Ycur = Yp(:,2:end);
    [U0,~,~] = svd(Yp,'econ');
    V = U0(:,1:r);
    max_iter = 30;
    tol = 1e-4;
    trace_hist = nan(max_iter,1);
    for iter = 1:max_iter
        Xlag = V'*Ylag;
        Xcur = V'*Ycur;
        Aivr = (Xlag*Xlag'+1e-8*eye(r)) \ (Xlag*Xcur');
        Xpred = Aivr'*Xlag;
        trace_hist(iter) = trace(Xpred*Xpred')/max(size(Xpred,2),1);
        M = Ycur*Xpred'/(Xpred*Xpred'+1e-8*eye(r))*Xpred*Ycur';
        M = (M+M')/2;
        [Ui,Di] = eig(M);
        [~,ord] = sort(real(diag(Di)),'descend');
        Vnew = Ui(:,ord(1:r));
        for j = 1:r
            if V(:,j)'*Vnew(:,j) < 0
                Vnew(:,j) = -Vnew(:,j);
            end
        end
        delta = norm(Vnew*Vnew'-V*V','fro');
        V = Vnew;
        if iter > 1 && abs(trace_hist(iter)-trace_hist(iter-1)) < ...
                tol*max(abs(trace_hist(iter)),1) && delta < tol
            break
        end
    end
    trace_hist = trace_hist(1:iter);

    % Empirical total free-coordinate covariance metric (not Sigma_n).
    C_emp = (Yp*Yp')/max(T-1,1);
    C_emp = (C_emp+C_emp')/2;
    ridge = 1e-8*max(trace(C_emp)/max(size(C_emp,1),1),1) + 1e-12;
    C_reg = C_emp + ridge*eye(size(C_emp));
    SinvV = C_reg \ V;
    Gram = V'*SinvV;
    Gram = (Gram+Gram')/2;
    Wfull = SinvV / Gram;
    % alpha=1: empirical-metric conditional optimum; alpha<1: interpolation only.
    Wfree = V + oblique_alpha*(Wfull - V);
end

Pfree = Nperp*V;
Rfree = Nperp*Wfree;
P = [E,Pfree];
R = [E,Rfree];

% Coordinate-consistent coupled VARX refit after choosing R.
z = R'*yc;
zn = z(:,2:end);
zc = z(:,1:end-1);
ur = uc(:,1:end-1);
Phi = [zc;ur];
Theta = (Phi*Phi'+1e-8*eye(ell+m)) \ (Phi*zn');
Ahat = Theta(1:ell,:)';
Bhat = Theta(ell+1:end,:)';
Eps = zn-Ahat*zc-Bhat*ur;
% Zero-mean innovation model: conditional Gaussian ML denominator Nres.
% Do not use Nres-1/T-2 without an explicit mean-estimation argument.
Nres = size(Eps,2);
Sigma_eps = Eps*Eps'/max(Nres,1);
Sigma_eps = (Sigma_eps+Sigma_eps')/2;

stats.y_mean = y_mean;
stats.u_mean = u_mean;
stats.tracked_projection_error = norm(P*P'*E-E,'fro');
stats.tracked_right_error = norm(P*R'*E-E,'fro');
stats.tracked_left_error = norm(E'*P*R'-E','fro');
stats.tracked_column_error = norm(R(:,1:q)-E,'fro');
stats.dual_error = norm(R'*P-eye(ell),'fro');
stats.pr_asymmetry = norm(P*R'-(P*R')','fro');
stats.cond_dual_gram = cond(R'*P);
stats.free_oblique_norm = norm(Rfree-Pfree,'fro');
stats.reconstruction_residual = norm(yc-P*(R'*yc),'fro')/max(norm(yc,'fro'),eps);
stats.ivr_iter = iter;
stats.ivr_trace = trace_hist;
stats.ivr_subspace_delta = delta;
stats.oblique_alpha = oblique_alpha;
stats.uses_true_Sigma_n = false;
stats.N_residual = Nres;
stats.Sigma_eps_denom_primary = max(Nres,1);
stats.metric_name = 'empirical-total-covariance metric';
stats.C_emp_ridge = ridge;
if r == 0
    stats.free_emp_cov_baseline = 0;
    stats.free_emp_cov_objective = 0;
    stats.free_emp_cov_improvement = 0;
    stats.C_emp = [];
else
    % Objective under free-block empirical total covariance metric.
    stats.C_emp = C_emp;
    stats.free_emp_cov_baseline = trace(V'*C_reg*V);
    stats.free_emp_cov_objective = trace(Wfree'*C_reg*Wfree);
    stats.free_emp_cov_improvement = ...
        stats.free_emp_cov_baseline - stats.free_emp_cov_objective;
end
end
