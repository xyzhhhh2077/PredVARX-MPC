function [Ahat,Bhat,P,R,Sigma_eps,stats] = split_control_free_ivr_varx(y,u,ell,tracked,Sigma_n)
% SPLIT_CONTROL_FREE_IVR_VARX
% Control-aware split construction:
%   P = [E, Nperp*V]   fixed tracked axes + IVR free loading axes
%   R = [E, Nperp*W]   fixed tracked extractors + true-noise-optimal free dual
%
% Sigma_n is the known/declared sensor-noise covariance. The free extractor
% minimizes trace(Rf'*Sigma_n*Rf) subject to Rf'*Pf=I. The complete latent
% coordinate is then re-extracted and the coupled VARX model is refitted.

p = size(y,1);
m = size(u,1);
q = numel(tracked);
assert(ell >= q,'ell must be at least the number of tracked outputs.');
assert(size(y,2) == size(u,2),'y and u must have the same sample count.');
assert(isequal(size(Sigma_n),[p p]),'Sigma_n must be p-by-p.');
Sigma_n = (Sigma_n+Sigma_n')/2;
assert(min(eig(Sigma_n)) > 0,'Sigma_n must be positive definite.');

y_mean = mean(y,2);
u_mean = mean(u,2);
yc = y-y_mean;
uc = u-u_mean;

E = zeros(p,q);
E(tracked,:) = eye(q);
Nperp = null(E');
r = ell-q;

if r == 0
    V = zeros(size(Nperp,2),0);
    Wfree = V;
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

    Sigma_perp = Nperp'*Sigma_n*Nperp;
    Sigma_perp = (Sigma_perp+Sigma_perp')/2;
    ridge = 1e-12*max(trace(Sigma_perp)/max(size(Sigma_perp,1),1),1);
    SinvV = (Sigma_perp+ridge*eye(size(Sigma_perp)))\V;
    Gram = V'*SinvV;
    Wfree = SinvV/((Gram+Gram')/2);
end

Pfree = Nperp*V;
Rfree = Nperp*Wfree;
P = [E,Pfree];
R = [E,Rfree];

% Coordinate-consistent coupled VARX refit.
z = R'*yc;
zn = z(:,2:end);
zc = z(:,1:end-1);
ur = uc(:,1:end-1);
Phi = [zc;ur];
Theta = (Phi*Phi'+1e-8*eye(ell+m)) \ (Phi*zn');
Ahat = Theta(1:ell,:)';
Bhat = Theta(ell+1:end,:)';
Eps = zn-Ahat*zc-Bhat*ur;
Sigma_eps = Eps*Eps'/max(size(Eps,2)-1,1);
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
stats.sensor_noise_covariance = Sigma_n;
if r == 0
    stats.free_noise_baseline = 0;
    stats.free_noise_objective = 0;
else
    stats.free_noise_baseline = trace(Pfree'*Sigma_n*Pfree);
    stats.free_noise_objective = trace(Rfree'*Sigma_n*Rfree);
end
stats.free_noise_improvement = stats.free_noise_baseline-stats.free_noise_objective;
end
