function [Ahat,Bhat,P,R,Sigma_eps,stats] = control_aware_iterative_ivr_varx(y,u,ell,tracked)
% CONTROL_AWARE_ITERATIVE_IVR_VARX
% Control-aware reduced model with fixed tracked axes and iterative IVR
% selection of the remaining predictable directions.
%
% The tracked axes are kept exactly as in copyU. IVR is applied only in
% their orthogonal complement so the control-facing geometry is unchanged.

p = size(y,1);
m = size(u,1);
q = numel(tracked);
if ell < q
    error('ell must be at least the number of tracked outputs.');
end
if size(y,2) ~= size(u,2)
    error('y and u must have the same number of samples.');
end

y_mean = mean(y,2);
u_mean = mean(u,2);
yc = y - y_mean;
uc = u - u_mean;

E = zeros(p,q);
E(tracked,:) = eye(q);
Nperp = null(E');
r = ell-q;
if r == 0
    V = zeros(size(Nperp,2),0);
else
    Yp = Nperp' * yc;
    Ylag = Yp(:,1:end-1);
    Ycur = Yp(:,2:end);

    % Initial predictable subspace from one SVD.
    [U0,~,~] = svd(Yp,'econ');
    V = U0(:,1:min(r,size(U0,2)));

    max_iter = 30;
    tol = 1e-4;
    trace_hist = nan(max_iter,1);
    for iter = 1:max_iter
        Xlag = V' * Ylag;
        Xcur = V' * Ycur;
        Aivr = (Xlag*Xlag' + 1e-8*eye(r)) \ (Xlag*Xcur');
        Xpred = Aivr' * Xlag;
        trace_hist(iter) = trace(Xpred*Xpred') / max(size(Xpred,2),1);

        % Re-estimate directions that best explain the predictable part.
        M = Ycur * Xpred' / (Xpred*Xpred' + 1e-8*eye(r)) * Xpred * Ycur';
        M = (M+M')/2;
        [Ui,Di] = eig(M);
        [~,ord] = sort(real(diag(Di)),'descend');
        Vnew = Ui(:,ord(1:r));

        % Align signs to avoid artificial oscillation between iterations.
        for j = 1:r
            if V(:,j)'*Vnew(:,j) < 0
                Vnew(:,j) = -Vnew(:,j);
            end
        end
        delta = norm(Vnew*Vnew' - V*V','fro');
        V = Vnew;
        if iter > 1 && abs(trace_hist(iter)-trace_hist(iter-1)) < tol*max(abs(trace_hist(iter)),1) && delta < tol
            break
        end
    end
    trace_hist = trace_hist(1:iter);
end

P = [E, Nperp*V];
R = P;
z = R' * yc;
zn = z(:,2:end);
zc = z(:,1:end-1);
ur = uc(:,1:end-1);
Phi = [zc;ur];
Theta = (Phi*Phi' + 1e-8*eye(ell+m)) \ (Phi*zn');
Ahat = Theta(1:ell,:)';
Bhat = Theta(ell+1:end,:)';
Eps = zn-Ahat*zc-Bhat*ur;
Sigma_eps = (Eps*Eps')/max(size(Eps,2)-1,1);
Sigma_eps = (Sigma_eps+Sigma_eps')/2;

stats.y_mean = y_mean;
stats.u_mean = u_mean;
stats.tracked_projection_error = norm(P*P'*E-E,'fro');
stats.reconstruction_residual = norm(yc-P*(P'*yc),'fro')/max(norm(yc,'fro'),eps);
stats.ivr_iter = exist('iter','var')*iter;
if r == 0
    stats.ivr_trace = [];
    stats.ivr_subspace_delta = 0;
else
    stats.ivr_trace = trace_hist;
    stats.ivr_subspace_delta = delta;
end
stats.ivr_max_iter = 30;
stats.ivr_tol = 1e-4;
end
