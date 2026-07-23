function [Ahat,Bhat,P,R,Sigma_eps,stats] = split_control_free_ivr_varx(y,u,ell,tracked,Sigma_n,varargin)
% SPLIT_CONTROL_FREE_IVR_VARX
% Control-aware split construction:
%   P = [E, Nperp*V]   fixed tracked axes + IVR free loading axes
%   R = [E, Nperp*W]   fixed tracked extractors + true-noise-optimal free dual
%
% Sigma_n is the known/declared sensor-noise covariance. The free extractor
% minimizes trace(Rf'*Sigma_n*Rf) subject to Rf'*Pf=I. The complete latent
% coordinate is then re-extracted and the coupled VARX model is refitted.
%
% Default (input_residualize=false):
%   IVR free-direction selection uses only projected output lags (not U).
%   Therefore the free span is predictable-in-output-history, not a proved
%   input-conditional PredVARX subspace. U enters only in the VARX step.
%
% Optional (input_residualize=true) - candidate scheme 1:
%   In Nperp coordinates, ridge-regress Yp on U, run IVR on the residual.
%   This is only a *candidate* input-conditionalization; it is NOT a proved
%   optimal input-conditional PredVARX free subspace. Flagged in
%   stats.ivr_input_conditional.
%
% Name-value options:
%   'input_residualize'  logical, default false

p = size(y,1);
m = size(u,1);
q = numel(tracked);
assert(ell >= q,'ell must be at least the number of tracked outputs.');
assert(size(y,2) == size(u,2),'y and u must have the same sample count.');
assert(isequal(size(Sigma_n),[p p]),'Sigma_n must be p-by-p.');
Sigma_n = (Sigma_n+Sigma_n')/2;
assert(min(eig(Sigma_n)) > 0,'Sigma_n must be positive definite.');

% --- optional name-value: input residualization for free IVR (default off)
input_residualize = false;
enforce_geometry = true;  % Opinion 1: geometry asserts default ON
if ~isempty(varargin)
    assert(mod(numel(varargin),2) == 0, ...
        'Optional arguments must be name-value pairs.');
    for k = 1:2:numel(varargin)
        name = varargin{k};
        val  = varargin{k+1};
        switch lower(char(name))
            case 'input_residualize'
                input_residualize = logical(val);
            case 'enforce_geometry'
                enforce_geometry = logical(val);
            otherwise
                error('split_control_free_ivr_varx:UnknownOption', ...
                    'Unknown option: %s', char(name));
        end
    end
end

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
    % Candidate scheme 1 (optional): residualize Yp on U before IVR.
    % Even when enabled this is only a candidate conditionalization -
    % not a proved optimal input-conditional free subspace.
    if input_residualize
        ridge_u = 1e-8*max(trace(uc*uc')/max(m,1),1);
        Bu = (uc*uc'+ridge_u*eye(m)) \ (uc*Yp.');  % m x d
        Yp_ivr = Yp - Bu.'*uc;
    else
        Yp_ivr = Yp;
    end
    Ylag = Yp_ivr(:,1:end-1);
    Ycur = Yp_ivr(:,2:end);
    [U0,~,~] = svd(Yp_ivr,'econ');
    V = U0(:,1:r);
    max_iter = 30;
    tol = 1e-4;
    trace_hist = nan(max_iter,1);
    % Opinion 3 - distinguish M vs tau = trace_hist(iter):
    %   Xpred  = Aivr'*Xlag : free-latent one-step predictor (coords z=V'y).
    %   tau_i  = trace_hist(iter) = tr(Xpred*Xpred')/T_lag
    %            = mean *latent* prediction energy of Xpred.
    %            NOT output predictable energy, and NOT tr(M)/anything.
    %   M      = Ycur*Xpred'*(S+lambda*I)^{-1}*Xpred*Ycur'  (lambda=1e-8)
    %            = regularized prediction *association* matrix used only to
    %            re-estimate free loading directions via eig(M).
    %            With ridge L_* = Ycur*Xpred'/(S+lambda*I), M equals the
    %            ridge objective-drop Gram, not Yhat*Yhat' in general:
    %              Yhat = L_* * Xpred
    %              Yhat*Yhat' = Ycur*Xpred'*(S+lambdaI)^{-1}*S*(S+lambdaI)^{-1}*Xpred*Ycur'
    %            Equality M = Yhat*Yhat' holds only when lambda=0 and S is
    %            invertible (or both use the same Moore-Penrose inverse).
    %            For the code lambda=1e-8 they are numerically close if S is
    %            well-conditioned (see tests/test_opinion03_M_vs_tau.m).
    for iter = 1:max_iter
        Xlag = V'*Ylag;
        Xcur = V'*Ycur;
        Aivr = (Xlag*Xlag'+1e-8*eye(r)) \ (Xlag*Xcur');
        Xpred = Aivr'*Xlag;
        % tau_i: latent prediction energy (free coords), not output energy
        trace_hist(iter) = trace(Xpred*Xpred')/max(size(Xpred,2),1);
        % M: regularized prediction association matrix (ridge Gram / drop)
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
        % Opinion 4: stop = numerical stationarity (trace + subspace delta),
        % not a proof of monotone convergence of the IVR objective.
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
% Opinion 5: residual covariance uses the conditional Gaussian ML scaling.
% The innovation mean is modeled as zero, so the primary denominator is
% Nres (not Nres-1).  OLS residual DOF is retained only as a diagnostic.
Nres = size(Eps,2);                       % = T-1 (one lag lost)
Gram_eps = Eps*Eps';
denom_ml  = max(Nres,1);                  % Gaussian ML / population MLE
denom_ols = max(Nres-(ell+m),1);          % OLS residual DOF after ell+m params
Sigma_eps     = (Gram_eps/denom_ml);       % primary, theory-backed zero-mean ML
Sigma_eps     = (Sigma_eps+Sigma_eps')/2;
Sigma_eps_ml  = Sigma_eps;
Sigma_eps_ols = (Gram_eps/denom_ols);
Sigma_eps_ols = (Sigma_eps_ols+Sigma_eps_ols')/2;

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
% Opinion 4: always export IVR iteration diagnostics (also when r==0).
% stop criterion above is numerical stationarity, not proven monotone.
stats.ivr_iter = iter;                 % iterations executed (0 if r==0)
stats.ivr_trace = trace_hist;          % per-iter predictable-energy trace
stats.ivr_subspace_delta = delta;      % final Frobenius projector change
stats.sensor_noise_covariance = Sigma_n;
% Explicit naming: true means candidate residualize-on-U IVR was used;
% it does NOT claim a proved optimal input-conditional free subspace.
stats.ivr_input_conditional = logical(input_residualize);
% Opinion 5 multi-denom stats (primary return is zero-mean Gaussian ML).
stats.N_residual = Nres;
stats.Sigma_eps_denom_ml  = denom_ml;
stats.Sigma_eps_denom_ols = denom_ols;
stats.Sigma_eps_denom_primary = denom_ml;
stats.Sigma_eps     = Sigma_eps;
stats.Sigma_eps_ml  = Sigma_eps_ml;
stats.Sigma_eps_ols = Sigma_eps_ols;
% Opinion 2: free-block noise objective is w.r.t. the *declared* Sigma_n
% (argument), not an estimated residual covariance. Report baseline R=P
% free block, optimized free dual, and improvement.
if r == 0
    stats.free_noise_baseline = 0;
    stats.free_noise_objective = 0;
else
    stats.free_noise_baseline = trace(Pfree'*Sigma_n*Pfree);
    stats.free_noise_objective = trace(Rfree'*Sigma_n*Rfree);
end
stats.free_noise_improvement = stats.free_noise_baseline-stats.free_noise_objective;
stats.free_noise_is_declared_Sigma_n = true;
% Degenerate case: when Sigma_n is isotropic on the free complement,
% the noise-optimal free dual collapses to R=P (orthogonal special case).
stats.isotropic_collapsed_to_R_eq_P = (norm(R-P,'fro') < 1e-8);
% Opinion 1: reconstruction-layer geometry enforcement.
% left / right / column / dual residuals must vanish numerically (tol=1e-8).
% Optional switch stats.enforce_geometry (default true). Also settable via
% name-value pair 'enforce_geometry', true|false.
stats.enforce_geometry = logical(enforce_geometry);
if stats.enforce_geometry
    geom_tol = 1e-8;
    assert(stats.tracked_left_error <= geom_tol, ...
        'geometry fail: tracked_left_error=%.3e > tol=%.1e', ...
        stats.tracked_left_error, geom_tol);
    assert(stats.tracked_right_error <= geom_tol, ...
        'geometry fail: tracked_right_error=%.3e > tol=%.1e', ...
        stats.tracked_right_error, geom_tol);
    assert(stats.tracked_column_error <= geom_tol, ...
        'geometry fail: tracked_column_error=%.3e > tol=%.1e', ...
        stats.tracked_column_error, geom_tol);
    assert(stats.dual_error <= geom_tol, ...
        'geometry fail: dual_error=%.3e > tol=%.1e', ...
        stats.dual_error, geom_tol);
end
end

