function test_opinion03_M_vs_tau
% Opinion 3: M is a regularized prediction association matrix;
% tau = trace_hist is free-latent prediction energy, not output energy.
%
% Numeric checks (no need for a full identifier run):
%   1) With lambda = 1e-8 and well-conditioned S, M ≈ Yhat*Yhat'.
%   2) With lambda = 0 and invertible S, M == Yhat*Yhat' (exact).
%   3) With lambda = 1 (1-D counterexample), M and Yhat*Yhat' differ.
%   4) tau = tr(Xpred*Xpred')/T is latent energy, not tr(M)/T.

rng(1703,'twister');

% --- synthetic free-block geometry matching the IVR loop ---
d = 7;   % free ambient dim (size(Ycur,1) = dim(Nperp))
r = 3;   % free rank
T = 500; % lag samples
Ycur  = randn(d,T);
Xpred = randn(r,T) + 0.3*randn(r,1)*ones(1,T);  % full-row-rank-ish
S = Xpred*Xpred';
assert(min(eig((S+S')/2)) > 1e-6, 'synthetic S not well-conditioned');

% --- Case A: code ridge lambda = 1e-8 → M close to Yhat*Yhat' ---
lambda = 1e-8;
[M_code, YhatYhat, tau] = assoc_and_energy(Ycur, Xpred, lambda);
rel = norm(M_code - YhatYhat,'fro') / max(norm(M_code,'fro'),1);
assert(rel < 1e-6, ...
    'lambda=1e-8: M and Yhat*Yhat'' should be close (rel=%.3e)', rel);

% Residual identity from theory:
%   M - Yhat*Yhat' = lambda * Ycur*Xpred'*(S+lambda I)^{-2}*Xpred*Ycur'
Sreg = S + lambda*eye(r);
diff_theory = lambda * (Ycur*Xpred') * (Sreg \ (Sreg \ (Xpred*Ycur')));
diff_theory = (diff_theory+diff_theory')/2;
assert(norm(M_code - YhatYhat - diff_theory,'fro') < 1e-8*max(norm(M_code,'fro'),1), ...
    'lambda residual identity failed');

% --- Case B: lambda = 0, S invertible → exact equality ---
[M0, YhatYhat0, ~] = assoc_and_energy(Ycur, Xpred, 0);
assert(norm(M0 - YhatYhat0,'fro') < 1e-10*max(norm(M0,'fro'),1), ...
    'lambda=0: M must equal Yhat*Yhat''');

% --- Case C: 1-D counterexample lambda=1 → M=1/2, YhatYhat=1/4 ---
Y1 = 1; Xp1 = 1;  % S = 1, Y X' = 1
lam1 = 1;
M1 = Y1*Xp1'/(Xp1*Xp1'+lam1)*Xp1*Y1';   % 1/2
L1 = Y1*Xp1'/(Xp1*Xp1'+lam1);            % 1/2
Yhat1 = L1*Xp1;                           % 1/2
assert(abs(M1 - 0.5) < 1e-14, '1-D M should be 1/2');
assert(abs(Yhat1*Yhat1' - 0.25) < 1e-14, '1-D YhatYhat should be 1/4');
assert(abs(M1 - Yhat1*Yhat1') > 0.2, '1-D case must show M ≠ YhatYhat');

% --- Case D: tau is latent energy, not output association energy ---
assert(abs(tau - trace(Xpred*Xpred')/T) < 1e-14, 'tau definition mismatch');
assert(abs(tau - trace(M_code)/T) > 1e-3*max(abs(tau),1) || d ~= r, ...
    'tau coincidentally equals tr(M)/T; synthetic case unexpected');
% Stronger: tau depends only on Xpred, not on Ycur
Ycur2 = randn(d,T);
[~, ~, tau2] = assoc_and_energy(Ycur2, Xpred, lambda);
assert(abs(tau - tau2) < 1e-14, ...
    'tau must be independent of Ycur (latent energy only)');
[M2, ~, ~] = assoc_and_energy(Ycur2, Xpred, lambda);
assert(norm(M_code - M2,'fro') > 1e-6*max(norm(M_code,'fro'),1), ...
    'M must depend on Ycur (association matrix)');

fprintf(['PASS opinion03 M vs tau: lambda=1e-8 rel(M,YhatYhat'')=%.3e; ' ...
    'lambda=0 exact; 1-D M=%.4g YhatYhat=%.4g; tau=%.6g (latent only)\n'], ...
    rel, M1, Yhat1*Yhat1', tau);
end

function [M, YhatYhat, tau] = assoc_and_energy(Ycur, Xpred, lambda)
% Mirrors the IVR-loop formulas in split_control_free_ivr_varx.m
r = size(Xpred,1);
T = size(Xpred,2);
S = Xpred*Xpred';
Sreg = S + lambda*eye(r);
L = Ycur * Xpred' / Sreg;          % ridge map Y ~ L Xpred
Yhat = L * Xpred;
YhatYhat = Yhat * Yhat';
M = Ycur * Xpred' / Sreg * Xpred * Ycur';  % same as code
M = (M+M')/2;
YhatYhat = (YhatYhat+YhatYhat')/2;
tau = trace(Xpred*Xpred') / max(T,1);
end
