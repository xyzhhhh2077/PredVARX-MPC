function test_split_control_free_ivr_varx
% RED/GREEN regression for split tracked/free extractor design.
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
rng(1707,'twister');

p = 9; m = 2; ell = 5; tracked = [1 2]; T = 650;
F = diag([0.91 0.77 0.58]);
G = [0.40 -0.15; 0.18 0.32; -0.08 0.24];
C = randn(p,3);
C(1,:) = [1 0 0];
C(2,:) = [0 1 0];
u = randn(m,T);
x = zeros(3,T+1);
y_clean = zeros(p,T);
for k = 1:T
    y_clean(:,k) = C*x(:,k);
    x(:,k+1) = F*x(:,k) + G*u(:,k) + 0.035*randn(3,1);
end

% Case 1: isotropic true sensor noise. The strict noise-optimal extractor
% must collapse to the orthogonal special case R=P.
sigma = 0.06;
Sigma_iso = sigma^2*eye(p);
y_iso = y_clean + sigma*randn(p,T);
[Ai,Bi,Pi,Ri,Si,sti] = split_control_free_ivr_varx( ...
    y_iso,u,ell,tracked,Sigma_iso);
E = zeros(p,numel(tracked)); E(tracked,:) = eye(numel(tracked));
assert(norm(Ri'*Pi-eye(ell),'fro') < 1e-9,'isotropic dual identity failed');
assert(norm(Pi*Ri'*E-E,'fro') < 1e-9,'isotropic right coverage failed');
assert(norm(E'*Pi*Ri'-E','fro') < 1e-9,'isotropic value preservation failed');
assert(norm(Ri-Pi,'fro') < 1e-8,'isotropic optimum must reduce to R=P');
assert(sti.free_noise_improvement >= -1e-10,'isotropic objective regressed');
assert(isequal(size(Ai),[ell ell]) && isequal(size(Bi),[ell m]));
assert(min(eig((Si+Si')/2)) > -1e-8,'isotropic Sigma is not PSD');

% Case 2: heteroscedastic and correlated true noise. The tracked columns
% remain exact while the free extractor should become genuinely oblique and
% strictly improve the true free-subspace noise objective over R=P.
scales = [0.025 0.030 0.045 0.080 0.130 0.210 0.330 0.480 0.700]';
Corr = eye(p);
Corr(3,4)=0.35; Corr(4,3)=0.35;
Corr(5,6)=-0.25; Corr(6,5)=-0.25;
D = diag(scales);
Sigma_het = D*Corr*D;
L = chol(Sigma_het,'lower');
y_het = y_clean + L*randn(p,T);
[Ah,Bh,Ph,Rh,Sh,sth] = split_control_free_ivr_varx( ...
    y_het,u,ell,tracked,Sigma_het);
assert(norm(Rh'*Ph-eye(ell),'fro') < 1e-8,'heteroscedastic dual identity failed');
assert(norm(Ph*Rh'*E-E,'fro') < 1e-8,'heteroscedastic right coverage failed');
assert(norm(E'*Ph*Rh'-E','fro') < 1e-8,'heteroscedastic value preservation failed');
assert(norm(Rh(:,1:numel(tracked))-E,'fro') < 1e-12,'tracked extractor columns changed');
assert(norm(Rh(:,numel(tracked)+1:end)-Ph(:,numel(tracked)+1:end),'fro') > 1e-4, ...
    'heteroscedastic free extractor did not become oblique');
assert(sth.free_noise_objective <= sth.free_noise_baseline + 1e-10, ...
    'true free-noise objective is worse than orthogonal baseline');
assert(sth.free_noise_improvement > 1e-6, ...
    'heteroscedastic case should strictly improve true free-noise objective');
assert(isequal(size(Ah),[ell ell]) && isequal(size(Bh),[ell m]));
assert(min(eig((Sh+Sh')/2)) > -1e-8,'heteroscedastic Sigma is not PSD');

fprintf(['PASS split extractor: iso ||R-P||=%.3e; het free-oblique=%.3e; ' ...
    'noise objective %.6g -> %.6g (improvement %.3e)\n'], ...
    norm(Ri-Pi,'fro'),sth.free_oblique_norm, ...
    sth.free_noise_baseline,sth.free_noise_objective,sth.free_noise_improvement);
end
