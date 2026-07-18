function cal = calibrate_alpha_term(Ahat, Pterm, varargin)
% CALIBRATE_ALPHA_TERM Scan alpha levels for ellipsoid z'Pterm z <= alpha.
%
% Numerical calibration only — does NOT prove positive invariance of Xf
% under constraints, chance rows, or noise.
%
% Name-value:
%   'Nmc' (800) Monte Carlo unit/random samples
%   'quantiles' ([0.5 0.8 0.9 0.95 0.99]) of free V(z)=z'P z on unit ball
%   'seed' (19)
%
% Free-dynamics one-step map: if V(Az) <= V(z) for all samples (Schur+Lyap),
% any alpha>0 level set is positively invariant under z+=A z (unconstrained).
% That is NOT Xf under input/output/chance constraints.

p = inputParser;
addParameter(p, 'Nmc', 800, @(x)isnumeric(x)&&isscalar(x)&&x>=10);
addParameter(p, 'quantiles', [0.5 0.8 0.9 0.95 0.99], @(x)isnumeric(x)&&isvector(x));
addParameter(p, 'seed', 19, @(x)isnumeric(x)&&isscalar(x));
parse(p, varargin{:});
Nmc = p.Results.Nmc;
qs = p.Results.quantiles(:)';
seed = p.Results.seed;

cal = struct();
cal.claim_level = 'numerical_calibration_only';
cal.is_Xf_invariance_proof = false;

nz = size(Ahat,1);
Pt = (Pterm+Pterm')/2;
cal.Pterm_min_eig = min(eig(Pt));
cal.rho_A = max(abs(eig(Ahat)));

rng(seed, 'twister');
V0 = zeros(Nmc,1);
V1 = zeros(Nmc,1);
for i = 1:Nmc
    z = randn(nz,1); z = z / (norm(z)+eps);
    V0(i) = z'*Pt*z;
    z1 = Ahat*z;
    V1(i) = z1'*Pt*z1;
end
cal.V_unit_mean = mean(V0);
cal.V_unit_quantiles = qs;
cal.V_unit_qvals = quantile(V0, qs);
cal.V_decrease_ok_ratio = mean(V1 <= V0 + 1e-9);
cal.free_level_set_invariant_proxy = (cal.rho_A < 1-1e-10) && (cal.V_decrease_ok_ratio >= 1-1e-12);

% Recommended alpha_term candidates from free unit-ball V quantiles
% (scale-free directions). For CL, runner may rescale by typical ||z||^2.
cal.alpha_candidates = cal.V_unit_qvals;
cal.alpha_recommend = cal.V_unit_qvals(max(1, find(qs>=0.9,1,'first')));
if isempty(cal.alpha_recommend) || ~isfinite(cal.alpha_recommend)
    cal.alpha_recommend = median(V0);
end

cal.note = ['alpha from free-dynamics unit-ball V quantiles. ' ...
    'Not a certified terminal set under constraints/chance/noise.'];
end
