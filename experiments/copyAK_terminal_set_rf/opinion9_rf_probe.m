function out = opinion9_rf_probe(y_off, u_off, ell, tracked, Sigma_n)
% OPINION9_RF_PROBE Terminal-set / RF numerical skeleton (NOT a theorem).
%
% Documents the minimal RF condition chain and runs numerical probes:
%   (C1) rho(Ahat)<1  (Schur open-loop latent model)
%   (C2) Pterm = dlyap(A',Qf) exists, residual tiny
%   (C3) V(z)=z'Pterm z decreases along free dynamics z+=A z
%   (C4) terminal law kappa(z)=u_mean (or 0 in centered v) keeps Xf
%        positively invariant  --- NOT PROVED here under constraints
%   (C5) stage-feasible => predicted z_N in Xf --- only enforced when
%        use_terminal_set=true as a coarse QP constraint (box approx)
%   (C6) chance/noise RF under joint risk --- NOT PROVED
%
% Explicit non-claims:
%   - No recursive feasibility theorem under chance constraints
%   - No closed-loop stability theorem under process/sensor noise
%   - Spectral-box outer approx is NOT exact ellipsoid Xf

out = struct();
out.claim_level = 'numerical_probe_only';
out.is_recursive_feasibility_proof = false;
out.is_stability_proof = false;
out.is_Xf_positive_invariance_proof = false;

[Ahat,Bhat,P,R,Sigma_eps,stats] = split_control_free_ivr_varx( ...
    y_off, u_off, ell, tracked, Sigma_n);
out.dual_error = stats.dual_error;
out.left_error = stats.tracked_left_error;
out.Ahat = Ahat;
out.Bhat = Bhat;
out.P = P;
out.R = R;
out.Sigma_eps = Sigma_eps;

ev = eig(Ahat);
out.A_eigs = ev;
out.rho_A = max(abs(ev));
out.is_schur = out.rho_A < 1 - 1e-10;

nz = size(Ahat,1);
Qf = eye(nz);
out.Qf = Qf;
out.Pterm = [];
out.lyap_residual_fro = NaN;
out.V_decrease_ok_ratio = NaN;
out.dlyap_ok = false;

if out.is_schur
    try
        Pt = dlyap(Ahat', Qf);
        Pt = (Pt+Pt')/2;
        out.Pterm = Pt;
        out.dlyap_ok = true;
        Res = Ahat'*Pt*Ahat - Pt + Qf;
        out.lyap_residual_fro = norm(Res,'fro');
        out.Pterm_min_eig = min(eig(Pt));
        out.Pterm_trace = trace(Pt);
        out.Pterm_logdet = sum(log(max(eig(Pt), realmin)));
        rng(9,'twister');
        Nmc = 400;
        ok = 0;
        for i=1:Nmc
            z = randn(nz,1); z = z/(norm(z)+eps);
            V0 = z'*Pt*z;
            z1 = Ahat*z;
            V1 = z1'*Pt*z1;
            if V1 <= V0 - 0.5*(z'*Qf*z) + 1e-9
                ok = ok+1;
            end
        end
        out.V_decrease_ok_ratio = ok/Nmc;
    catch ME
        out.dlyap_ok = false;
        out.dlyap_error = ME.message;
    end
end

% Controllability rank proxy
Co = Bhat;
for k=1:nz-1
    Co = [Co, Ahat^k*Bhat]; %#ok<AGROW>
end
out.ctrb_rank = rank(Co, 1e-8);

% Alpha calibration on free dynamics
if out.dlyap_ok
    out.alpha_cal = calibrate_alpha_term(Ahat, out.Pterm, 'Nmc', 800, 'seed', 19);
else
    out.alpha_cal = struct('alpha_recommend', NaN, 'alpha_candidates', [], ...
        'free_level_set_invariant_proxy', false);
end

% RF condition chain status (boolean probes, not proofs)
out.C1_schur = logical(out.is_schur);
out.C2_dlyap = logical(out.dlyap_ok) && isfinite(out.lyap_residual_fro) ...
    && out.lyap_residual_fro < 1e-8;
out.C3_Vdec_free = isfinite(out.V_decrease_ok_ratio) && out.V_decrease_ok_ratio >= 0.999;
out.C4_Xf_invar_constrained = false;   % NOT proved
out.C5_stage_to_Xf_theorem = false;    % NOT proved (only optional QP rows)
out.C6_chance_RF = false;              % NOT proved

out.rf_chain = { ...
    'C1 rho(Ahat)<1: PROBED (Schur numeric)'; ...
    'C2 Pterm=dlyap(A'',Qf), residual small: PROBED'; ...
    'C3 V-decrease on free z+=Az: PROBED'; ...
    'C4 terminal law kappa keeps Xf +invariant under constraints: UNPROVED'; ...
    'C5 stage feasible => z_N in Xf (exact ellipsoid + constraints): UNPROVED'; ...
    'C6 chance-constraint recursive feasibility under noise: UNPROVED'};

out.note = ['Opinion9 RF probe: free Schur+Lyapunov is necessary-style ' ...
    'numeric evidence only. Constraints, chance rows, and noise break ' ...
    'automatic RF/stability. Terminal box is a coarse QP probe, not Xf theorem.'];
end
