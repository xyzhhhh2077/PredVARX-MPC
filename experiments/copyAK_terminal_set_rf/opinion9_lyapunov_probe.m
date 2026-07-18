function out = opinion9_lyapunov_probe(y_off, u_off, ell, tracked, Sigma_n)
% OPINION9_LYAPUNOV_PROBE Numerical probes only — NOT recursive feasibility proof.
%
% Checks:
%   1) spectral radius of Ahat (Schur?)
%   2) discrete Lyapunov Pterm = dlyap(A',Qf) residual
%   3) simple level-set volume proxy trace/det of Pterm
%   4) one-step nominal decrease of V(z)=z'Pterm z along free dynamics z+=A z
%
% Explicit non-claims:
%   - No terminal invariant set X_f under constraints
%   - No chance-constraint recursive feasibility
%   - No closed-loop stability theorem under noise

out = struct();
out.claim_level = 'numerical_probe_only';
out.is_recursive_feasibility_proof = false;
out.is_stability_proof = false;

[Ahat,Bhat,P,R,Sigma_eps,stats] = split_control_free_ivr_varx( ...
    y_off, u_off, ell, tracked, Sigma_n);
out.dual_error = stats.dual_error;
out.left_error = stats.tracked_left_error;

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
        % Monte Carlo free dynamics decrease of V
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

% Controllability-ish proxy: rank of [B AB ...]
Co = Bhat;
for k=1:nz-1
    Co = [Co, Ahat^k*Bhat]; %#ok<AGROW>
end
out.ctrb_rank = rank(Co, 1e-8);
out.note = ['Opinion9 probe: Schur+Lyapunov decrease on free dynamics is necessary-style ' ...
    'numerical evidence only. Constraints, chance rows, and noise break automatic RF/stability.'];
end
