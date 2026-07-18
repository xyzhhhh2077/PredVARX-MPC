function cert = certify_input_original_alpha(uk, y, model, opt, tol)
% CERTIFY_INPUT_ORIGINAL_ALPHA Re-check a candidate first-step input under
% the ORIGINAL alpha_joint chance rows (constant-input extension over horizon).
%
% Certification here means:
%   max(A_ch * U_ext - b_ch) <= tol
% where U_ext stacks uk for all N steps (open-loop constant hold of uk).
%
% This is a model-based certificate for the ORIGINAL Boole allocation.
% It is NOT a recursive-feasibility or closed-loop stability proof.

if nargin < 5 || isempty(tol), tol = 1e-7; end
nu = size(model.B,2);
uk = uk(:);
assert(numel(uk)==nu, 'uk dimension mismatch');

z = model.R' * (y - model.y_mean);
opt_chk = opt; % keep original alpha_joint / N
[A_ch, b_ch, meta] = build_chance_rows(z, model, opt_chk);
N = meta.N;
U_ext = repmat(uk, N, 1);

% Also enforce input bounds on the extension
if isscalar(opt.u_min), lb = opt.u_min*ones(nu,1); else, lb = opt.u_min(:); end
if isscalar(opt.u_max), ub = opt.u_max*ones(nu,1); else, ub = opt.u_max(:); end
bound_ok = all(uk >= lb-1e-9) && all(uk <= ub+1e-9);

resid = A_ch*U_ext - b_ch;
max_viol = max(resid);
pass = bound_ok && (max_viol <= tol);

cert = struct();
cert.pass = logical(pass);
cert.max_cc_violation = max_viol;
cert.bound_ok = logical(bound_ok);
cert.tol = tol;
cert.alpha_joint = opt.alpha_joint;
cert.N = N;
cert.risk_each = meta.risk_each;
cert.U_ext = U_ext;
cert.note = 'original-alpha constant-hold certificate (model-based; not RF/stability)';
if pass
    cert.level = 'original_alpha_certified';
else
    cert.level = 'not_original_alpha';
end
end
