function info = fallback_certify_step(y, model, opt, uk)
% FALLBACK_CERTIFY_STEP  Heuristic fallback bookkeeping after QP failure.
%
% After applying a saturated mean (or other heuristic) input, optionally
% compute a one-step *deterministic mean* prediction residual against the
% hard quality bound opt.h. This is a simplified diagnostic only.
%
% IMPORTANT — chance-constraint status:
%   This function ALWAYS returns
%       info.cc_cert_level = 'uncertified_fallback'
%   A mean prediction that happens to lie under y_max does NOT restore the
%   Boole/Gaussian chance-constraint guarantee of the original QP. Do not
%   treat fallback steps as "constraint-guaranteed success".
%
% Inputs
%   y      current measured output (p x 1)
%   model  struct with A,B,P,R,y_mean,u_mean (and optionally Sigma_*)
%   opt    struct with H,h (tracked quality rows) and optional u_min/u_max
%   uk     applied fallback input (nu x 1), already saturated if desired
%
% Outputs (info fields)
%   cc_cert_level            always 'uncertified_fallback'
%   exitflag                 -1 (not a successful QP solve)
%   max_cc_violation         NaN  (no A_ch*U-b_ch residual available)
%   mean_pred                one-step mean y prediction (p x 1), or []
%   det_mean_violation       max(H*mean_pred - h); NaN if prediction fails
%   det_mean_safe            true iff det_mean_violation <= 0 (when known)
%   uncertified_fallback     true
%   note                     short human-readable status

info = struct();
info.cc_cert_level = 'uncertified_fallback';
info.exitflag = -1;
info.max_cc_violation = NaN;   % deliberately not a QP residual
info.uncertified_fallback = true;
info.mean_pred = [];
info.det_mean_violation = NaN;
info.det_mean_safe = false;
info.note = 'uncertified_fallback: no chance-constraint guarantee';

uk = uk(:);
try
    z = model.R' * (y(:) - model.y_mean(:));
    z1 = model.A * z + model.B * (uk - model.u_mean(:));
    mean_pred = model.y_mean(:) + model.P * z1;
    info.mean_pred = mean_pred;

    if isfield(opt, 'H') && isfield(opt, 'h') && ~isempty(opt.H)
        residual = opt.H * mean_pred - opt.h(:);
        info.det_mean_violation = max(residual);
        info.det_mean_safe = (info.det_mean_violation <= 0);
        if info.det_mean_safe
            info.note = ['uncertified_fallback: one-step mean under hard y bound, ' ...
                'but NOT chance-constraint certified'];
        else
            info.note = sprintf([ ...
                'uncertified_fallback: one-step mean violates hard y bound ' ...
                '(det_mean_viol=%.3e); NOT chance-constraint certified'], ...
                info.det_mean_violation);
        end
    else
        info.note = 'uncertified_fallback: no H/h for mean check; NOT chance-constraint certified';
    end
catch ME
    info.note = sprintf( ...
        'uncertified_fallback: mean prediction failed (%s); NOT chance-constraint certified', ...
        ME.message);
end
end
