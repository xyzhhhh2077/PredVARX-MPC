function info = original_alpha_recovery(y, r, model, opt, u_prev)
% ORIGINAL_ALPHA_RECOVERY Opinion-8 guarantee-layer recovery.
% Only accepts candidates that re-pass ORIGINAL alpha_joint chance rows
% under constant-hold extension of the first-step input.
%
% Ladder (first success wins):
%   1) shorter horizon QP with SAME alpha_joint, then re-certify uk
%   2) one-step original-alpha chance QP
%   3) previous input u_prev if certified
%   4) u_mean if certified
% Else: uncertified sat(u_mean)
%
% Explicit non-claims:
% - Not recursive feasibility
% - Not infinite-horizon joint chance certificate
% - Short-horizon solve is only a candidate generator; acceptance needs
%   original-alpha re-check on the applied uk.

if nargin < 5, u_prev = []; end

info = struct();
info.success = false;
info.cc_cert_level = 'uncertified_fallback';
info.mode = 'none';
info.stage_class = 'none';
info.uk = [];
info.exitflag = -1;
info.max_cc_violation = NaN;
info.cost = NaN;
info.note = 'original-alpha recovery failed';
info.cert = struct();

nu = size(model.B,2);
if isscalar(opt.u_min), lb = opt.u_min*ones(nu,1); else, lb = opt.u_min(:); end
if isscalar(opt.u_max), ub = opt.u_max*ones(nu,1); else, ub = opt.u_max(:); end
u_mean = model.u_mean(:);

% --- 1) short horizon, SAME alpha ---
N0 = opt.N;
Ntries = unique([max(2,floor(N0/2)), max(2,floor(N0/3)), 2]);
for Ntry = Ntries
    if Ntry >= N0, continue; end
    opt_try = opt;
    opt_try.N = Ntry;
    try
        [~, ypred, U, out] = centered_smpc_step(y, r, model, opt_try);
        if out.exitflag > 0
            uk = U(1:nu);
            cert = certify_input_original_alpha(uk, y, model, opt);
            if cert.pass
                info.success = true;
                info.cc_cert_level = 'original_alpha_certified';
                info.mode = sprintf('short_horizon_N%d_recheck', Ntry);
                info.stage_class = 'short_horizon_original_alpha';
                info.uk = uk;
                info.y_pred = ypred;
                info.exitflag = out.exitflag;
                info.max_cc_violation = cert.max_cc_violation;
                info.cost = out.cost;
                info.cert = cert;
                info.note = ['accepted after original-alpha recheck: ' info.mode];
                return;
            end
        end
    catch
    end
end

% --- 2) one-step original-alpha chance QP ---
try
    opt1 = opt;
    opt1.N = 1;
    [~, ypred, U, out] = centered_smpc_step(y, r, model, opt1);
    if out.exitflag > 0
        uk = U(1:nu);
        cert = certify_input_original_alpha(uk, y, model, opt);
        if cert.pass
            info.success = true;
            info.cc_cert_level = 'original_alpha_certified';
            info.mode = 'one_step_original_alpha';
            info.stage_class = 'one_step_original_alpha';
            info.uk = uk;
            info.y_pred = ypred;
            info.exitflag = out.exitflag;
            info.max_cc_violation = cert.max_cc_violation;
            info.cost = out.cost;
            info.cert = cert;
            info.note = 'accepted one-step original-alpha QP after recheck';
            return;
        end
    end
catch
end

% --- 3) previous input ---
if ~isempty(u_prev)
    uk = min(max(u_prev(:), lb), ub);
    cert = certify_input_original_alpha(uk, y, model, opt);
    if cert.pass
        info.success = true;
        info.cc_cert_level = 'original_alpha_certified';
        info.mode = 'hold_u_prev_certified';
        info.stage_class = 'backup_prev';
        info.uk = uk;
        info.exitflag = 1;
        info.max_cc_violation = cert.max_cc_violation;
        info.cert = cert;
        info.note = 'previous input passes original-alpha constant-hold certificate';
        z = model.R'*(y-model.y_mean);
        info.y_pred = model.y_mean + model.P*(model.A*z) + model.P*model.B*(uk-u_mean);
        return;
    end
end

% --- 4) u_mean backup ---
uk = min(max(u_mean, lb), ub);
cert = certify_input_original_alpha(uk, y, model, opt);
if cert.pass
    info.success = true;
    info.cc_cert_level = 'original_alpha_certified';
    info.mode = 'u_mean_certified';
    info.stage_class = 'backup_umean';
    info.uk = uk;
    info.exitflag = 1;
    info.max_cc_violation = cert.max_cc_violation;
    info.cert = cert;
    info.note = 'u_mean passes original-alpha constant-hold certificate';
    z = model.R'*(y-model.y_mean);
    info.y_pred = model.y_mean + model.P*(model.A*z) + model.P*model.B*(uk-u_mean);
    return;
end

% --- fail: uncertified hold ---
info.success = false;
info.cc_cert_level = 'uncertified_fallback';
info.mode = 'uncertified_u_mean';
info.stage_class = 'uncertified';
info.uk = uk;
info.exitflag = -1;
info.max_cc_violation = cert.max_cc_violation;
info.cert = cert;
info.note = 'no original-alpha-certified recovery; applied sat(u_mean) uncertified';
z = model.R'*(y-model.y_mean);
info.y_pred = model.y_mean + model.P*(model.A*z);
end
