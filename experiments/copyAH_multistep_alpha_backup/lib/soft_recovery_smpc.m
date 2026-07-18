function info = soft_recovery_smpc(y, r, model, opt)
% SOFT_RECOVERY_SMPC Opinion 8 L2 candidate: recover a QP-feasible input
% after primary chance-constrained QP fails.
%
% Recovery ladder (first success wins):
%   1) inflate joint risk (looser Boole allocation)
%   2) shorten horizon N
%   3) drop chance rows, keep input bounds only
%
% On success: info.cc_cert_level = 'soft_recovery' and info.mode names the step.
% On total failure: info.cc_cert_level = 'uncertified_fallback' (caller may still
% apply sat(u_mean)). Soft recovery is NOT the original risk certificate.
%
% Default non-claim: recovered alpha'/N/no-CC is a different certificate.

info = struct();
info.success = false;
info.cc_cert_level = 'uncertified_fallback';
info.mode = 'none';
info.stage_class = 'none';
info.stage_param = NaN;
info.uk = [];
info.U = [];
info.y_pred = [];
info.exitflag = -1;
info.max_cc_violation = NaN;
info.cost = NaN;
info.note = 'soft recovery failed';

nu = size(model.B, 2);
if isfield(opt, 'u_min'), u_min = opt.u_min; else, u_min = -inf; end
if isfield(opt, 'u_max'), u_max = opt.u_max; else, u_max = inf; end

% --- attempt 1: risk inflation ---
risk_list = [0.20, 0.35, 0.50];
for ii = 1:numel(risk_list)
    opt_try = opt;
    opt_try.alpha_joint = risk_list(ii);
    try
        [~, ypred, U, out] = centered_smpc_step(y, r, model, opt_try);
        if out.exitflag > 0
            info.success = true;
            info.cc_cert_level = 'soft_recovery';
            info.mode = sprintf('risk_inflate_%.2f', risk_list(ii));
            info.stage_class = 'risk_inflate';
            info.stage_param = risk_list(ii);
            info.uk = U(1:nu);
            info.U = U;
            info.y_pred = ypred;
            info.exitflag = out.exitflag;
            info.max_cc_violation = max(out.A_ch*U - out.b_ch);
            info.cost = out.cost;
            info.note = ['soft_recovery: ' info.mode ...
                ' (NOT original alpha_joint certificate)'];
            return;
        end
    catch
    end
end

% --- attempt 2: shorter horizon ---
N0 = opt.N;
for Ntry = unique([max(2, floor(N0/2)), max(2, floor(N0/3)), 2])
    if Ntry >= N0, continue; end
    opt_try = opt;
    opt_try.N = Ntry;
    try
        [~, ypred, U, out] = centered_smpc_step(y, r, model, opt_try);
        if out.exitflag > 0
            info.success = true;
            info.cc_cert_level = 'soft_recovery';
            info.mode = sprintf('short_horizon_N%d', Ntry);
            info.stage_class = 'short_horizon';
            info.stage_param = Ntry;
            info.uk = U(1:nu);
            info.U = U;
            info.y_pred = ypred;
            info.exitflag = out.exitflag;
            info.max_cc_violation = max(out.A_ch*U - out.b_ch);
            info.cost = out.cost;
            info.note = ['soft_recovery: ' info.mode ...
                ' (NOT original horizon certificate)'];
            return;
        end
    catch
    end
end

% --- attempt 3: input-bound only QP (no chance rows) ---
try
    opt_try = opt;
    opt_try.alpha_joint = 0.99;  % almost no tightening if used
    % Direct bound-only solve on first-step tracking linearization
    z = model.R' * (y - model.y_mean);
    nu = size(model.B, 2);
    G1 = model.P * model.B;
    e0 = model.y_mean + model.P*(model.A*z) - model.P*model.B*model.u_mean - r;
    % one-step: minimize ||e0 + G1*(u-u_mean)||_Q^2 + ||u-u_mean||_Ru^2
    Hraw = G1'*opt.Q*G1 + opt.Ru;
    fraw = G1'*opt.Q*e0;
    Hqp = 2*((Hraw+Hraw')/2) + 1e-9*eye(nu);
    fqp = 2*fraw;
    if isscalar(u_min), lb = u_min*ones(nu,1); else, lb = u_min(:); end
    if isscalar(u_max), ub = u_max*ones(nu,1); else, ub = u_max(:); end
    [uk,~,ef] = quadprog(Hqp, fqp, [], [], [], [], lb, ub, [], optimset('Display','off'));
    if ef > 0
        info.success = true;
        info.cc_cert_level = 'soft_recovery';
        info.mode = 'input_bound_only';
        info.stage_class = 'bound_only';
        info.stage_param = 0;
        info.uk = uk;
        info.U = uk;
        info.y_pred = model.y_mean + model.P*(model.A*z) + G1*(uk-model.u_mean);
        info.exitflag = ef;
        info.max_cc_violation = NaN;
        info.cost = 0.5*uk'*(Hqp-1e-9*eye(nu))*uk + fqp'*uk;
        info.note = 'soft_recovery: input_bound_only (NO chance-constraint certificate)';
        return;
    end
catch
end
end
