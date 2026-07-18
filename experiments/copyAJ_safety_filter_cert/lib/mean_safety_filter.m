function info = mean_safety_filter(y, r, model, opt, u_prev)
% MEAN_SAFETY_FILTER Opinion-8 L2' candidate after primary chance-QP fail.
%
% Object change (NOT another original-alpha online heuristic):
%   Seek first-step uk in input bounds such that the ONE-STEP DETERMINISTIC
%   mean prediction of tracked outputs satisfies hard bounds
%       |h_q' * y_mean_pred(uk)|  <=  h_q - margin
%   (equivalently H*y_mean_pred <= h - margin, componentwise).
%
% Ladder (first success wins):
%   1) one-step QP: track r + Ru, s.t. mean hard bounds + input bounds
%   2) check u_prev against mean hard bounds (if provided)
%   3) check u_mean against mean hard bounds
%
% Certificate level on success:
%   info.cc_cert_level = 'mean_safety_filter'
% This is EXPLICITLY NOT an original-alpha chance certificate, NOT
% recursive feasibility, and NOT a multi-step joint chance recovery.
%
% On total failure: success=false, level='uncertified_fallback'.

if nargin < 5, u_prev = []; end

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
info.max_mean_violation = NaN;
info.cost = NaN;
info.note = 'mean_safety_filter failed';
info.is_original_alpha = false;
info.cert_object = 'one_step_deterministic_mean_hard_bound';

nu = size(model.B, 2);
z = model.R' * (y - model.y_mean);
G1 = model.P * model.B;                 % dy/du at one step (centered v)
mu_free = model.y_mean + model.P*(model.A*z) - G1*model.u_mean;  % y_mean at u=0 abs? see below
% y_pred_mean(u) = y_mean + P*A*z + P*B*(u - u_mean) = mu0 + G1*u
% with mu0 = y_mean + P*A*z - G1*u_mean
mu0 = mu_free;

if isfield(opt,'u_min'), u_min = opt.u_min; else, u_min = -inf; end
if isfield(opt,'u_max'), u_max = opt.u_max; else, u_max = inf; end
if isscalar(u_min), lb = u_min*ones(nu,1); else, lb = u_min(:); end
if isscalar(u_max), ub = u_max*ones(nu,1); else, ub = u_max(:); end

H = opt.H; h = opt.h(:);
nq = size(H,1);
% Optional mean margin (shrink hard bound slightly for numerical room)
if isfield(opt,'mean_bound_margin') && ~isempty(opt.mean_bound_margin)
    margin = opt.mean_bound_margin;
else
    margin = 0;
end
h_eff = h - margin;

% Mean hard-bound rows: H * (mu0 + G1*u) <= h_eff
%   (H*G1) u  <=  h_eff - H*mu0
A_mean = H * G1;
b_mean = h_eff - H * mu0;

% Also lower-side if constraints are one-sided upper on tracked axes only.
% Primary SMPC uses both +H and -H chance rows (box). Mirror for mean filter.
A_mean = [A_mean; -H * G1];
b_mean = [b_mean; h_eff + H * mu0];

tol = 1e-7;
qpopt = optimset('Display','off');

% --- attempt 1: one-step mean-constrained tracking QP ---
try
    e0 = mu0 - r;  % when u=0 contribution already in mu0; cost on uk absolute via G1*(uk)
    % minimize ||mu0 + G1*uk - r||_Q^2 + ||uk - u_mean||_Ru^2
    Hraw = G1'*opt.Q*G1 + opt.Ru;
    fraw = G1'*opt.Q*e0 - opt.Ru*model.u_mean;
    Hqp = 2*((Hraw+Hraw')/2) + 1e-9*eye(nu);
    fqp = 2*fraw;
    [uk, ~, ef] = quadprog(Hqp, fqp, A_mean, b_mean, [], [], lb, ub, [], qpopt);
    if ef > 0
        ypred = mu0 + G1*uk;
        mv = max([H*ypred - h_eff; -H*ypred - h_eff]);
        if mv <= tol
            info.success = true;
            info.cc_cert_level = 'mean_safety_filter';
            info.mode = 'one_step_mean_qp';
            info.stage_class = 'one_step_mean_qp';
            info.stage_param = margin;
            info.uk = uk;
            info.U = uk;
            info.y_pred = ypred;
            info.exitflag = ef;
            info.max_mean_violation = mv;
            info.cost = 0.5*uk'*(Hqp-1e-9*eye(nu))*uk + fqp'*uk;
            info.note = ['L2_prime mean_safety_filter: one-step deterministic mean ' ...
                'hard bound (NOT original-alpha chance cert; NOT RF)'];
            return;
        end
    end
catch
end

% --- attempt 2: backup candidates screened by mean hard bound ---
cands = {};
clabels = {};
if ~isempty(u_prev)
    up = u_prev(:);
    if numel(up) == nu
        cands{end+1} = min(max(up, lb), ub); %#ok<AGROW>
        clabels{end+1} = 'backup_u_prev_mean'; %#ok<AGROW>
    end
end
um = min(max(model.u_mean(:), lb), ub);
cands{end+1} = um;
clabels{end+1} = 'backup_u_mean_mean';

for ic = 1:numel(cands)
    uk = cands{ic};
    ypred = mu0 + G1*uk;
    resid = [H*ypred - h_eff; -H*ypred - h_eff];
    mv = max(resid);
    bound_ok = all(uk >= lb-1e-9) && all(uk <= ub+1e-9);
    if bound_ok && mv <= tol
        info.success = true;
        info.cc_cert_level = 'mean_safety_filter';
        info.mode = clabels{ic};
        info.stage_class = clabels{ic};
        info.stage_param = margin;
        info.uk = uk;
        info.U = uk;
        info.y_pred = ypred;
        info.exitflag = 1;
        info.max_mean_violation = mv;
        info.cost = NaN;
        info.note = ['L2_prime mean_safety_filter: ' clabels{ic} ...
            ' passed one-step mean hard bound (NOT original-alpha)'];
        return;
    end
end

% failure: suggest sat(u_mean) to caller but do not certify
info.uk = um;
info.U = um;
info.y_pred = mu0 + G1*um;
info.max_mean_violation = max([H*info.y_pred - h_eff; -H*info.y_pred - h_eff]);
info.note = 'mean_safety_filter: no candidate satisfied one-step mean hard bounds';
end
