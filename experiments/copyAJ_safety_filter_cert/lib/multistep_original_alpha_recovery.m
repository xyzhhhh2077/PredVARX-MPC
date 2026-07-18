function info = multistep_original_alpha_recovery(y, r, model, opt, u_prev)
% MULTISTEP_ORIGINAL_ALPHA_RECOVERY Opinion-8 stronger backup.
% Accept only FULL multi-step plans U that satisfy ORIGINAL alpha_joint rows:
%   max(A_ch(alpha)*U - b_ch(alpha)) <= tol
%
% Ladder (first success wins):
%   1) shorter horizon QP with SAME alpha; accept full U if original rows pass
%   2) one-step original-alpha QP; accept if rows pass
%   3) pure feasibility QP at original alpha (min ||U-U0||^2 s.t. chance+bounds)
%   4) reduced tracking weight Q/10 at original alpha
%   5) certified constant-hold backups (u_prev / u_mean) as last resort
% Else uncertified sat(u_mean).
%
% Non-claims: not recursive feasibility; not infinite-horizon joint chance;
% not closed-loop stability.

if nargin < 5, u_prev = []; end
tol = 1e-7;
nu = size(model.B,2);
if isscalar(opt.u_min), lb1 = opt.u_min*ones(nu,1); else, lb1 = opt.u_min(:); end
if isscalar(opt.u_max), ub1 = opt.u_max*ones(nu,1); else, ub1 = opt.u_max(:); end

info = struct();
info.success = false;
info.cc_cert_level = 'uncertified_fallback';
info.mode = 'none';
info.stage_class = 'none';
info.uk = [];
info.U = [];
info.exitflag = -1;
info.max_cc_violation = NaN;
info.cost = NaN;
info.note = 'multistep original-alpha recovery failed';
info.cert_tol = tol;

% --- helpers as nested local functions via subfunctions at end of file style ---
% 1) short horizon same alpha
N0 = opt.N;
Ntries = unique([max(2,floor(N0/2)), max(2,floor(N0/3)), max(2,floor(N0/4)), 2]);
for Ntry = Ntries
    if Ntry >= N0, continue; end
    opt_try = opt; opt_try.N = Ntry;
    try
        [~, ypred, U, out] = centered_smpc_step(y, r, model, opt_try);
        if out.exitflag > 0
            % Re-evaluate this U under ORIGINAL N chance rows by padding/trunc?
            % Safer: rebuild original-N rows and embed U as first Ntry steps + hold last
            U_full = local_embed_plan(U, nu, N0, Ntry);
            cert = local_cert_full_U(U_full, y, model, opt, tol);
            if cert.pass
                info = local_accept(info, U_full(1:nu), U_full, out, cert, ypred, ...
                    sprintf('short_horizon_N%d_fullU', Ntry), 'short_horizon_fullU');
                return;
            end
            % Also accept if the short plan's own rows are original-alpha with N=Ntry
            % AND we apply only first step while re-certifying first step under full N hold.
            % Already covered by full embed. Keep short-plan residual as diagnostic.
        end
    catch
    end
end

% --- 2) one-step original alpha ---
try
    opt1 = opt; opt1.N = 1;
    [~, ypred, U, out] = centered_smpc_step(y, r, model, opt1);
    if out.exitflag > 0
        U_full = local_embed_plan(U, nu, N0, 1);
        cert = local_cert_full_U(U_full, y, model, opt, tol);
        if cert.pass
            info = local_accept(info, U(1:nu), U_full, out, cert, ypred, ...
                'one_step_embed_fullU', 'one_step_fullU');
            return;
        end
        % accept pure one-step certificate at N=1 under original alpha, then embed hold
        mv = max(out.A_ch*U - out.b_ch);
        if mv <= tol
            % still need full-N original rows for applied policy certificate
            cert = local_cert_full_U(U_full, y, model, opt, tol);
            if cert.pass
                info = local_accept(info, U(1:nu), U_full, out, cert, ypred, ...
                    'one_step_original_alpha', 'one_step_fullU');
                return;
            end
        end
    end
catch
end

% --- 3) feasibility QP at original alpha: min ||U-U0||^2 s.t. rows ---
try
    [U_feas, feas_ok, feas_mv, ypred] = local_feasibility_qp(y, model, opt, tol);
    if feas_ok
        fake_out = struct('exitflag',1,'cost',norm(U_feas)^2);
        cert = struct('pass',true,'max_cc_violation',feas_mv,'tol',tol);
        info = local_accept(info, U_feas(1:nu), U_feas, fake_out, cert, ypred, ...
            'feasibility_qp_original_alpha', 'feasibility_fullU');
        return;
    end
catch
end

% --- 4) reduced Q tracking at original alpha ---
try
    optQ = opt; optQ.Q = opt.Q / 10;
    [~, ypred, U, out] = centered_smpc_step(y, r, model, optQ);
    if out.exitflag > 0
        mv = max(out.A_ch*U - out.b_ch);
        if mv <= tol
            cert = struct('pass',true,'max_cc_violation',mv,'tol',tol);
            info = local_accept(info, U(1:nu), U, out, cert, ypred, ...
                'reduced_Q_original_alpha', 'reduced_Q_fullU');
            return;
        end
    end
catch
end

% --- 5) constant-hold certified backups ---
cands = {};
if ~isempty(u_prev), cands{end+1} = {'u_prev', min(max(u_prev(:),lb1),ub1)}; end %#ok<*AGROW>
cands{end+1} = {'u_mean', min(max(model.u_mean(:),lb1),ub1)};
for ic = 1:numel(cands)
    tag = cands{ic}{1}; uk = cands{ic}{2};
    cert = certify_input_original_alpha(uk, y, model, opt, tol);
    if cert.pass
        U_full = repmat(uk, opt.N, 1);
        z = model.R'*(y-model.y_mean);
        ypred = model.y_mean + model.P*(model.A*z) + model.P*model.B*(uk-model.u_mean);
        fake_out = struct('exitflag',1,'cost',0);
        info = local_accept(info, uk, U_full, fake_out, cert, ypred, ...
            ['hold_' tag '_certified'], ['backup_' tag]);
        return;
    end
end

% fail
uk = min(max(model.u_mean(:),lb1),ub1);
info.success = false;
info.cc_cert_level = 'uncertified_fallback';
info.mode = 'uncertified_u_mean';
info.stage_class = 'uncertified';
info.uk = uk;
info.U = repmat(uk, opt.N, 1);
info.exitflag = -1;
cert = certify_input_original_alpha(uk, y, model, opt, tol);
info.max_cc_violation = cert.max_cc_violation;
info.note = 'no multistep original-alpha plan found; uncertified sat(u_mean)';
z = model.R'*(y-model.y_mean);
info.y_pred = model.y_mean + model.P*(model.A*z);
end

function U_full = local_embed_plan(U_short, nu, Nfull, Nshort)
U_short = U_short(:);
nush = numel(U_short);
assert(mod(nush,nu)==0);
Ns = nush/nu;
assert(Ns==Nshort);
U_full = zeros(Nfull*nu,1);
U_full(1:nush) = U_short;
uk_last = U_short(end-nu+1:end);
for j = Nshort+1:Nfull
    U_full((j-1)*nu+1:j*nu) = uk_last;
end
end

function cert = local_cert_full_U(U, y, model, opt, tol)
z = model.R' * (y - model.y_mean);
[A_ch, b_ch] = build_chance_rows(z, model, opt);
U = U(:);
assert(numel(U)==size(A_ch,2), 'U dim mismatch for chance rows');
% bounds
nu = size(model.B,2); N = opt.N;
if isscalar(opt.u_min), lb = repmat(opt.u_min*ones(nu,1),N,1); else, lb = repmat(opt.u_min(:),N,1); end
if isscalar(opt.u_max), ub = repmat(opt.u_max*ones(nu,1),N,1); else, ub = repmat(opt.u_max(:),N,1); end
bound_ok = all(U>=lb-1e-9) && all(U<=ub+1e-9);
mv = max(A_ch*U - b_ch);
cert = struct();
cert.pass = bound_ok && (mv <= tol);
cert.max_cc_violation = mv;
cert.bound_ok = bound_ok;
cert.tol = tol;
end

function info = local_accept(info, uk, U, out, cert, ypred, mode, stage)
info.success = true;
info.cc_cert_level = 'original_alpha_certified';
info.mode = mode;
info.stage_class = stage;
info.uk = uk;
info.U = U;
info.exitflag = out.exitflag;
info.max_cc_violation = cert.max_cc_violation;
if isfield(out,'cost'), info.cost = out.cost; else, info.cost = NaN; end
info.y_pred = ypred;
info.note = ['accepted multistep original-alpha plan: ' mode];
info.cert = cert;
end

function [U, ok, mv, ypred] = local_feasibility_qp(y, model, opt, tol)
% min 0.5||U-U0||^2 s.t. original chance rows + bounds
z = model.R'*(y-model.y_mean);
[A_ch, b_ch, meta] = build_chance_rows(z, model, opt);
N = meta.N; nu = meta.nu; U0 = meta.U0;
Hqp = eye(N*nu);
fqp = -U0; % 0.5||U-U0||^2 => H=I, f=-U0 in 0.5 U'HU + f'U form with H=I -> use quadprog 0.5U'U - U0'U
Hqp = eye(N*nu);
fqp = -U0;
if isscalar(opt.u_min), lb = repmat(opt.u_min*ones(nu,1),N,1); else, lb = repmat(opt.u_min(:),N,1); end
if isscalar(opt.u_max), ub = repmat(opt.u_max*ones(nu,1),N,1); else, ub = repmat(opt.u_max(:),N,1); end
opts = optimset('Display','off');
[U,~,ef] = quadprog(Hqp, fqp, A_ch, b_ch, [], [], lb, ub, [], opts);
ok = false; mv = NaN; ypred = [];
if ef > 0
    mv = max(A_ch*U - b_ch);
    ok = mv <= tol;
    if ok
        G1 = model.P*model.B;
        ypred = model.y_mean + model.P*(model.A*z) + G1*(U(1:nu)-model.u_mean);
    end
end
end
