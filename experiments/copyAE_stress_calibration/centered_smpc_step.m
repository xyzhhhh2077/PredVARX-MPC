function [z, y_pred, U, out] = centered_smpc_step(y, r, model, opt)
% CENTERED_SMPC_STEP One-step centered linear-Gaussian SMPC QP.
% Model coordinates are z=R'*(y-y_mean), v=u-u_mean, y=y_mean+P*z.
%
% Opinion-trial options (all default OFF / baseline):
%   opt.use_cross_cov (false): if true and model.Sigma_zo provided, use
%       Sigma_y = P Sigma_z P' + P Sigma_zo + Sigma_zo' P' + Sigma_obs
%   opt.use_terminal_cost (false): optional terminal quadratic on z_N

nq = size(opt.H,1);
N = opt.N;
ny = size(model.P,1);
nu = size(model.B,2);
z = model.R' * (y - model.y_mean);
U0 = repmat(model.u_mean, N, 1);

M = cell(1,N); G = cell(1,N);
for j = 1:N
    M{j} = model.P * (model.A^j);
    G{j} = zeros(ny, N*nu);
    for i = 0:j-1
        G{j}(:,i*nu+1:(i+1)*nu) = model.P*(model.A^(j-1-i))*model.B;
    end
end

RuBar = kron(eye(N), opt.Ru);
Hraw = RuBar;
fraw = -RuBar*U0;
J_const = U0'*RuBar*U0;
for j = 1:N
    ej0 = M{j}*z + model.y_mean - G{j}*U0 - r;
    Hraw = Hraw + G{j}'*opt.Q*G{j};
    fraw = fraw + G{j}'*opt.Q*ej0;
    J_const = J_const + ej0'*opt.Q*ej0;
end

use_terminal_cost = isfield(opt,'use_terminal_cost') && logical(opt.use_terminal_cost);
terminal_cost_applied = false;
Pterm = [];
if use_terminal_cost
    try
        nz = size(model.A,1);
        if isfield(opt,'Qf') && ~isempty(opt.Qf)
            Qf = (opt.Qf + opt.Qf')/2;
        else
            Qf = model.P' * opt.Q * model.P;
            Qf = (Qf + Qf')/2;
            if min(eig(Qf + 1e-12*eye(nz))) <= 0
                Qf = eye(nz);
            end
        end
        Pterm = dlyap(model.A', Qf);
        Pterm = (Pterm + Pterm')/2;
        if any(~isfinite(Pterm(:))) || min(eig(Pterm + 1e-12*eye(nz))) <= 0
            error('centered_smpc_step:BadPterm','Pterm not finite PSD.');
        end
        Gz = zeros(nz, N*nu);
        for i = 0:N-1
            Gz(:, i*nu+1:(i+1)*nu) = (model.A^(N-1-i)) * model.B;
        end
        zN0 = (model.A^N)*z - Gz*U0;
        Hraw = Hraw + Gz' * Pterm * Gz;
        fraw = fraw + Gz' * Pterm * zN0;
        J_const = J_const + zN0' * Pterm * zN0;
        terminal_cost_applied = true;
    catch
        terminal_cost_applied = false;
        Pterm = [];
    end
end

Hqp = 2*((Hraw+Hraw')/2) + 1e-9*eye(N*nu);
fqp = 2*fraw;

risk_each = opt.alpha_joint/(2*nq*N);
z_quantile = norminv(1-risk_each);
Sigma_z = zeros(size(model.A));
A_ch=[]; b_ch=[];
use_cross = isfield(opt,'use_cross_cov') && logical(opt.use_cross_cov) ...
    && isfield(model,'Sigma_zo') && ~isempty(model.Sigma_zo);
Szo = [];
if use_cross
    Szo = model.Sigma_zo;
    if size(Szo,1) ~= size(model.A,1) || size(Szo,2) ~= ny
        % expect ell x p; if p x ell, transpose
        if size(Szo,1) == ny && size(Szo,2) == size(model.A,1)
            Szo = Szo';
        else
            use_cross = false;
        end
    end
end

for j = 1:N
    Sigma_z = model.A*Sigma_z*model.A' + model.Sigma_eps;
    Sigma_y = model.P*Sigma_z*model.P' + model.Sigma_obs;
    if use_cross
        % Full block: P Sigma_z P' + P Sigma_zo + Sigma_zo' P' + Sigma_obs
        % with Sigma_zo = Cov(z, o) and o residual part absorbed in Sigma_obs
        % when Sigma_obs stands for Sigma_o; here we add cross only.
        Sigma_y = Sigma_y + model.P*Szo + Szo'*model.P';
        Sigma_y = (Sigma_y + Sigma_y')/2;
    end
    % Numerical PSD floor for quantile sqrt
    [Vsy, Dsy] = eig((Sigma_y+Sigma_y')/2);
    dsy = max(real(diag(Dsy)), 0);
    Sigma_y = Vsy*diag(dsy)*Vsy';
    mu0 = model.y_mean + M{j}*z - G{j}*U0;
    for q = 1:nq
        hq = opt.H(q,:)';
        tight = z_quantile*sqrt(max(hq'*Sigma_y*hq,1e-12));
        A_ch = [A_ch; hq'*G{j}; -hq'*G{j}];
        b_ch = [b_ch; opt.h(q)-hq'*mu0-tight; opt.h(q)+hq'*mu0-tight];
    end
end

if isscalar(opt.u_min)
    lb_step = opt.u_min*ones(nu,1);
else
    lb_step = opt.u_min(:);
end
if isscalar(opt.u_max)
    ub_step = opt.u_max*ones(nu,1);
else
    ub_step = opt.u_max(:);
end
assert(numel(lb_step)==nu && numel(ub_step)==nu, ...
    'Input bounds must be scalar or have one entry per input channel.');
lb = repmat(lb_step,N,1); ub = repmat(ub_step,N,1);
qpopt = optimset('Display','off');
[U,~,exitflag] = quadprog(Hqp,fqp,A_ch,b_ch,[],[],lb,ub,[],qpopt);
if exitflag <= 0
    error('centered_smpc_step:Infeasible','Chance-constrained QP is infeasible.');
end
y_pred = model.y_mean + M{1}*z + G{1}*(U-U0);
out.A_ch=A_ch; out.b_ch=b_ch; out.risk_each=risk_each;
out.z_quantile=z_quantile; out.exitflag=exitflag;
out.lb=lb; out.ub=ub; out.U0=U0;
out.use_terminal_cost = use_terminal_cost;
out.terminal_cost_applied = terminal_cost_applied;
out.Pterm = Pterm;
out.use_cross_cov = use_cross;
out.cost = 0.5*U'*(Hqp-1e-9*eye(N*nu))*U + fqp'*U + J_const;
out.estimated_sigma_eps = sqrt(trace(model.Sigma_eps)/size(model.Sigma_eps,1));
out.estimated_sigma_obs = sqrt(trace(model.Sigma_obs)/size(model.Sigma_obs,1));
end
