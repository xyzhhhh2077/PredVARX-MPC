function [z, y_pred, U, out] = centered_smpc_step(y, r, model, opt)
% CENTERED_SMPC_STEP One-step centered linear-Gaussian SMPC QP.
% Model coordinates are z=R'*(y-y_mean), v=u-u_mean, y=y_mean+P*z.

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

% Optimize absolute U, but predict and penalize centered input DU=U-U0.
% Expand the raw objective J=sum||e0_j+G_j U||_Q^2+||U-U0||_RuBar^2
% into quadprog form 0.5*U'*Hqp*U+fqp'*U.
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
Hqp = 2*((Hraw+Hraw')/2) + 1e-9*eye(N*nu);
fqp = 2*fraw;

risk_each = opt.alpha_joint/(2*nq*N);
z_quantile = norminv(1-risk_each);
Sigma_z = zeros(size(model.A));
A_ch=[]; b_ch=[];
for j = 1:N
    Sigma_z = model.A*Sigma_z*model.A' + model.Sigma_eps;
    % Tracked-only covariance path: chance rows act on tracked axes, and
    % their projection-out residual is exactly zero by the AL geometry.
    % Therefore do not add the empirical residual proxy to these rows.
    Sigma_y_model = model.P*Sigma_z*model.P';
    Sigma_y = Sigma_y_model;
    mu0 = model.y_mean + M{j}*z - G{j}*U0;
    for q = 1:nq
        hq = opt.H(q,:)';
        tight = z_quantile*sqrt(max(hq'*Sigma_y_model*hq,1e-12));
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
% Remove only the tiny numerical Hessian regularizer from the reported raw cost.
out.cost = 0.5*U'*(Hqp-1e-9*eye(N*nu))*U + fqp'*U + J_const;
out.estimated_sigma_eps = sqrt(trace(model.Sigma_eps)/size(model.Sigma_eps,1));
out.estimated_sigma_obs = sqrt(trace(model.Sigma_obs)/size(model.Sigma_obs,1));
end
