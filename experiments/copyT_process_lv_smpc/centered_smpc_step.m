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

Hqp = kron(eye(N), opt.Ru);
fqp = zeros(N*nu,1);
for j = 1:N
    ej = M{j}*z + model.y_mean - r;
    Hqp = Hqp + G{j}'*opt.Q*G{j};
    fqp = fqp + G{j}'*opt.Q*ej;
end
Hqp = (Hqp+Hqp')/2 + 1e-9*eye(N*nu);
J_const = 0;
for j = 1:N
    ej = M{j}*z + model.y_mean - r;
    J_const = J_const + ej'*opt.Q*ej;
end

risk_each = opt.alpha_joint/(2*nq*N);
z_quantile = norminv(1-risk_each);
Sigma_z = zeros(size(model.A));
A_ch=[]; b_ch=[];
for j = 1:N
    Sigma_z = model.A*Sigma_z*model.A' + model.Sigma_eps;
    Sigma_y = model.P*Sigma_z*model.P' + model.Sigma_obs;
    mu0 = model.y_mean + M{j}*z - G{j}*U0;
    for q = 1:nq
        hq = opt.H(q,:)';
        tight = z_quantile*sqrt(max(hq'*Sigma_y*hq,1e-12));
        A_ch = [A_ch; hq'*G{j}; -hq'*G{j}];
        b_ch = [b_ch; opt.h(q)-hq'*mu0-tight; opt.h(q)+hq'*mu0-tight];
    end
end

lb = repmat(opt.u_min,N,1); ub = repmat(opt.u_max,N,1);
qpopt = optimset('Display','off');
[U,~,exitflag] = quadprog(Hqp,fqp,A_ch,b_ch,[],[],lb,ub,[],qpopt);
if exitflag <= 0
    error('centered_smpc_step:Infeasible','Chance-constrained QP is infeasible.');
end
y_pred = model.y_mean + M{1}*z + G{1}*(U-U0);
out.A_ch=A_ch; out.b_ch=b_ch; out.risk_each=risk_each;
out.z_quantile=z_quantile; out.exitflag=exitflag;
out.cost = 0.5*U'*Hqp*U + fqp'*U + J_const;
out.estimated_sigma_eps = sqrt(trace(model.Sigma_eps)/size(model.Sigma_eps,1));
out.estimated_sigma_obs = sqrt(trace(model.Sigma_obs)/size(model.Sigma_obs,1));
end
