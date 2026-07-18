function [A_ch, b_ch, meta] = build_chance_rows(z, model, opt)
% BUILD_CHANCE_ROWS Construct Boole-tightened chance rows at given alpha_joint.
% Uses the same law as centered_smpc_step (drop-cross unless opt.use_cross_cov).

N = opt.N;
nq = size(opt.H,1);
ny = size(model.P,1);
nu = size(model.B,2);
U0 = repmat(model.u_mean, N, 1);

M = cell(1,N); G = cell(1,N);
for j = 1:N
    M{j} = model.P * (model.A^j);
    G{j} = zeros(ny, N*nu);
    for i = 0:j-1
        G{j}(:, i*nu+1:(i+1)*nu) = model.P*(model.A^(j-1-i))*model.B;
    end
end

risk_each = opt.alpha_joint/(2*nq*N);
z_quantile = norminv(1-risk_each);
Sigma_z = zeros(size(model.A));
A_ch = []; b_ch = [];

use_cross = isfield(opt,'use_cross_cov') && logical(opt.use_cross_cov) ...
    && isfield(model,'Sigma_zo') && ~isempty(model.Sigma_zo);
Szo = [];
if use_cross
    Szo = model.Sigma_zo;
    if size(Szo,1) ~= size(model.A,1) || size(Szo,2) ~= ny
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
        Sigma_y = Sigma_y + model.P*Szo + Szo'*model.P';
        Sigma_y = (Sigma_y + Sigma_y')/2;
    end
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

meta = struct();
meta.risk_each = risk_each;
meta.z_quantile = z_quantile;
meta.use_cross = use_cross;
meta.U0 = U0;
meta.N = N;
meta.nu = nu;
end
