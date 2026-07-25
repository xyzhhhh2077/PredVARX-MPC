function [zCurrent,yPred,U,out] = companionSmpcStep(xi,reference,model,opt)
%COMPANIONSMPCSTEP Chance-constrained MPC using the exact ARX companion state.
% The decision is absolute U, while the companion input is v=U-u_mean.
% Nominal prediction and covariance propagation both use (Ac,Bc,Cc,Qc).

stateDimension = size(model.Ac,1);
nu = size(model.Bc,2);
ny = size(model.Cc,1);
latentDimension = model.companion.latent_dimension;
assert(numel(xi)==stateDimension,'xi does not match the companion dimension.');
assert(isequal(size(model.Ac),[stateDimension stateDimension]));
assert(isequal(size(model.Qc),[stateDimension stateDimension]));
zCurrent = xi(1:latentDimension);
N = opt.N;
constraintCount = size(opt.H,1);
U0 = repmat(model.u_mean,N,1);

M = cell(1,N);
G = cell(1,N);
for horizonIndex = 1:N
    M{horizonIndex} = model.Cc*(model.Ac^horizonIndex);
    G{horizonIndex} = zeros(ny,N*nu);
    for inputIndex = 0:horizonIndex-1
        columns = inputIndex*nu+(1:nu);
        G{horizonIndex}(:,columns) = model.Cc* ...
            (model.Ac^(horizonIndex-1-inputIndex))*model.Bc;
    end
end

RuBar = kron(eye(N),opt.Ru);
Hraw = RuBar;
fraw = -RuBar*U0;
costConstant = U0'*RuBar*U0;
for horizonIndex = 1:N
    baseError = M{horizonIndex}*xi+model.y_mean- ...
        G{horizonIndex}*U0-reference;
    Hraw = Hraw+G{horizonIndex}'*opt.Q*G{horizonIndex};
    fraw = fraw+G{horizonIndex}'*opt.Q*baseError;
    costConstant = costConstant+baseError'*opt.Q*baseError;
end
Hqp = 2*((Hraw+Hraw')/2)+1e-9*eye(N*nu);
fqp = 2*fraw;

riskEach = opt.alpha_joint/(2*constraintCount*N);
zQuantile = norminv(1-riskEach);
SigmaXi = zeros(stateDimension);
Ach = [];
bch = [];
for horizonIndex = 1:N
    SigmaXi = model.Ac*SigmaXi*model.Ac'+model.Qc;
    SigmaY = model.Cc*SigmaXi*model.Cc';
    nominalOffset = model.y_mean+M{horizonIndex}*xi- ...
        G{horizonIndex}*U0;
    for constraintIndex = 1:constraintCount
        direction = opt.H(constraintIndex,:)';
        tightening = zQuantile*sqrt(max(direction'*SigmaY*direction,1e-12));
        Ach = [Ach;direction'*G{horizonIndex}; ...
            -direction'*G{horizonIndex}]; %#ok<AGROW>
        bch = [bch;opt.h(constraintIndex)-direction'*nominalOffset-tightening; ...
            opt.h(constraintIndex)+direction'*nominalOffset-tightening]; %#ok<AGROW>
    end
end

if isscalar(opt.u_min)
    lowerStep = opt.u_min*ones(nu,1);
else
    lowerStep = opt.u_min(:);
end
if isscalar(opt.u_max)
    upperStep = opt.u_max*ones(nu,1);
else
    upperStep = opt.u_max(:);
end
assert(numel(lowerStep)==nu && numel(upperStep)==nu, ...
    'Input bounds must be scalar or have one entry per channel.');
lb = repmat(lowerStep,N,1);
ub = repmat(upperStep,N,1);
qpOptions = optimset('Display','off');
[U,~,exitflag] = callMathWorksQuadprog( ...
    Hqp,fqp,Ach,bch,lb,ub,qpOptions);
if exitflag<=0
    error('companionSmpcStep:Infeasible', ...
        'Chance-constrained companion QP is infeasible.');
end
centeredFirstInput = U(1:nu)-model.u_mean;
xiNext = model.Ac*xi+model.Bc*centeredFirstInput;
yPred = model.y_mean+model.Cc*xiNext;

out.A_ch = Ach;
out.b_ch = bch;
out.risk_each = riskEach;
out.z_quantile = zQuantile;
out.exitflag = exitflag;
out.lb = lb;
out.ub = ub;
out.U0 = U0;
out.cost = 0.5*U'*(Hqp-1e-9*eye(N*nu))*U+fqp'*U+costConstant;
out.prediction_state_dimension = stateDimension;
out.used_companion_prediction = true;
out.companion_order = model.companion.order;
out.estimated_sigma_eps = sqrt(trace( ...
    model.Qc(1:latentDimension,1:latentDimension))/latentDimension);
out.estimated_sigma_obs = sqrt(trace(model.Sigma_obs)/ny);
end

function [solution,objective,exitflag] = callMathWorksQuadprog( ...
    Hessian,linearTerm,Aineq,bineq,lb,ub,options)
% Avoid a broken MOSEK quadprog shim without changing the persistent path.
resolved = which('quadprog');
shadowFolder = '';
if contains(lower(resolved),'mosek')
    shadowFolder = fileparts(resolved);
    rmpath(shadowFolder);
    cleanup = onCleanup(@() addpath(shadowFolder,'-begin')); %#ok<NASGU>
    rehash path;
    resolved = which('quadprog');
end
assert(~contains(lower(resolved),'mosek'), ...
    'companionSmpcStep:QuadprogShadow', ...
    'Unable to resolve the MathWorks Optimization Toolbox quadprog.');
[solution,objective,exitflag] = feval('quadprog', ...
    Hessian,linearTerm,Aineq,bineq,[],[],lb,ub,[],options);
end
