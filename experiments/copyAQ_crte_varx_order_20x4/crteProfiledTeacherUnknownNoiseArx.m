function [ABlocks,BBlocks,P,R,SigmaEps,stats] = crteProfiledTeacherUnknownNoiseArx(y,u,ell,tracked,order,opt)
%CRTEPROFILEDTEACHERUNKNOWNNOISEARX Profile the CRTE teacher with ARX(s).
% Every candidate rebuilds its metric dual, re-extracts latent coordinates,
% and fits a genuine s-lag latent ARX model.  True sensor Sigma_n is neither
% accepted nor used; the noise term is a two-fold cross-fitted ARX residual
% covariance proxy and is not a probability-safe sensor covariance.

if nargin < 6 || isempty(opt)
    opt = struct();
end
defaults = struct('mu_grid',[0 0.25 0.5 0.75 1], ...
    'alpha',1,'beta',1,'omega',[],'prediction_horizon',18, ...
    'Ru',eye(size(u,1)),'val_fraction',0.25,'ridge',1e-8, ...
    'rank_tol',1e-9,'reach_tau',1e-10,'num_random_subspaces',30, ...
    'G',[],'seed',20260725);
names = fieldnames(defaults);
for index = 1:numel(names)
    name = names{index};
    if ~isfield(opt,name) || isempty(opt.(name))
        opt.(name) = defaults.(name);
    end
end

validateattributes(order,{'numeric'},{'scalar','integer','positive'});
p = size(y,1);
nu = size(u,1);
sampleCount = size(y,2);
q = numel(tracked);
freeDimension = ell-q;
assert(size(u,2)==sampleCount,'y and u sample counts differ.');
assert(freeDimension>0 && ell<=p,'Need 0 < ell-q and ell <= p.');
assert(all(opt.mu_grid>=0 & opt.mu_grid<=1),'mu must lie in [0,1].');
assert(isequal(size(opt.Ru),[nu nu]) && ...
    min(eig((opt.Ru+opt.Ru')/2))>0,'Ru must be SPD.');
predictionHorizon = opt.prediction_horizon;
if isempty(opt.omega)
    omega = ones(1,predictionHorizon)/predictionHorizon;
else
    omega = opt.omega(:)'/sum(opt.omega);
end
assert(numel(omega)==predictionHorizon && all(omega>=0), ...
    'omega must be nonnegative and match the horizon.');

validationCount = max(50,round(opt.val_fraction*sampleCount));
trainingCount = sampleCount-validationCount;
trainingSamples = 1:trainingCount;
validationSamples = trainingCount+1:sampleCount;
assert(trainingCount>predictionHorizon+order+ell+nu+20, ...
    'Training segment is too short for ARX(s).');
trainingTransitions = order:trainingCount-1;

yMean = mean(y(:,trainingSamples),2);
uMean = mean(u(:,trainingSamples),2);
yCentered = y-yMean;
uCentered = u-uMean;
trackedBasis = zeros(p,q);
trackedBasis(tracked,:) = eye(q);
freeAmbientBasis = null(trackedBasis');
ambientFreeDimension = size(freeAmbientBasis,2);
freeOutput = freeAmbientBasis'*yCentered;
SigmaPerp = cov(freeOutput(:,trainingSamples)',1);
SigmaPerp = (SigmaPerp+SigmaPerp')/2;
scale = max(trace(SigmaPerp)/ambientFreeDimension,1e-12);
ridge = opt.ridge*max(scale,1)+1e-12;
SigmaPerp = SigmaPerp+ridge*eye(ambientFreeDimension);
if isempty(opt.G)
    G = eye(ambientFreeDimension);
else
    G = (opt.G+opt.G')/2;
end
assert(isequal(size(G),[ambientFreeDimension ambientFreeDimension]) && ...
    min(eig(G))>0,'G must be d-by-d SPD.');
tauG = trace(G\SigmaPerp)/ambientFreeDimension;

% Exact OLS FWL residual maker using the same s-lag information set.
base = buildLagRegressor(yCentered(tracked,:),uCentered,order,trainingTransitions);
W = freeOutput(:,trainingTransitions);
trackedFuture = yCentered(tracked,trainingTransitions+1);
H0 = base'*pinv(base*base')*base;
H0 = (H0+H0')/2;
M0 = eye(size(H0))-H0;
M0 = (M0+M0')/2;
Z0 = W*M0;
BT = Z0*Z0';
taskGram = trackedFuture'*trackedFuture;
AT = Z0*((taskGram+taskGram')/2)*Z0';
BT = (BT+BT')/2;
AT = (AT+AT')/2;
projectorIdempotency = norm(H0*H0-H0,'fro');

[Us,Ss,~] = svd(Z0,'econ');
singularValues = diag(Ss);
svdTolerance = opt.rank_tol*max(max(singularValues),1);
fwlRank = sum(singularValues>svdTolerance);
if fwlRank<freeDimension
    error('crteProfiledTeacherUnknownNoiseArx:InsufficientFWLRank', ...
        'FWL support rank %d is smaller than requested ell_f=%d.', ...
        fwlRank,freeDimension);
end
Us = Us(:,1:fwlRank);
singularValues = singularValues(1:fwlRank);
Qsupport = Us*diag(1./singularValues);
Psupport = Us*Us';
supportIdentityError = norm(Qsupport'*BT*Qsupport-eye(fwlRank),'fro');

SigmaNoiseProxy = crossfitFreeResidualCovArx( ...
    freeOutput,uCentered,order,trainingTransitions,ridge);
predictionContent = fixedPredictionContentArx( ...
    freeOutput,uCentered,order,trainingTransitions,ridge);
normalizedPrediction = normalizeFrobenius(predictionContent);
normalizedTask = normalizeFrobenius(AT);
normalizedNoise = normalizeFrobenius(SigmaNoiseProxy);

pool = {};
candidateCount = 0;
for metricIndex = 1:numel(opt.mu_grid)
    mu = opt.mu_grid(metricIndex);
    [metricSqrt,~,metric] = metricRoots(SigmaPerp,G,tauG,mu,ridge);
    supportMetric = Qsupport'*metric*Qsupport;
    supportMetric = (supportMetric+supportMetric')/2;
    [~,supportMetricInvSqrt] = spdRoots(supportMetric,ridge);
    for predictionWeight = [0 0.5 1]
        for taskWeight = [0 0.5 1]
            score = Qsupport'*(normalizedPrediction + ...
                predictionWeight*normalizedTask-taskWeight*normalizedNoise)*Qsupport;
            whitened = supportMetricInvSqrt*((score+score')/2)*supportMetricInvSqrt;
            whitened = (whitened+whitened')/2;
            [vectors,values] = eig(whitened);
            [~,sortOrder] = sort(real(diag(values)),'descend');
            zeta = fixColumnSigns(vectors(:,sortOrder(1:freeDimension)));
            readout = Qsupport*(supportMetricInvSqrt*zeta);
            readout = metricNormalize(readout,metric,ridge);
            loading = metric*readout;
            X = metricSqrt*readout;
            candidateCount = candidateCount+1;
            pool{candidateCount} = makeCandidate(mu,X,readout,loading, ...
                'spectral',predictionWeight,taskWeight); %#ok<AGROW>
        end
    end
end

rng(opt.seed,'twister');
for metricIndex = 1:numel(opt.mu_grid)
    mu = opt.mu_grid(metricIndex);
    [metricSqrt,~,metric] = metricRoots(SigmaPerp,G,tauG,mu,ridge);
    for randomIndex = 1:opt.num_random_subspaces
        [zeta,~] = qr(randn(fwlRank,freeDimension),0);
        zeta = fixColumnSigns(zeta);
        readout = metricNormalize(Qsupport*zeta,metric,ridge);
        loading = metric*readout;
        X = metricSqrt*readout;
        candidateCount = candidateCount+1;
        pool{candidateCount} = makeCandidate(mu,X,readout,loading, ...
            'random',NaN,NaN); %#ok<AGROW>
    end
end

rows = struct([]);
for candidateIndex = 1:candidateCount
    candidate = pool{candidateIndex};
    [candidateA,candidateB,candidateP,candidateR,~,detail] = refitCandidate( ...
        yCentered,uCentered,trackedBasis,freeAmbientBasis, ...
        candidate.loading,candidate.readout,order,trainingTransitions,ridge);
    [companionA,companionB,~,~,companionMeta] = buildArxCompanion( ...
        candidateA,candidateB,candidateP,eye(ell));

    predictionTerm = multistepFreeResidualArx( ...
        yCentered,uCentered,candidateR,candidateA,candidateB,q, ...
        trainingSamples,predictionHorizon,omega);
    denominator = candidate.readout'*BT*candidate.readout;
    denominator = (denominator+denominator')/2;
    numerator = candidate.readout'*AT*candidate.readout;
    numerator = (numerator+numerator')/2;
    denominatorEigenvalues = eig(denominator);
    fwlValid = min(denominatorEigenvalues) > ...
        opt.rank_tol*max(max(denominatorEigenvalues),1);
    if fwlValid
        taskTerm = real(trace(denominator\numerator))/size(W,2);
    else
        taskTerm = NaN;
    end
    noiseTerm = real(trace(candidate.readout'*SigmaNoiseProxy*candidate.readout));

    authority = zeros(size(companionA));
    Ru = (opt.Ru+opt.Ru')/2;
    for horizonIndex = 0:predictionHorizon-1
        reachMap = (companionA^horizonIndex)*companionB;
        authority = authority+(reachMap/Ru)*reachMap';
    end
    freeAuthority = diag(authority(q+1:ell,q+1:ell));
    minimumAuthority = min(freeAuthority);
    teacherObjective = predictionTerm-opt.alpha*taskTerm+opt.beta*noiseTerm;
    validationNrmse = validationMultistepNrmseArx( ...
        yCentered,uCentered,tracked,candidateA,candidateB, ...
        candidateP,candidateR,order,validationSamples,predictionHorizon,omega);
    feasible = fwlValid && minimumAuthority>=opt.reach_tau && ...
        detail.dual_error<1e-8 && detail.spectral_radius<1.05 && ...
        isfinite(teacherObjective);

    rows(candidateIndex).mu = candidate.mu;
    rows(candidateIndex).source = candidate.source;
    rows(candidateIndex).initializer_alpha = candidate.alpha;
    rows(candidateIndex).initializer_beta = candidate.beta;
    rows(candidateIndex).prediction_term = predictionTerm;
    rows(candidateIndex).task_term = taskTerm;
    rows(candidateIndex).noise_term = noiseTerm;
    rows(candidateIndex).teacher_objective = teacherObjective;
    rows(candidateIndex).reach_min = minimumAuthority;
    rows(candidateIndex).fwl_min_eig = min(denominatorEigenvalues);
    rows(candidateIndex).support_residual = norm( ...
        (eye(ambientFreeDimension)-Psupport)*candidate.readout,'fro');
    rows(candidateIndex).fwl_valid = fwlValid;
    rows(candidateIndex).validation_nrmse = validationNrmse;
    rows(candidateIndex).spectral_radius = detail.spectral_radius;
    rows(candidateIndex).dual_error = detail.dual_error;
    rows(candidateIndex).feasible = feasible;
    rows(candidateIndex).arx_order = order;
    rows(candidateIndex).companion_state_dimension = companionMeta.state_dimension;
    rows(candidateIndex).X = candidate.X;
    rows(candidateIndex).readout = candidate.readout;
end

feasibleMask = [rows.feasible];
assert(any(feasibleMask),'No feasible profiled-teacher ARX candidate.');
feasibleIndices = find(feasibleMask);
[~,bestPosition] = min([rows(feasibleIndices).teacher_objective]);
bestIndex = feasibleIndices(bestPosition);
best = rows(bestIndex);

% Freeze selected subspace/metric and refit ARX(s) on all offline data.
yMeanFull = mean(y,2);
uMeanFull = mean(u,2);
yCenteredFull = y-yMeanFull;
uCenteredFull = u-uMeanFull;
freeOutputFull = freeAmbientBasis'*yCenteredFull;
SigmaPerpFull = cov(freeOutputFull',1);
SigmaPerpFull = (SigmaPerpFull+SigmaPerpFull')/2+ridge*eye(ambientFreeDimension);
tauFull = trace(G\SigmaPerpFull)/ambientFreeDimension;
[~,~,metricFull] = metricRoots(SigmaPerpFull,G,tauFull,best.mu,ridge);
readout = metricNormalize(best.readout,metricFull,ridge);
loading = metricFull*readout;
allTransitions = order:sampleCount-1;
[ABlocks,BBlocks,P,R,SigmaEps,detail] = refitCandidate( ...
    yCenteredFull,uCenteredFull,trackedBasis,freeAmbientBasis, ...
    loading,readout,order,allTransitions,ridge);
Pbar = null(R');
completePrimal = [P,Pbar];
assert(rcond(completePrimal)>1e-12,'Complete dual basis is ill-conditioned.');
completeDual = inv(completePrimal');
Rbar = completeDual(:,ell+1:end);
[companionA,~,~,~,companionMeta] = buildArxCompanion( ...
    ABlocks,BBlocks,P,SigmaEps);

stats.y_mean = yMeanFull;
stats.u_mean = uMeanFull;
stats.rows = rows;
stats.best_index = bestIndex;
stats.selected_mu = best.mu;
stats.selected_source = best.source;
stats.selected_teacher_objective = best.teacher_objective;
stats.selected_prediction_term = best.prediction_term;
stats.selected_task_term = best.task_term;
stats.selected_noise_term = best.noise_term;
stats.selected_reach_min = best.reach_min;
stats.selected_validation_nrmse = best.validation_nrmse;
stats.num_candidates = candidateCount;
stats.num_feasible = sum(feasibleMask);
stats.alpha = opt.alpha;
stats.beta = opt.beta;
stats.omega = omega;
stats.horizon = predictionHorizon;
stats.arx_order = order;
stats.companion_state_dimension = companionMeta.state_dimension;
stats.uses_true_Sigma_n=false;
stats.noise_object = 'two-fold cross-fitted free-output ARX residual covariance proxy';
stats.search_claim = 'finite candidate verification; no Stiefel global-optimum claim';
stats.H0_idempotency_error = projectorIdempotency;
stats.B_T = BT;
stats.A_T = AT;
stats.Z0 = Z0;
stats.fwl_support_rank = fwlRank;
stats.fwl_support_tol = svdTolerance;
stats.fwl_support_basis = Us;
stats.fwl_support_singular_values = singularValues;
stats.support_B_identity_error = supportIdentityError;
stats.max_candidate_support_residual = max([rows.support_residual]);
stats.Sigma_noise_proxy = SigmaNoiseProxy;
stats.C_mu = metricFull;
stats.dual_error = norm(R'*P-eye(ell),'fro');
stats.tracked_right_error = norm(P*R'*trackedBasis-trackedBasis,'fro');
stats.tracked_left_error = norm(trackedBasis'*P*R'-trackedBasis','fro');
stats.Pi_idempotency_error = norm(P*R'*P*R'-P*R','fro');
stats.dual_errors_4piece = [stats.dual_error,norm(R'*Pbar,'fro'), ...
    norm(Rbar'*P,'fro'),norm(Rbar'*Pbar-eye(p-ell),'fro')];
stats.spectral_radius = max(abs(eig(companionA)));
stats.detail = detail;
end

function candidate = makeCandidate(mu,X,readout,loading,source,alpha,beta)
candidate = struct('mu',mu,'X',X,'readout',readout, ...
    'loading',loading,'source',source,'alpha',alpha,'beta',beta);
end

function regressor = buildLagRegressor(z,u,order,transitions)
ell = size(z,1);
nu = size(u,1);
regressor = zeros(order*(ell+nu),numel(transitions));
for lag = 1:order
    sampleAtLag = transitions-lag+1;
    regressor((lag-1)*ell+(1:ell),:) = z(:,sampleAtLag);
    inputRows = order*ell+(lag-1)*nu+(1:nu);
    regressor(inputRows,:) = u(:,sampleAtLag);
end
end

function covariance = crossfitFreeResidualCovArx(z,u,order,transitions,ridge)
fold = mod(1:numel(transitions),2)+1;
residuals = [];
for foldIndex = 1:2
    fitTransitions = transitions(fold~=foldIndex);
    testTransitions = transitions(fold==foldIndex);
    [A,B] = fitLatentArx(z,u,order,ridge,fitTransitions);
    phi = buildLagRegressor(z,u,order,testTransitions);
    theta = blocksToTheta(A,B);
    errors = z(:,testTransitions+1)-theta*phi;
    residuals = [residuals,errors]; %#ok<AGROW>
end
covariance = residuals*residuals'/max(size(residuals,2),1);
covariance = (covariance+covariance')/2+ridge*eye(size(covariance));
end

function content = fixedPredictionContentArx(z,u,order,transitions,ridge)
[A,B] = fitLatentArx(z,u,order,ridge,transitions);
phi = buildLagRegressor(z,u,order,transitions);
prediction = blocksToTheta(A,B)*phi;
content = prediction*prediction'/size(prediction,2);
content = (content+content')/2;
end

function theta = blocksToTheta(ABlocks,BBlocks)
order = size(ABlocks,3);
theta = [];
for lag = 1:order
    theta = [theta,ABlocks(:,:,lag)]; %#ok<AGROW>
end
for lag = 1:order
    theta = [theta,BBlocks(:,:,lag)]; %#ok<AGROW>
end
end

function [A,B,P,R,SigmaEps,detail] = refitCandidate( ...
    yCentered,uCentered,trackedBasis,freeAmbientBasis,loading,readout, ...
    order,transitions,ridge)
q = size(trackedBasis,2);
ell = q+size(loading,2);
P = [trackedBasis,freeAmbientBasis*loading];
R = [trackedBasis,freeAmbientBasis*readout];
assert(norm(R'*P-eye(ell),'fro')<1e-7,'Metric dual failed.');
z = R'*yCentered;
[A,B,SigmaEps,fitStats] = fitLatentArx( ...
    z,uCentered,order,ridge,transitions);
[companionA,~,~,~,meta] = buildArxCompanion(A,B,P,SigmaEps);
detail.dual_error = norm(R'*P-eye(ell),'fro');
detail.spectral_radius = max(abs(eig(companionA)));
detail.fit = fitStats;
detail.companion = meta;
end

function objective = multistepFreeResidualArx( ...
    yCentered,uCentered,R,A,B,q,samples,horizon,omega)
z = R'*yCentered;
order = size(A,3);
Pidentity = eye(size(z,1));
[Ac,Bc,~,~,~] = buildArxCompanion(A,B,Pidentity,eye(size(z,1)));
lastTransition = samples(end)-horizon;
startTransition = max(samples(1)+order-1,order);
total = 0;
originCount = 0;
for transition = startTransition:lastTransition
    xi = makeCompanionState(z,uCentered,order,transition);
    for step = 1:horizon
        xi = Ac*xi+Bc*uCentered(:,transition+step-1);
        errorFree = z(q+1:end,transition+step)-xi(q+1:size(z,1));
        total = total+omega(step)*(errorFree'*errorFree);
    end
    originCount = originCount+1;
end
objective = total/max(originCount,1);
end

function value = validationMultistepNrmseArx( ...
    yCentered,uCentered,tracked,A,B,P,R,order,samples,horizon,omega)
z = R'*yCentered;
[Ac,Bc,Cc,~,~] = buildArxCompanion(A,B,P,eye(size(A,1)));
lastTransition = samples(end)-horizon;
startTransition = max(samples(1)+order-1,order);
squaredError = zeros(numel(tracked),1);
originCount = 0;
for transition = startTransition:lastTransition
    xi = makeCompanionState(z,uCentered,order,transition);
    for step = 1:horizon
        xi = Ac*xi+Bc*uCentered(:,transition+step-1);
        prediction = Cc*xi;
        errorTracked = yCentered(tracked,transition+step)-prediction(tracked);
        squaredError = squaredError+omega(step)*(errorTracked.^2);
    end
    originCount = originCount+1;
end
rmse = sqrt(squaredError/max(originCount,1));
scale = std(yCentered(tracked,samples),0,2);
value = mean(rmse./max(scale,1e-8));
end

function xi = makeCompanionState(z,u,order,transition)
latentHistory = [];
for lag = 0:order-1
    latentHistory = [latentHistory;z(:,transition-lag)]; %#ok<AGROW>
end
inputHistory = [];
for lag = 1:order-1
    inputHistory = [inputHistory;u(:,transition-lag)]; %#ok<AGROW>
end
xi = [latentHistory;inputHistory];
end

function normalized = normalizeFrobenius(matrix)
matrix = (matrix+matrix')/2;
scale = norm(matrix,'fro');
if scale>1e-12
    normalized = matrix/scale;
else
    normalized = zeros(size(matrix));
end
end

function normalized = metricNormalize(candidate,metric,ridge)
gram = candidate'*metric*candidate;
gram = (gram+gram')/2;
[~,inverseSqrt] = spdRoots(gram,ridge);
normalized = candidate*inverseSqrt;
end

function [sqrtMatrix,inverseSqrt] = spdRoots(matrix,ridge)
[vectors,values] = eig((matrix+matrix')/2);
eigenvalues = max(real(diag(values)),ridge);
sqrtMatrix = vectors*diag(sqrt(eigenvalues))*vectors';
inverseSqrt = vectors*diag(1./sqrt(eigenvalues))*vectors';
end

function [sqrtMetric,inverseSqrtMetric,metric] = metricRoots( ...
    covariance,G,tau,mu,ridge)
metric = (1-mu)*covariance+mu*tau*G;
metric = (metric+metric')/2+ridge*eye(size(metric));
[sqrtMetric,inverseSqrtMetric] = spdRoots(metric,ridge);
end

function matrix = fixColumnSigns(matrix)
for column = 1:size(matrix,2)
    [~,row] = max(abs(matrix(:,column)));
    if matrix(row,column)<0
        matrix(:,column) = -matrix(:,column);
    end
end
end
