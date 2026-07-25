function [ABlocks, BBlocks, SigmaEps, stats] = fitLatentArx(z, u, order, ridge, transitionIndices)
%FITLATENTARX Ridge-regularized latent ARX(s) regression.
% Fits
%   z(k+1)=sum_{i=1}^s A_i z(k+1-i)+sum_{i=1}^s B_i v(k+1-i)+eps(k)
% on explicitly supplied transition indices k.  Inputs must already be
% centered when centering is required by the caller.

if nargin < 4 || isempty(ridge)
    ridge = 1e-8;
end
sampleCount = size(z,2);
if nargin < 5 || isempty(transitionIndices)
    transitionIndices = order:sampleCount-1;
end
validateattributes(order, {'numeric'}, {'scalar','integer','positive'});
assert(size(u,2) == sampleCount, 'z and u sample counts must match.');
transitionIndices = transitionIndices(:)';
assert(~isempty(transitionIndices), 'At least one transition is required.');
assert(all(transitionIndices >= order & transitionIndices <= sampleCount-1), ...
    'Transitions must satisfy order <= k <= sampleCount-1.');

ell = size(z,1);
nu = size(u,1);
numTransitions = numel(transitionIndices);
PhiZ = zeros(order*ell,numTransitions);
PhiU = zeros(order*nu,numTransitions);
for lag = 1:order
    sampleAtLag = transitionIndices-lag+1;
    PhiZ((lag-1)*ell+(1:ell),:) = z(:,sampleAtLag);
    PhiU((lag-1)*nu+(1:nu),:) = u(:,sampleAtLag);
end
Phi = [PhiZ;PhiU];
Target = z(:,transitionIndices+1);
Gram = Phi*Phi';
scale = max(trace(Gram)/max(size(Gram,1),1),1);
regularizer = ridge*scale;
Theta = (Gram+regularizer*eye(size(Gram)))\(Phi*Target');

ABlocks = zeros(ell,ell,order);
BBlocks = zeros(ell,nu,order);
for lag = 1:order
    zRows = (lag-1)*ell+(1:ell);
    uRows = order*ell+(lag-1)*nu+(1:nu);
    ABlocks(:,:,lag) = Theta(zRows,:)';
    BBlocks(:,:,lag) = Theta(uRows,:)';
end
Residual = Target-Theta'*Phi;
SigmaEps = (Residual*Residual')/max(numTransitions,1);
SigmaEps = (SigmaEps+SigmaEps')/2;

stats.order = order;
stats.num_transitions = numTransitions;
stats.regressor_dimension = order*(ell+nu);
stats.design_rank = rank(Phi);
stats.transition_indices = transitionIndices;
stats.ridge = ridge;
stats.regularizer = regularizer;
stats.residual_rms = sqrt(mean(Residual(:).^2));
stats.residuals = Residual;
end
