function [Ac, Bc, Cc, Qc, meta] = buildArxCompanion(ABlocks, BBlocks, P, SigmaEps)
%BUILDARXCOMPANION Exact state-space realization of a latent ARX(s) model.
% The ARX convention is
%   z(k+1) = sum_i A_i z(k+1-i) + sum_i B_i v(k+1-i),
% where v=u-u_mean.  The companion state is
%   xi(k)=[z(k);...;z(k-s+1);v(k-1);...;v(k-s+1)].

if ismatrix(ABlocks)
    ABlocks = reshape(ABlocks, size(ABlocks,1), size(ABlocks,2), 1);
end
if ismatrix(BBlocks)
    BBlocks = reshape(BBlocks, size(BBlocks,1), size(BBlocks,2), 1);
end
ell = size(ABlocks,1);
order = size(ABlocks,3);
nu = size(BBlocks,2);
ny = size(P,1);

assert(size(ABlocks,2) == ell, 'Each A block must be ell-by-ell.');
assert(size(BBlocks,1) == ell && size(BBlocks,3) == order, ...
    'B blocks must be ell-by-nu-by-order.');
assert(size(P,2) == ell, 'P must have ell columns.');
assert(isequal(size(SigmaEps), [ell ell]), ...
    'SigmaEps must be ell-by-ell.');

stateDimension = order*ell + (order-1)*nu;
Ac = zeros(stateDimension);
Bc = zeros(stateDimension, nu);

for lag = 1:order
    zCols = (lag-1)*ell + (1:ell);
    Ac(1:ell,zCols) = ABlocks(:,:,lag);
end
Bc(1:ell,:) = BBlocks(:,:,1);

if order > 1
    inputStart = order*ell;
    for lag = 2:order
        inputCols = inputStart + (lag-2)*nu + (1:nu);
        Ac(1:ell,inputCols) = BBlocks(:,:,lag);
    end

    % Shift latent histories: z(k), z(k-1), ...
    latentShiftRows = ell + (1:(order-1)*ell);
    latentShiftCols = 1:(order-1)*ell;
    Ac(latentShiftRows,latentShiftCols) = eye((order-1)*ell);

    % Insert v(k) and shift previous centered inputs.
    firstInputRows = inputStart + (1:nu);
    Bc(firstInputRows,:) = eye(nu);
    if order > 2
        lowerInputRows = inputStart + nu + (1:(order-2)*nu);
        previousInputCols = inputStart + (1:(order-2)*nu);
        Ac(lowerInputRows,previousInputCols) = eye((order-2)*nu);
    end
end

Cc = [P, zeros(ny,stateDimension-ell)];
Qc = zeros(stateDimension);
Qc(1:ell,1:ell) = (SigmaEps+SigmaEps')/2;

meta.order = order;
meta.latent_dimension = ell;
meta.input_dimension = nu;
meta.output_dimension = ny;
meta.state_dimension = stateDimension;
meta.state_layout = '[z(k:-1:k-s+1); v(k-1:-1:k-s+1)]';
end
