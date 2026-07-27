function [E,stats] = learn_output_directions(y,u,q,opt)
%LEARN_OUTPUT_DIRECTIONS Learn q-dimensional output directions from old data.
% mode='supervised': preserve the span declared by task_outputs, then use
%                     data covariance only to normalize its coordinates.
% mode='authority' : learn directions with maximum finite-horizon input
%                     authority under the output covariance metric.
if nargin < 4 || isempty(opt), opt = struct(); end
if ~isfield(opt,'mode'), opt.mode = 'authority'; end
if ~isfield(opt,'task_outputs'), opt.task_outputs = []; end
if ~isfield(opt,'reach_horizon'), opt.reach_horizon = 18; end
if ~isfield(opt,'Ru') || isempty(opt.Ru), opt.Ru = eye(size(u,1)); end
if ~isfield(opt,'ridge'), opt.ridge = 1e-8; end

[p,T] = size(y); m = size(u,1);
assert(size(u,2)==T,'y and u must have equal sample counts.');
assert(q>=1 && q<=p,'q must satisfy 1 <= q <= p.');
Ru = (opt.Ru+opt.Ru')/2;
assert(isequal(size(Ru),[m m]) && min(eig(Ru))>0,'Ru must be SPD.');

yc = y-mean(y,2); uc = u-mean(u,2);
scale = max(trace(yc*yc'/T)/p,1e-12);
ridge = opt.ridge*max(scale,1)+1e-12;
Cy = yc*yc'/T + ridge*eye(p); Cy = (Cy+Cy')/2;

Ylag = yc(:,1:end-1); Ycur = yc(:,2:end); Ulag = uc(:,1:end-1);
Phi = [Ylag;Ulag];
Theta = (Phi*Phi'+ridge*eye(p+m))\(Phi*Ycur');
Ay = Theta(1:p,:)'; By = Theta(p+1:end,:)';
Wy = zeros(p); Rui = Ru\eye(m);
for h = 0:opt.reach_horizon-1
    Gh = (Ay^h)*By;
    Wy = Wy+Gh*Rui*Gh';
end
Wy = (Wy+Wy')/2;

switch lower(opt.mode)
    case 'supervised'
        F = opt.task_outputs;
        assert(isequal(size(F),[p q]) && rank(F)==q, ...
            'task_outputs must be p-by-q and full rank.');
        [Qf,~] = qr(F,0);
        E = Qf/real_sqrtm(Qf'*Cy*Qf,ridge);
        score_values = eig((E'*Wy*E+E'*Wy'*E)/2);
        method = 'supervised output-span learner';
    case 'authority'
        [Uc,Dc] = eig(Cy); dc = max(real(diag(Dc)),ridge);
        Cih = Uc*diag(1./sqrt(dc))*Uc';
        M = Cih*Wy*Cih; M = (M+M')/2;
        [Z,D] = eig(M); [all_values,ord] = sort(real(diag(D)),'descend');
        Z = deterministic_sign(Z(:,ord));
        E = Cih*Z(:,1:q);
        E = E/real_sqrtm(E'*Cy*E,ridge);
        score_values = all_values(1:q);
        method = 'finite-horizon input-authority learner';
    otherwise
        error('Unknown mode: %s',opt.mode);
end

Panchor = E;
Ranchor = Cy*E;
dual_error = norm(Ranchor'*Panchor-eye(q),'fro');
assert(dual_error<1e-7,'Learned output anchor is not dual.');
[~,Dall] = eig((Wy+Wy')/2,Cy);
all_authority = sort(max(real(diag(Dall)),0),'descend');
captured = sum(max(real(eig((E'*Wy*E+E'*Wy'*E)/2)),0));
total = sum(all_authority);
if total<=1e-12, authority_fraction=0; else, authority_fraction=captured/total; end

stats = struct();
stats.method = method;
stats.mode = lower(opt.mode);
stats.C_y = Cy; stats.A_y = Ay; stats.B_y = By; stats.W_y = Wy;
stats.P_anchor = Panchor; stats.R_anchor = Ranchor;
stats.dual_error = dual_error; stats.score_values = score_values;
stats.authority_fraction = authority_fraction;
stats.uses_new_training_data = false;
stats.reach_horizon = opt.reach_horizon;
end

function S = real_sqrtm(A,ridge)
A=(A+A')/2; [U,D]=eig(A); d=max(real(diag(D)),ridge);
S=U*diag(sqrt(d))*U';
end

function X = deterministic_sign(X)
for j=1:size(X,2)
    [~,i]=max(abs(X(:,j)));
    if X(i,j)<0, X(:,j)=-X(:,j); end
end
end
