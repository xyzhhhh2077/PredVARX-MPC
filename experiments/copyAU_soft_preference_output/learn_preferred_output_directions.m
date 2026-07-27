function [E,stats] = learn_preferred_output_directions(y,u,q,opt)
%LEARN_PREFERRED_OUTPUT_DIRECTIONS Weighted output preference plus authority.
if ~isfield(opt,'weights'), error('weights are required'); end
if ~isfield(opt,'preference_strength'), opt.preference_strength=0.7; end
if ~isfield(opt,'reach_horizon'), opt.reach_horizon=18; end
if ~isfield(opt,'Ru'), opt.Ru=eye(size(u,1)); end
if ~isfield(opt,'ridge'), opt.ridge=1e-8; end
[p,T]=size(y); m=size(u,1); w=opt.weights(:);
assert(numel(w)==p && all(w>=0) && any(w>0),'weights must be nonnegative p-vector.');
assert(opt.preference_strength>=0 && opt.preference_strength<=1,'strength must be in [0,1].');
yc=y-mean(y,2); uc=u-mean(u,2);
scale=max(trace(yc*yc'/T)/p,1e-12); ridge=opt.ridge*max(scale,1)+1e-12;
Cy=yc*yc'/T+ridge*eye(p); Cy=(Cy+Cy')/2;
Phi=[yc(:,1:end-1);uc(:,1:end-1)]; Ycur=yc(:,2:end);
Theta=(Phi*Phi'+ridge*eye(p+m))\(Phi*Ycur');
Ay=Theta(1:p,:)'; By=Theta(p+1:end,:)'; Wy=zeros(p);
for h=0:opt.reach_horizon-1
    G=(Ay^h)*By; Wy=Wy+G*(opt.Ru\G');
end
Wy=(Wy+Wy')/2; Wpref=diag(w/max(w));
% Normalize both criteria before convex combination.
Wyn=Wy/max(norm(Wy,'fro'),ridge); Wpn=Wpref/max(norm(Wpref,'fro'),ridge);
M=opt.preference_strength*Wpn+(1-opt.preference_strength)*Wyn; M=(M+M')/2;
[Uc,Dc]=eig(Cy); dc=max(real(diag(Dc)),ridge); Cih=Uc*diag(1./sqrt(dc))*Uc';
K=Cih*M*Cih; K=(K+K')/2; [Z,D]=eig(K);
[val,ord]=sort(real(diag(D)),'descend'); Z=Z(:,ord);
for j=1:size(Z,2), [~,i]=max(abs(Z(:,j))); if Z(i,j)<0, Z(:,j)=-Z(:,j); end, end
Eraw=Cih*Z(:,1:q); [E,~]=qr(Eraw,0);
contribution=sum(E.^2,2); preference_capture=sum(w.*contribution)/sum(w);
stats=struct('weights',w,'preference_strength',opt.preference_strength, ...
    'contribution',contribution,'preference_capture',preference_capture, ...
    'score_values',val(1:q),'C_y',Cy,'W_y',Wy,'uses_new_training_data',false);
end
