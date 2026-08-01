function [E,stats] = learn_segmented_preferred_output_directions(y,u,run_id,q,opt)
%LEARN_SEGMENTED_PREFERRED_OUTPUT_DIRECTIONS copyAU objective on run segments.
if ~isfield(opt,'weights'), error('weights are required'); end
if ~isfield(opt,'preference_strength'), opt.preference_strength=0.7; end
if ~isfield(opt,'reach_horizon'), opt.reach_horizon=18; end
if ~isfield(opt,'Ru'), opt.Ru=eye(size(u,1)); end
if ~isfield(opt,'ridge'), opt.ridge=1e-8; end
[p,T]=size(y); m=size(u,1); w=opt.weights(:); run_id=run_id(:)';
assert(numel(run_id)==T,'run_id must match the sample count.');
assert(numel(w)==p && all(w>=0) && any(w>0),'weights must be nonnegative p-vector.');
assert(opt.preference_strength>=0 && opt.preference_strength<=1,'strength must be in [0,1].');
valid=find(run_id(1:end-1)==run_id(2:end));
yc=y-mean(y,2); uc=u-mean(u,2);
scale=max(trace(yc*yc'/T)/p,1e-12); ridge=opt.ridge*max(scale,1)+1e-12;
Cy=yc*yc'/T+ridge*eye(p); Cy=(Cy+Cy')/2;
Phi=[yc(:,valid);uc(:,valid)]; Ycur=yc(:,valid+1);
Theta=(Phi*Phi'+ridge*eye(p+m))\(Phi*Ycur');
Ay=Theta(1:p,:)'; By=Theta(p+1:end,:)'; Wy=zeros(p);
for h=0:opt.reach_horizon-1
    G=(Ay^h)*By; Wy=Wy+G*(opt.Ru\G');
end
Wy=(Wy+Wy')/2; Wpref=diag(w/max(w));
Wyn=Wy/max(norm(Wy,'fro'),ridge); Wpn=Wpref/max(norm(Wpref,'fro'),ridge);
M=opt.preference_strength*Wpn+(1-opt.preference_strength)*Wyn; M=(M+M')/2;
[Uc,Dc]=eig(Cy); dc=max(real(diag(Dc)),ridge); Cih=Uc*diag(1./sqrt(dc))*Uc';
K=Cih*M*Cih; K=(K+K')/2; [Z,D]=eig(K);
[val,ord]=sort(real(diag(D)),'descend'); Z=deterministic_sign(Z(:,ord));
Eraw=Cih*Z(:,1:q); [E,~]=qr(Eraw,0); E=deterministic_sign(E);
contribution=sum(E.^2,2);
stats=struct('weights',w,'preference_strength',opt.preference_strength, ...
    'contribution',contribution,'preference_capture',sum(w.*contribution)/sum(w), ...
    'score_values',val(1:q),'C_y',Cy,'W_y',Wy,'uses_new_training_data',true, ...
    'transition_count',numel(valid),'segment_count',numel(unique(run_id)));
end

function X=deterministic_sign(X)
for j=1:size(X,2)
    [~,i]=max(abs(X(:,j)));
    if X(i,j)<0, X(:,j)=-X(:,j); end
end
end