function [E,stats]=select_hard_preference_outputs(weights,q)
%SELECT_HARD_PREFERENCE_OUTPUTS Select q original output axes by weight.
w=weights(:); p=numel(w); assert(q>=1 && q<=p && all(isfinite(w)) && all(w>=0));
[~,idx]=sort(w,'descend'); idx=idx(1:q); E=zeros(p,q);
for j=1:q, E(idx(j),j)=1; end
stats=struct('weights',w,'selected_indices',idx(:)','hard_locked',true, ...
    'uses_new_training_data',false);
end
