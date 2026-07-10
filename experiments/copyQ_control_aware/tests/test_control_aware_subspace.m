function test_control_aware_subspace
% A reduced control model must retain the tracked output directions exactly.
rng(11,'twister'); y=randn(8,300); u=randn(2,300); tracked=[1 2]; ell=4;
[A,B,P,R,S,stats] = control_aware_subspace_varx(y,u,ell,tracked);
E=zeros(8,2); E(tracked,:)=eye(2);
assert(isequal(size(P),[8,ell]));
assert(norm(P'*P-eye(ell),'fro') < 1e-10);
assert(norm(P*P'*E-E,'fro') < 1e-10, ...
    'Tracked output axes must lie in the retained subspace.');
assert(norm(R-P,'fro') < 1e-12);
assert(all(eig((S+S')/2)>-1e-10));
assert(stats.tracked_projection_error < 1e-10);
fprintf('PASS: control-aware subspace contains tracked outputs.\n');
end
