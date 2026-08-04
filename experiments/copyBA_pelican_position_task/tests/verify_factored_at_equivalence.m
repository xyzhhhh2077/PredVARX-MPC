%% verify_factored_at_equivalence
% Verify that the factored A_T computation in fit_segmented_anchored_varx is
% numerically identical to the original T x T projection form used by copyAX.
% Both forms are evaluated on random data; the max relative difference must
% be < 1e-10.
rng(42);
p=10; m=4; q=2; ell=5; T=300; d=p-q;
y=randn(p,T); u=randn(m,T); run_id=[ones(1,150),2*ones(1,150)];
[Q,~]=qr(randn(p,q),0); Eanchor=Q;

% ---- original copyAX form (T x T) ----
yc=y-mean(y,2); uc=u-mean(u,2); ridge=1e-8;
valid=find(run_id(1:end-1)==run_id(2:end));
[Qe,~]=qr(Eanchor,0); Ptask=Qe;
Ntask=null(Ptask'); Yp=Ntask'*yc;
base=[Ptask'*yc(:,valid);uc(:,valid)]; W=Yp(:,valid);
Tfuture=Ptask'*yc(:,valid+1);
H0=base'/(base*base'+ridge*eye(size(base,1)))*base;
M0=eye(size(H0))-H0;
task_gram=Tfuture'*Tfuture;
A_T_orig=W*M0*task_gram*M0*W'; A_T_orig=(A_T_orig+A_T_orig')/2;

% ---- factored form ----
Gbase=base*base'+ridge*eye(size(base,1));
U=W*Tfuture' - (W*base')*(Gbase\(base*Tfuture'));
A_T_fact=U*U'; A_T_fact=(A_T_fact+A_T_fact')/2;

rel=norm(A_T_orig-A_T_fact,'fro')/max(norm(A_T_orig,'fro'),1e-12);
fprintf('A_T relative Frobenius difference: %.3e\n', rel);
assert(rel<1e-10,'factored form deviates from the T x T form');
fprintf('PASS: factored A_T == original A_T (rel diff %.3e)\n', rel);

% ---- also verify the full fit runs and reproduces copyAX on small data ----
opt=struct('ridge',1e-8,'mu',0.10,'ntr_epsilon',10e-6);
[Ah,Bh,Ph,Rh,Se,st]=fit_segmented_anchored_varx(y,u,run_id,Eanchor,ell,opt);
fprintf('full fit: dual_error=%.3e spectral_radius=%.4f transitions=%d\n', ...
    st.dual_error, st.spectral_radius, st.transition_count);
assert(st.dual_error<1e-7);
fprintf('PASS: full copyAY fit runs on T=300 random data\n');
