function test_crte_profiled_teacher_unknown_noise
here=fileparts(mfilename('fullpath')); addpath(fileparts(here)); rng(2207,'twister');
p=10; m=2; ell=5; tracked=[1 2]; T=650;
F=diag([0.92 0.79 0.61]); Gp=[0.40 -0.12;0.15 0.30;-0.05 0.22]; C=randn(p,3); C(1,:)=[1 0 0]; C(2,:)=[0 1 0];
u=randn(m,T); x=zeros(3,T+1); y=zeros(p,T);
D=diag(linspace(.03,.30,p)); L=chol(D*D,'lower');
for k=1:T, y(:,k)=C*x(:,k)+L*randn(p,1); x(:,k+1)=F*x(:,k)+Gp*u(:,k)+.03*randn(3,1); end
opt=struct('mu_grid',[0 .5 1],'alpha',.5,'beta',.5,'prediction_horizon',5, ...
    'Ru',diag([.2 .8]),'num_random_subspaces',4,'seed',77,'reach_tau',1e-12);
[A,B,P,R,S,st]=crte_profiled_teacher_unknown_noise(y,u,ell,tracked,opt);

%% Contract / no oracle Sigma_n
src=fileread(fullfile(fileparts(here),'crte_profiled_teacher_unknown_noise.m'));
assert(~contains(regexp(src,'function[^\n]+','match','once'),'Sigma_n'),'Signature must not accept Sigma_n.');
assert(st.uses_true_Sigma_n==false);
assert(contains(st.noise_object,'cross-fitted'));

%% Exact OLS FWL projector, compact-SVD support, and domain
assert(st.H0_idempotency_error<1e-8,'H0 is not an OLS projector.');
assert(st.fwl_support_rank >= ell-numel(tracked));
assert(st.support_B_identity_error < 1e-8,'Qsupport''*B_T*Qsupport is not identity.');
assert(st.max_candidate_support_residual < 1e-8,'Candidate left the compact-SVD FWL support.');
assert(all(arrayfun(@(r) r.fwl_valid || isnan(r.task_term),st.rows)));
assert(all([st.rows.fwl_valid]),'Support-generated candidates must have positive FWL denominator.');

%% Complete candidate-dependent reconstruction/refit
assert(norm(R'*P-eye(ell),'fro')<1e-8);
assert(max(st.dual_errors_4piece)<1e-8);
assert(st.num_candidates==numel(opt.mu_grid)*(9+opt.num_random_subspaces));
assert(st.num_feasible>0);
vals=[st.rows.teacher_objective]; vals=vals(isfinite(vals));
assert(range(vals)>1e-8,'Teacher objective does not vary across candidates.');
rhos=[st.rows.spectral_radius]; assert(range(rhos)>1e-8,'VARX refit appears candidate-invariant.');

%% Teacher decomposition identity
br=st.rows(st.best_index);
assert(abs(br.teacher_objective-(br.prediction_term-opt.alpha*br.task_term+opt.beta*br.noise_term))<1e-10);
assert(br.reach_min>=opt.reach_tau);
assert(min(eig((S+S')/2))>-1e-8);
assert(st.spectral_radius<1.05);

%% Ru^{-1} must affect authority scale
opt2=opt; opt2.Ru=2*opt.Ru;
[~,~,~,~,~,st2]=crte_profiled_teacher_unknown_noise(y,u,ell,tracked,opt2);
% Same candidate pool/VARX; doubling Ru should halve every finite-horizon authority.
r1=[st.rows.reach_min]; r2=[st2.rows.reach_min];
assert(max(abs(r2-.5*r1)./max(abs(r1),1e-12))<1e-6,'Ru^{-1} is missing or incorrectly applied.');

%% Rank-deficient support must reject ell_f larger than effective FWL rank
% Construct free data with only one residualized free direction, but request
% ell_f=3. Sec. 5.3 requires an explicit rank rejection, not denominator ridge.
p2=7; tracked2=[1 2]; ell2=5; T2=300; u2=randn(2,T2);
t=(1:T2); latent=sin(.03*t)+.1*randn(1,T2);
y2=zeros(p2,T2); y2(1,:)=.2*latent; y2(2,:)=.1*cos(.02*t);
y2(3:end,:)=repmat(latent,p2-2,1); % free block rank 1
caught=false;
try
    crte_profiled_teacher_unknown_noise(y2,u2,ell2,tracked2, ...
        struct('mu_grid',0,'prediction_horizon',3,'num_random_subspaces',0,'rank_tol',1e-8));
catch ME
    caught=strcmp(ME.identifier,'crte_profiled_teacher_unknown_noise:InsufficientFWLRank');
end
assert(caught,'Rank-deficient FWL support did not reject ell_f > effective rank.');

fprintf('PASS profiled teacher unknown-noise + Sec5.3 SVD support: rank=%d support-I=%.2e support-res=%.2e candidates=%d feasible=%d selected=%d J=%.6g [pred=%.6g task=%.6g noise=%.6g] reach=%.3e dual=%.2e rho=%.4f val=%.4f\n', ...
    st.fwl_support_rank,st.support_B_identity_error,st.max_candidate_support_residual, ...
    st.num_candidates,st.num_feasible,st.best_index,st.selected_teacher_objective, ...
    st.selected_prediction_term,st.selected_task_term,st.selected_noise_term, ...
    st.selected_reach_min,st.dual_error,st.spectral_radius,st.selected_validation_nrmse);
end
