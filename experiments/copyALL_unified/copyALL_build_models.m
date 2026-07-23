function models=copyALL_build_models(D,cfg)
% Restore H/K/O/P/Q/R/S/T/U/V/X/Y/Z on one common dataset.
y=D.y_off; u=D.u_off; E=eye_cols(cfg.p,cfg.tracked);
[A,B,P,~,R,~,~,S]=predvarx_identify_oblique(y,u,cfg.ell,1,1,cfg.A,cfg.B,D.C,cfg.n,cfg.m,cfg.p);
st=base_stats(y,u,P,R,E); models.O=pack(A,B,P,R,S,st,'O: free oblique IVR');
[A,B,P,R,S,st]=control_ready_subspace_varx(y,u,cfg.ell); st=complete_stats(st,y,u,P,R,E); models.P=pack(A,B,P,R,S,st,'P: centered PCA/SVD orthogonal');
[A,B,P,R,S,st]=control_aware_iterative_ivr_varx(y,u,cfg.ell,cfg.tracked); st=complete_stats(st,y,u,P,R,E); models.Q=pack(A,B,P,R,S,st,'Q: tracked-axis control-aware IVR');
[A,B,P,~,R,~,~,S,~,~,~,~,info]=predvarx_identify_strict_whitened(y,u,cfg.ell); st=base_stats(y,u,P,R,E); st.whitening_applied=info.whitening_applied; st.complete_algorithm1=info.complete_algorithm1; models.R=pack(A,B,P,R,S,st,'R: rank-full s=1 whitened IVR');
models.H=models.Q; models.H.description='H: sliding-window full structure re-identification'; models.H.online_reidentify=true; models.H.online_covariance_update=true; models.H.reidentify_period=cfg.hk_reidentify_period; models.H.reidentify_recent=cfg.hk_reidentify_window; models.H.reidentify_mode='sliding'; models.H.covariance_mode='rolling'; models.H.bias_correction=false; models.H.bias_gamma=0; models.H.dinkla_window=cfg.noise_window;
models.K=models.Q; models.K.description='K: cumulative re-identification + Dinkla covariance + bias correction'; models.K.online_reidentify=true; models.K.online_covariance_update=true; models.K.reidentify_period=cfg.hk_reidentify_period; models.K.reidentify_recent=inf; models.K.reidentify_mode='cumulative'; models.K.covariance_mode='dinkla'; models.K.bias_correction=true; models.K.bias_gamma=cfg.hk_bias_gamma; models.K.dinkla_window=cfg.hk_dinkla_window;
models.S=models.Q; models.S.description='S: periodic offline+recent structure re-identification of Q geometry'; models.S.online_reidentify=true; models.S.online_covariance_update=true; models.S.reidentify_period=cfg.reidentify_period; models.S.reidentify_recent=cfg.reidentify_recent; models.S.reidentify_mode='offline_plus_recent'; models.S.covariance_mode='rolling'; models.S.bias_correction=false; models.S.bias_gamma=0; models.S.dinkla_window=cfg.noise_window;
models.T=models.Q; models.T.description='T: fixed Q geometry with online covariance updates only'; models.T.online_reidentify=false; models.T.online_covariance_update=true; models.T.reidentify_period=inf; models.T.reidentify_recent=0; models.T.reidentify_mode='fixed'; models.T.covariance_mode='rolling'; models.T.bias_correction=false; models.T.bias_gamma=0; models.T.dinkla_window=cfg.noise_window;
for n={'O','P','Q','R'}, models.(n{1}).online_reidentify=false; models.(n{1}).online_covariance_update=true; models.(n{1}).reidentify_period=inf; models.(n{1}).reidentify_recent=0; models.(n{1}).reidentify_mode='fixed'; models.(n{1}).covariance_mode='rolling'; models.(n{1}).bias_correction=false; models.(n{1}).bias_gamma=0; models.(n{1}).dinkla_window=cfg.noise_window; end

% U: one-shot control-aware SVD/PCA complement. This is the smooth-noise
% process copy's identifier placed on the exact same common data/evaluator.
[A,B,P,R,S,st]=control_aware_subspace_varx(y,u,cfg.ell,cfg.tracked);
st=complete_stats(st,y,u,P,R,E);
models.U=pack(A,B,P,R,S,st,'U: control-aware one-shot SVD complement');

% V: iterative predictable directions in the tracked-axis complement.
% Numerically this is the same identifier as Q here; both labels are kept
% because Q is the control-aware reduction branch and V is the later
% process/smooth-noise lineage. Equality is an expected audit result.
models.V=models.Q;
models.V.description='V: control-aware iterative IVR process lineage';

% X: preserve V loading P but use a mildly oblique covariance-weighted
% dual extractor. Held-out screening in the original copy selected alpha=.02.
[A,B,P,R,S,st]=control_aware_oblique_ivr_varx(y,u,cfg.ell,cfg.tracked,cfg.oblique_alpha);
st=complete_stats(st,y,u,P,R,E);
models.X=pack(A,B,P,R,S,st,'X: control-aware mild oblique dual, alpha=0.02');

% Y: explicit no-whitening/direct-update label. The audited X identifier is
% already a raw-centered direct-update implementation, so X/Y are expected
% to be numerically identical under this fixture; keep both labels rather
% than inventing an artificial difference.
[A,B,P,R,S,st]=control_aware_direct_update_varx(y,u,cfg.ell,cfg.tracked,cfg.oblique_alpha);
st=complete_stats(st,y,u,P,R,E);
models.Y=pack(A,B,P,R,S,st,'Y: no-whitening direct eigenspace update');

% Z: same strict rank-full s=1 whitened ablation as R. R is retained for
% the Mo-Qin baseline branch; Z records the later copyX-plant whitening
% ablation lineage. Equality is expected and documented.
models.Z=models.R;
models.Z.description='Z: strict whitened rank-full s=1 ablation';

for n={'U','V','X','Y','Z'}
 models.(n{1}).online_reidentify=false; models.(n{1}).online_covariance_update=true;
 models.(n{1}).reidentify_period=inf; models.(n{1}).reidentify_recent=0;
 models.(n{1}).reidentify_mode='fixed'; models.(n{1}).covariance_mode='rolling';
 models.(n{1}).bias_correction=false; models.(n{1}).bias_gamma=0;
 models.(n{1}).dinkla_window=cfg.noise_window;
end
end
function E=eye_cols(p,t), E=zeros(p,numel(t)); E(t,:)=eye(numel(t)); end
function st=base_stats(y,u,P,R,E)
yc=y-mean(y,2); st.y_mean=mean(y,2); st.u_mean=mean(u,2); st.tracked_projection_error=norm(P*R'*E-E,'fro'); st.reconstruction_residual=norm(yc-P*(R'*yc),'fro')/max(norm(yc,'fro'),eps); st.dual_error=norm(R'*P-eye(size(P,2)),'fro');
end
function st=complete_stats(st,y,u,P,R,E)
base=base_stats(y,u,P,R,E); f=fieldnames(base); for i=1:numel(f), st.(f{i})=base.(f{i}); end
end
function c=pack(A,B,P,R,S,st,desc)
c.A=A; c.B=B; c.P=P; c.R=R; c.Sigma_eps=(S+S')/2; c.y_mean=st.y_mean; c.u_mean=st.u_mean; c.stats=st; c.description=desc;
end
