function models=copyOPQSTR_build_models(D,cfg)
% Restore Obsidian O/P/Q/R/S/T distinctions on one common dataset.
y=D.y_off; u=D.u_off; E=eye_cols(cfg.p,cfg.tracked);
[A,B,P,~,R,~,~,S]=predvarx_identify_oblique(y,u,cfg.ell,1,1,cfg.A,cfg.B,D.C,cfg.n,cfg.m,cfg.p);
st=base_stats(y,u,P,R,E); models.O=pack(A,B,P,R,S,st,'O: free oblique IVR');
[A,B,P,R,S,st]=control_ready_subspace_varx(y,u,cfg.ell); st=complete_stats(st,y,u,P,R,E); models.P=pack(A,B,P,R,S,st,'P: centered PCA/SVD orthogonal');
[A,B,P,R,S,st]=control_aware_iterative_ivr_varx(y,u,cfg.ell,cfg.tracked); st=complete_stats(st,y,u,P,R,E); models.Q=pack(A,B,P,R,S,st,'Q: tracked-axis control-aware IVR');
[A,B,P,~,R,~,~,S,~,~,~,~,info]=predvarx_identify_strict_whitened(y,u,cfg.ell); st=base_stats(y,u,P,R,E); st.whitening_applied=info.whitening_applied; st.complete_algorithm1=info.complete_algorithm1; models.R=pack(A,B,P,R,S,st,'R: rank-full s=1 whitened IVR');
models.S=models.Q; models.S.description='S: periodic online re-identification of Q geometry'; models.S.online_reidentify=true; models.S.online_covariance_update=true; models.S.reidentify_period=cfg.reidentify_period;
models.T=models.Q; models.T.description='T: fixed Q geometry with online covariance updates only'; models.T.online_reidentify=false; models.T.online_covariance_update=true; models.T.reidentify_period=inf;
for n={'O','P','Q','R'}, models.(n{1}).online_reidentify=false; models.(n{1}).online_covariance_update=true; models.(n{1}).reidentify_period=inf; end
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
