%% Structural test for full H/K/O/P/Q/R/S/T/U/V/X/Y/Z suite
here=fileparts(fileparts(mfilename('fullpath'))); addpath(here);
cfg=copyALL_common_config();
assert(cfg.seed==20260710 && cfg.n==6 && cfg.m==3 && cfg.p==30 && cfg.ell==5);
assert(cfg.T_off==1500 && cfg.T_cl==1200 && cfg.N==18 && isequal(cfg.tracked,[1 2]));
D=copyALL_generate_common_data(cfg); models=copyALL_build_models(D,cfg);
required={'H','K','O','P','Q','R','S','T','U','V','X','Y','Z'};
for i=1:numel(required), assert(isfield(models,required{i}),required{i}); end
E=zeros(cfg.p,numel(cfg.tracked)); E(cfg.tracked,:)=eye(numel(cfg.tracked));
assert(norm(models.Q.P*models.Q.R'*E-E,'fro')<1e-10);
assert(norm(models.U.P*models.U.R'*E-E,'fro')<1e-10);
assert(norm(models.V.P*models.V.R'*E-E,'fro')<1e-10);
assert(norm(models.X.P*models.X.R'*E-E,'fro')<1e-10);
assert(norm(models.Y.P*models.Y.R'*E-E,'fro')<1e-10);
assert(norm(models.X.R-models.X.P,'fro')>1e-6);
assert(norm(models.Y.R-models.Y.P,'fro')>1e-6);
assert(norm(models.R.R'*models.R.P-eye(cfg.ell),'fro')<1e-8);
assert(norm(models.Z.R'*models.Z.P-eye(cfg.ell),'fro')<1e-8);
assert(norm(models.Q.P-models.V.P,'fro')<1e-12);
assert(norm(models.X.P-models.Y.P,'fro')<1e-12);
assert(norm(models.X.R-models.Y.R,'fro')<1e-12);
assert(norm(models.R.P-models.Z.P,'fro')<1e-12);
assert(models.H.online_reidentify && strcmp(models.H.reidentify_mode,'sliding'));
assert(models.K.online_reidentify && strcmp(models.K.reidentify_mode,'cumulative'));
assert(models.S.online_reidentify && strcmp(models.S.reidentify_mode,'offline_plus_recent'));
assert(~models.T.online_reidentify);
fprintf('PASS copyALL unified structural test for 13 versions.\n');