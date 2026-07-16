%% Structural regression test for the Obsidian-restored OPQSTR suite
here=fileparts(fileparts(mfilename('fullpath'))); addpath(here);
cfg=copyOPQSTR_common_config();
assert(cfg.p==30 && cfg.ell==5 && isequal(cfg.tracked,[1 2]));
assert(cfg.T_off==1500 && cfg.T_cl==1200 && cfg.N==18);
assert(cfg.seed==20260710 && cfg.reidentify_period==30);
D=copyOPQSTR_generate_common_data(cfg);
models=copyOPQSTR_build_models(D,cfg);
required={'O','P','Q','R','S','T'};
for i=1:numel(required), assert(isfield(models,required{i})); end
E=zeros(cfg.p,numel(cfg.tracked)); E(cfg.tracked,:)=eye(numel(cfg.tracked));
assert(norm(models.O.R'*models.O.P-eye(cfg.ell),'fro')<1e-8);
assert(norm(models.P.R-models.P.P,'fro')<1e-12);
assert(norm(models.Q.P*models.Q.R'*E-E,'fro')<1e-10);
assert(norm(models.R.R'*models.R.P-eye(cfg.ell),'fro')<1e-8);
assert(models.R.stats.whitening_applied==true);
assert(models.S.online_reidentify==true && models.S.reidentify_period==cfg.reidentify_period);
assert(models.T.online_reidentify==false && models.T.online_covariance_update==true);
fprintf('PASS copyOPQSTR unified structural test.\n');
