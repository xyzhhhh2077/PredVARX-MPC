%% copyOPQSTR_unified — restore Obsidian O/P/Q/S/T/R differences on common parameters
clear; clc; close all; here=fileparts(mfilename('fullpath')); addpath(here);
cfg=copyOPQSTR_common_config(); D=copyOPQSTR_generate_common_data(cfg); models=copyOPQSTR_build_models(D,cfg);
names={'O','P','Q','R','S','T'}; results=cell(size(names));
for i=1:numel(names)
 n=names{i}; results{i}=run_copyOPQSTR_case(n,models.(n),D,cfg);
 r=results{i}; fprintf('%s: MAE=[%.4f %.4f] RMSE=[%.4f %.4f] pred=%.4f cover=%.3e dual=%.3e QP=%.3f fallback=%d reid=%d\n',n,r.MAE,r.RMSE,r.prediction_rmse,r.coverage_error,r.dual_error,r.qp_success_rate,r.fallbacks,r.reidentify_count);
end
results_dir=fullfile(here,'results'); if ~exist(results_dir,'dir'), mkdir(results_dir); end
fid=fopen(fullfile(results_dir,'copyOPQSTR_unified_metrics.csv'),'w');
fprintf(fid,'method,mae_y1,mae_y2,rmse_y1,rmse_y2,prediction_rmse,reconstruction_residual,coverage_error,dual_error,avg_cost,qp_success,fallbacks,reidentify_count,upper_viol_y1,upper_viol_y2,abs_viol_y1,abs_viol_y2,u_rms_1,u_rms_2,u_rms_3,max_constraint\n');
for i=1:numel(names), r=results{i}; fprintf(fid,'%s,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n',r.method,r.MAE,r.RMSE,r.prediction_rmse,r.reconstruction_residual,r.coverage_error,r.dual_error,r.avg_cost,r.qp_success_rate,r.fallbacks,r.reidentify_count,r.upper_violation_rate,r.abs_violation_rate,r.u_rms,r.max_constraint); end
fclose(fid);
fig=figure('Position',[80 80 1900 1100],'Color','w'); tl=tiledlayout(3,1,'TileSpacing','compact'); t=1:cfg.T_cl;
for q=1:2, ax=nexttile; plot(t,D.Rf(q,:),'k--','LineWidth',1.3); hold on; for i=1:numel(names), plot(t,results{i}.y(q,:),'LineWidth',0.65); end; grid on; ylabel(sprintf('y_%d',q)); if q==1, legend([{'reference'},names],'Location','eastoutside'); end; end
ax=nexttile; vals=zeros(numel(names),4); for i=1:numel(names), vals(i,:)=[mean(results{i}.MAE),results{i}.prediction_rmse,results{i}.qp_success_rate,results{i}.coverage_error]; end; bar(vals); set(gca,'XTickLabel',names); legend({'mean MAE','prediction RMSE','QP success','coverage error'}); grid on; title('Unified OPQSTR metrics'); title(tl,'Obsidian-restored OPQSTR on one common copyX parameter set');
print(fig,fullfile(results_dir,'copyOPQSTR_unified_fig'),'-dpng','-r160');
assert(results{5}.reidentify_count==floor(cfg.T_cl/cfg.reidentify_period)); assert(results{6}.reidentify_count==0);
assert(all(cellfun(@(r) all(isfinite(r.MAE)),results)));
assert(all(cellfun(@(r) r.max_constraint<=1e-7 || r.qp_success_rate<1,results)));
fprintf('Saved %s\n',fullfile(results_dir,'copyOPQSTR_unified_metrics.csv'));
