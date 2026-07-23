%% copyALL_unified — full H/K/O/P/Q/R/S/T/U/V/X/Y/Z comparison
clear; clc; close all; here=fileparts(mfilename('fullpath')); addpath(here);
cfg=copyALL_common_config();
D=copyALL_generate_common_data(cfg);
models=copyALL_build_models(D,cfg);
names={'H','K','O','P','Q','R','S','T','U','V','X','Y','Z'};
results=cell(size(names));
for i=1:numel(names)
 n=names{i};
 results{i}=run_copyALL_case(n,models.(n),D,cfg);
 r=results{i};
 fprintf('%s: MAE=[%.4f %.4f] RMSE=[%.4f %.4f] pred=%.4f cover=%.3e dual=%.3e QP=%.3f fallback=%d reid=%d\n', ...
  n,r.MAE,r.RMSE,r.prediction_rmse,r.coverage_error,r.dual_error,r.qp_success_rate,r.fallbacks,r.reidentify_count);
end

results_dir=fullfile(here,'results'); if ~exist(results_dir,'dir'), mkdir(results_dir); end
csv_path=fullfile(results_dir,'copyALL_unified_metrics.csv');
fid=fopen(csv_path,'w');
fprintf(fid,'method,mae_y1,mae_y2,rmse_y1,rmse_y2,prediction_rmse,reconstruction_residual,coverage_error,dual_error,avg_cost,qp_success,fallbacks,reidentify_count,upper_viol_y1,upper_viol_y2,abs_viol_y1,abs_viol_y2,u_rms_1,u_rms_2,u_rms_3,max_constraint\n');
for i=1:numel(names)
 r=results{i};
 fprintf(fid,'%s,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%d,%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n', ...
  r.method,r.MAE,r.RMSE,r.prediction_rmse,r.reconstruction_residual,r.coverage_error,r.dual_error,r.avg_cost,r.qp_success_rate,r.fallbacks,r.reidentify_count,r.upper_violation_rate,r.abs_violation_rate,r.u_rms,r.max_constraint);
end
fclose(fid);

fig=figure('Position',[60 40 2300 1350],'Color','w');
tl=tiledlayout(4,1,'TileSpacing','compact','Padding','compact'); t=1:cfg.T_cl;
colors=lines(numel(names));
for q=1:2
 ax=nexttile; plot(t,D.Rf(q,:),'k--','LineWidth',1.5,'DisplayName','reference'); hold on;
 for i=1:numel(names), plot(t,results{i}.y(q,:),'Color',colors(i,:),'LineWidth',0.60,'DisplayName',names{i}); end
 grid on; ylabel(sprintf('y_%d',q)); title(sprintf('Tracked output y_%d',q));
 if q==1, legend('Location','eastoutside','NumColumns',2); end
end
ax=nexttile; vals=zeros(numel(names),3);
for i=1:numel(names), vals(i,:)=[mean(results{i}.MAE),results{i}.prediction_rmse,results{i}.qp_success_rate]; end
bar(vals); grid on; set(gca,'XTick',1:numel(names),'XTickLabel',names);
legend({'mean MAE','prediction RMSE','QP success'},'Location','eastoutside'); title('Tracking, prediction, and feasibility');
ax=nexttile; vals2=zeros(numel(names),3);
for i=1:numel(names), vals2(i,:)=[results{i}.coverage_error,results{i}.reconstruction_residual,results{i}.fallbacks/cfg.T_cl]; end
bar(vals2); grid on; set(gca,'XTick',1:numel(names),'XTickLabel',names);
legend({'coverage error','reconstruction residual','fallback rate'},'Location','eastoutside'); title('Geometry and fallback diagnostics'); xlabel('copy version');
title(tl,'Full unified PredVARX-SMPC copy comparison: H K O P Q R S T U V X Y Z');
png_path=fullfile(results_dir,'copyALL_unified_fig.png'); print(fig,png_path,'-dpng','-r170');

% One canonical figure per runnable copy.  Each image is generated from the
% same result object used by the unified CSV, so it cannot silently drift to
% another plant, seed, reference, noise realization, or SMPC convention.
individual_dir=fullfile(results_dir,'individual');
if ~exist(individual_dir,'dir'), mkdir(individual_dir); end
for i=1:numel(names)
 n=names{i}; r=results{i};
 f=figure('Position',[80 50 1900 1200],'Color','w','Visible','off');
 tli=tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
 for q=1:2
  nexttile; plot(t,D.Rf(q,:),'k--','LineWidth',1.45,'DisplayName','reference'); hold on;
  plot(t,r.y(q,:),'Color',colors(i,:),'LineWidth',0.95,'DisplayName',['copy' n]);
  yline(cfg.y_max,'r:','LineWidth',0.85,'DisplayName','constraint');
  yline(-cfg.y_max,'r:','LineWidth',0.85,'HandleVisibility','off');
  grid on; ylabel(sprintf('y_%d',q)); title(sprintf('copy%s tracked output y_%d',n,q));
  if q==1, legend('Location','eastoutside'); end
 end
 nexttile;
 metric_values=[mean(r.MAE),r.prediction_rmse,r.qp_success_rate,r.fallbacks/cfg.T_cl,r.coverage_error,r.reconstruction_residual];
 bar(metric_values,'FaceColor',colors(i,:)); grid on;
 set(gca,'XTick',1:6,'XTickLabel',{'mean MAE','pred RMSE','QP success','fallback rate','coverage','reconstruction'},'XTickLabelRotation',15);
 title(sprintf('copy%s diagnostics: MAE=[%.4f %.4f], QP=%.2f%%, fallback=%d, reid=%d', ...
  n,r.MAE(1),r.MAE(2),100*r.qp_success_rate,r.fallbacks,r.reidentify_count));
 title(tli,sprintf('copy%s — unified PredVARX-SMPC run',n));
 individual_path=fullfile(individual_dir,sprintf('copy%s_unified.png',n));
 print(f,individual_path,'-dpng','-r170'); close(f);
end

idx=@(name)find(strcmp(names,name),1);
assert(results{idx('H')}.reidentify_count==floor(cfg.T_cl/cfg.hk_reidentify_period));
assert(results{idx('K')}.reidentify_count==floor(cfg.T_cl/cfg.hk_reidentify_period));
assert(results{idx('S')}.reidentify_count==floor(cfg.T_cl/cfg.reidentify_period));
assert(results{idx('T')}.reidentify_count==0);
assert(norm(results{idx('Q')}.MAE-results{idx('V')}.MAE)<1e-12);
assert(norm(results{idx('X')}.MAE-results{idx('Y')}.MAE)<1e-12);
assert(norm(results{idx('R')}.MAE-results{idx('Z')}.MAE)<1e-12);
assert(all(cellfun(@(r)all(isfinite(r.MAE)),results)));
fprintf('Saved %s\n',csv_path);
fprintf('Saved %s\n',png_path);
fprintf('Saved %d individual copy figures in %s\n',numel(names),individual_dir);