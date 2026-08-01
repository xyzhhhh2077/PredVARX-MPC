function run_cdr_au_closed_loop()
%RUN_CDR_AU_CLOSED_LOOP Formal 1200-step AU-SMPC experiment on true CDR.
here=fileparts(mfilename('fullpath')); results_dir=fullfile(here,'results');
if ~exist(results_dir,'dir'), mkdir(results_dir); end
cfg=build_cdr_closed_loop_config();
au=simulate_cdr_closed_loop(cfg);
baseline=simulate_cdr_zero_input_baseline(cfg,au.reference);
improvement=1-au.RMSE./baseline.RMSE;

Gau=cfg.H*cfg.model_control.P*((eye(cfg.ell)-cfg.model_control.A)\cfg.model_control.B);
SyInv=diag(1./cfg.y_scale); Su=diag(cfg.u_scale);
Gtrue=cfg.Etask'*SyInv*cfg.C_plant*((eye(size(cfg.A_plant))-cfg.A_plant)\cfg.B_plant)*Su;
gain_relative_error=norm(Gau-Gtrue,'fro')/norm(Gtrue,'fro');
physical_action_rms=sqrt(mean(au.u.^2,'all'));
input_violation_rate=mean((au.u<cfg.u_min-1e-9)|(au.u>cfg.u_max+1e-9),'all');

save(fullfile(results_dir,'copyAX_cdr_au_closed_loop_data.mat'), ...
    'cfg','au','baseline','improvement','Gau','Gtrue','gain_relative_error', ...
    'physical_action_rms','input_violation_rate');

metrics=fullfile(results_dir,'copyAX_cdr_au_closed_loop_metrics.txt');
fid=fopen(metrics,'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'experiment copyAX CDR AU closed loop\n');
fprintf(fid,'closed_loop_mpc_run 1\ncontroller AU centered_smpc_step\n');
fprintf(fid,'T %d\nwarmup %d\nN %d\nQ_task_scale 80\nRu_scale 0.18\n',cfg.T,cfg.warmup,cfg.N);
fprintf(fid,'reference_source identified_AU_model_only\ntrue_plant_used_for controller_design 0 reference_design 0 simulation 1 evaluation 1\n');
fprintf(fid,'reference_amplitude %.12g\ntask_limit %.12g\nmax_chance_tightening %.12g\n',cfg.reference_amplitude,cfg.task_limit,cfg.max_chance_tightening);
fprintf(fid,'AU_RMSE %s\nzero_input_RMSE %s\nRMSE_improvement_fraction %s\n',mat2str(au.RMSE',10),mat2str(baseline.RMSE',10),mat2str(improvement',10));
fprintf(fid,'AU_MAE %s\nzero_input_MAE %s\nAU_Bias %s\n',mat2str(au.MAE',10),mat2str(baseline.MAE',10),mat2str(au.Bias',10));
fprintf(fid,'QP_success %.12f\nfallback %d\nmax_qp_row %.12e\n',au.qp_success,au.fallback,au.max_qp);
fprintf(fid,'physical_action_rms %.12g\ninput_saturation_rate %.12g\ninput_violation_rate %.12g\n',physical_action_rms,au.input_saturation_rate,input_violation_rate);
fprintf(fid,'task_violation_rate %s\n',mat2str(au.task_violation_rate',10));
fprintf(fid,'AU_gain_singular_values %s\ntrue_gain_singular_values %s\ngain_relative_error %.12g\n',mat2str(svd(Gau)',10),mat2str(svd(Gtrue)',10),gain_relative_error);
fprintf(fid,'interpretation closed_loop_executed_and_both_task_RMSEs_improved_vs_same-noise_zero-input_baseline\n');
fprintf(fid,'boundary improvement_is_empirical_not_a_stability_or_recursive-feasibility_proof; identified_gain_mismatch_remains\n');
clear cleanup

fig=figure('Visible','off','Color','w','Position',[80 80 1220 820]);
tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
t=1:cfg.T;
for j=1:2
    ax=nexttile(tl,j); hold(ax,'on');
    plot(ax,t,au.reference(j,:),'k--','LineWidth',1.5,'DisplayName','reference');
    plot(ax,t,au.s(j,:),'Color',[0.00 0.42 0.68],'LineWidth',1.0,'DisplayName','AU-SMPC');
    plot(ax,t,baseline.s(j,:),'Color',[0.72 0.26 0.17],'LineWidth',0.9,'DisplayName','zero input');
    yline(ax,cfg.task_limit,':','Color',[0.35 0.35 0.35],'HandleVisibility','off');
    yline(ax,-cfg.task_limit,':','Color',[0.35 0.35 0.35],'HandleVisibility','off');
    grid(ax,'on'); xlabel(ax,'time step'); ylabel(ax,sprintf('task axis %d',j));
    title(ax,sprintf('Axis %d: RMSE improvement %.1f%%',j,100*improvement(j)));
    if j==1, legend(ax,'Location','best'); end
end
ax=nexttile(tl,3,[1 2]); plot(ax,t,au.u','LineWidth',0.75); grid(ax,'on');
input_span=max(abs(au.u),[],'all'); ylim(ax,1.15*[-input_span input_span]);
yline(ax,0,'k:','HandleVisibility','off'); xlabel(ax,'time step'); ylabel(ax,'physical CDR input');
title(ax,sprintf('Eight applied inputs (plant bounds [%.0f, %.0f], saturation %.1f%%)',cfg.u_min,cfg.u_max,100*au.input_saturation_rate));
ax=nexttile(tl,5); plot(ax,t,au.maxcc,'k','LineWidth',0.9); yline(ax,0,'r--'); grid(ax,'on'); xlabel(ax,'time step'); ylabel(ax,'max(AU-b)'); title(ax,sprintf('QP success %.3f, fallback %d',au.qp_success,au.fallback));
ax=nexttile(tl,6); bar(ax,100*[improvement(:), 1-au.MAE./baseline.MAE]); grid(ax,'on'); yline(ax,0,'k:'); xlabel(ax,'task axis'); ylabel(ax,'improvement (%)'); legend(ax,{'RMSE','MAE'},'Location','best'); title(ax,'Against same-noise zero-input baseline');
sgtitle(tl,'copyAX: AU-based SMPC control of ControlGym CDR');
exportgraphics(fig,fullfile(results_dir,'copyAX_cdr_au_closed_loop_fig.png'),'Resolution',180); close(fig);
fprintf('copyAX CDR AU closed loop: RMSE improvement=%s, QP=%.3f, fallback=%d\n',mat2str(improvement',4),au.qp_success,au.fallback);
end
