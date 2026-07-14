%% Clear visual comparison from persisted MAT artifacts
clear; clc; close all; here=fileparts(mfilename('fullpath')); repo=fileparts(fileparts(here));
v=load(fullfile(repo,'experiments','copyV_iterative_ivr','results','copyV_iterative_ivr_data.mat'));
x=load(fullfile(here,'results','copyX_control_aware_oblique_data.mat'));
assert(norm(v.Phat-x.Phat,'fro')<1e-12); assert(isequal(v.Rf,x.Rf));
t=1:v.T_cl; warm=151:v.T_cl; tracked=v.tracked; win=25;
errV=abs(v.y(tracked,:)-v.Rf(tracked,:)); errX=abs(x.y(tracked,:)-x.Rf(tracked,:));
rollV=movmean(errV,win,2); rollX=movmean(errX,win,2);

fig=figure('Position',[30 30 2600 1700],'Color','w');
tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');

ax1=nexttile(tl,1); plot(t,v.Rf(1,:),'k--','LineWidth',1.2); hold on; plot(t,v.y(1,:),'b','LineWidth',.8); plot(t,v.Rf(2,:),'Color',[.35 .35 .35],'LineStyle','--','LineWidth',1.2); plot(t,v.y(2,:),'Color',[0 .55 .75],'LineWidth',.8); yline(v.y_max,'m--','LineWidth',1.2); grid on; title('copyV orthogonal extractor: R = P'); ylabel('outputs'); legend('r_1','y_1','r_2','y_2','limit','Location','eastoutside');
ax2=nexttile(tl,2); plot(t,x.Rf(1,:),'k--','LineWidth',1.2); hold on; plot(t,x.y(1,:),'r','LineWidth',.8); plot(t,x.Rf(2,:),'Color',[.35 .35 .35],'LineStyle','--','LineWidth',1.2); plot(t,x.y(2,:),'Color',[.9 .45 0],'LineWidth',.8); yline(x.y_max,'m--','LineWidth',1.2); grid on; title(sprintf('copyX oblique extractor: P not equal R, alpha = %.2f',x.oblique_alpha)); ylabel('outputs'); legend('r_1','y_1','r_2','y_2','limit','Location','eastoutside');
yv_tracked=v.y(tracked,:); yx_tracked=x.y(tracked,:); rf_tracked=v.Rf(tracked,:);
all_output_values=[yv_tracked(:); yx_tracked(:); rf_tracked(:)];
ylims=[min(all_output_values)-.15, max([max(all_output_values),v.y_max])+.15]; ylim(ax1,ylims); ylim(ax2,ylims);

ax3=nexttile(tl,3); plot(t,x.y(1,:)-v.y(1,:),'r','LineWidth',.85); hold on; plot(t,x.y(2,:)-v.y(2,:),'Color',[.9 .45 0],'LineWidth',.85); yline(0,'k--'); grid on; ylabel('copyX - copyV'); title(sprintf('Direct output difference: max |dy| = %.4f',max(abs(x.y(:)-v.y(:))))); legend('dy_1','dy_2','Location','eastoutside');
ax4=nexttile(tl,4); plot(t,rollV(1,:),'b','LineWidth',1.0); hold on; plot(t,rollX(1,:),'r--','LineWidth',1.0); plot(t,rollV(2,:),'Color',[0 .55 .75],'LineWidth',1.0); plot(t,rollX(2,:),'Color',[.9 .45 0],'LineStyle','--','LineWidth',1.0); grid on; ylabel(sprintf('%d-step mean |error|',win)); title('Rolling tracking error makes the small trade-off visible'); legend('V y_1','X y_1','V y_2','X y_2','Location','eastoutside');

ax5=nexttile(tl,5); du=vecnorm(x.u-v.u,2,1); plot(t,du,'Color',[.15 .55 .2],'LineWidth',.85); grid on; ylabel('norm(u_X-u_V)'); xlabel('time step'); title(sprintf('Control difference: max %.4f, trajectory norm %.4f',max(du),norm(x.u-v.u,'fro')));

ax6=nexttile(tl,6); vals=[v.MAE(:),x.MAE(:); mean(v.costJ(warm),'omitnan')/100,mean(x.costJ(warm),'omitnan')/100; v.qp_success_rate,x.qp_success_rate];
% Rows: MAE y1, MAE y2, cost/100, QP success
b=bar(ax6,vals,'grouped'); b(1).FaceColor=[.1 .35 .85]; b(2).FaceColor=[.85 .2 .15]; grid(ax6,'on'); xticklabels(ax6,{'MAE y_1','MAE y_2','mean J / 100','QP success'}); ylabel(ax6,'scaled metric'); title(ax6,'Metric comparison (cost divided by 100 only for display)'); legend(ax6,'copyV R=P','copyX P not equal R','Location','eastoutside');

linkaxes([ax1 ax2 ax3 ax4 ax5],'x'); xlim([1 v.T_cl]);
title(tl,sprintf('Persisted-MAT comparison: same plant, same P, different R (||R_V-R_X||_F = %.4f)',norm(v.Rhat-x.Rhat,'fro')),'FontWeight','bold');
out=fullfile(here,'results','copyV_copyX_clear_comparison.png'); print(fig,out,'-dpng','-r180'); fprintf('saved %s\n',out);
