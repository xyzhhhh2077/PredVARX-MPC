%% Direct MAT audit: copyV orthogonal vs copyX oblique
clear; clc; here=fileparts(mfilename('fullpath'));
repo=fileparts(fileparts(here));
v=load(fullfile(repo,'experiments','copyV_iterative_ivr','results','copyV_iterative_ivr_data.mat'));
x=load(fullfile(here,'results','copyX_control_aware_oblique_data.mat'));

% Fairness assertions from persisted artifacts.
assert(isequal(v.A_true,x.A_true) && isequal(v.B_true,x.B_true) && isequal(v.C_true,x.C_true));
assert(isequal(v.Rf,x.Rf));
assert(isequal(v.sigma_w_profile,x.sigma_w_profile));
assert(isequal(v.sigma_e_profile,x.sigma_e_profile));
assert(norm(v.Phat-x.Phat,'fro')<1e-12,'P differs; comparison is not isolated to R');

warm=151:v.T_cl; seg_len=floor(v.T_cl/5); tracked=v.tracked;
errV=v.y(tracked,:)-v.Rf(tracked,:); errX=x.y(tracked,:)-x.Rf(tracked,:);
first_y=find(vecnorm(v.y-x.y,2,1)>1e-12,1,'first');
first_u=find(vecnorm(v.u-x.u,2,1)>1e-12,1,'first');
first_pred=find(vecnorm(v.yhat-x.yhat,2,1)>1e-12,1,'first');

fprintf('MAT FAIRNESS PASS: plant/reference/noise envelopes identical; ||P_V-P_X||=%.3e\n',norm(v.Phat-x.Phat,'fro'));
fprintf('||R_V-R_X||=%.6f, dualX=%.3e, asymX=%.6f\n',norm(v.Rhat-x.Rhat,'fro'),norm(x.Rhat'*x.Phat-eye(x.ell),'fro'),norm(x.Phat*x.Rhat'-(x.Phat*x.Rhat')','fro'));
fprintf('first differing prediction=%d, control=%d, plant output=%d\n',first_pred,first_u,first_y);
fprintf('trajectory norms: ||yV-yX||=%.6f, ||uV-uX||=%.6f, ||ypredV-ypredX||=%.6f\n',norm(v.y-x.y,'fro'),norm(v.u-x.u,'fro'),norm(v.yhat-x.yhat,'fro'));
fprintf('max abs differences: y=%.6f, u=%.6f, ypred=%.6f, cost=%.6f, qpResidual=%.3e\n', ...
 max(abs(v.y(:)-x.y(:))),max(abs(v.u(:)-x.u(:))),max(abs(v.yhat(:)-x.yhat(:))), ...
 max(abs(v.costJ(:)-x.costJ(:))),max(abs(v.max_cc_violation(:)-x.max_cc_violation(:))));

fprintf('\nPer-segment MAE (V then X; delta%%):\n');
seg=zeros(5,6);
for s=1:5
    ix=max((s-1)*seg_len+1,151):min(s*seg_len,v.T_cl);
    mv=mean(abs(errV(:,ix)),2); mx=mean(abs(errX(:,ix)),2);
    seg(s,:)=[mv(1) mx(1) 100*(mx(1)-mv(1))/mv(1), mv(2) mx(2) 100*(mx(2)-mv(2))/mv(2)];
    fprintf('seg%d y1 %.5f -> %.5f (%+.2f%%), y2 %.5f -> %.5f (%+.2f%%)\n',s,seg(s,:));
end

% Focused visual comparison generated only from persisted MAT arrays.
fig=figure('Position',[50 50 2200 1400],'Color','w'); t=1:v.T_cl;
tl=tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
ax1=nexttile(tl); plot(t,v.Rf(1,:),'k--','LineWidth',1.2); hold on; plot(t,v.y(1,:),'b','LineWidth',.75); plot(t,x.y(1,:),'r','LineWidth',.75); yline(v.y_max,'m--'); grid on; title('Tracked output y_1 from persisted MAT'); legend('reference','copyV R=P','copyX P not equal R','limit','Location','eastoutside');
ax2=nexttile(tl); plot(t,v.Rf(2,:),'k--','LineWidth',1.2); hold on; plot(t,v.y(2,:),'b','LineWidth',.75); plot(t,x.y(2,:),'r','LineWidth',.75); yline(v.y_max,'m--'); grid on; title('Tracked output y_2 from persisted MAT'); legend('reference','copyV R=P','copyX P not equal R','limit','Location','eastoutside');
ax3=nexttile(tl); plot(t,vecnorm(v.u-x.u,2,1),'Color',[.1 .5 .1]); grid on; ylabel('norm(u_V-u_X)'); title(sprintf('Control divergence: Frobenius norm %.3f; max step %.3f',norm(v.u-x.u,'fro'),max(vecnorm(v.u-x.u,2,1))));
ax4=nexttile(tl); plot(t,v.costJ,'b'); hold on; plot(t,x.costJ,'r'); grid on; ylabel('J'); xlabel('time'); title(sprintf('Full MPC cost: mean warm %.2f (V) vs %.2f (X)',mean(v.costJ(warm),'omitnan'),mean(x.costJ(warm),'omitnan'))); legend('copyV','copyX','Location','eastoutside');
linkaxes([ax1 ax2 ax3 ax4],'x'); title(tl,'MAT artifact audit: identical P, different R');
outpng=fullfile(here,'results','copyV_copyX_mat_comparison.png'); print(fig,outpng,'-dpng','-r180');

save(fullfile(here,'results','copyV_copyX_mat_audit.mat'),'seg','first_y','first_u','first_pred');
fprintf('saved %s\n',outpng);
