%% 副本 Z 四行图: C + Z1 + Z2
load('copyZ_data.mat');

% 实际噪声参考线
sw0_v=0.1; nf_v=nf; ni_v=ni;
act_sw=zeros(1,T_cl);
for seg=1:10, k1=(seg-1)*1000+1; k2=min(seg*1000,T_cl);
    sg=min(floor((k1-1)/ni_v)+1,length(nf_v));
    act_sw(k1:k2)=sw0_v*nf_v(sg);
end

figure('Position',[50,50,2000,1600],'Name','副本Z');

% 行1: y跟踪
subplot(5,1,1); hold on; grid on;
plot(1:T_cl,Rf(2,:),'k:','LineWidth',1.2);
plot(1:T_cl,yC(2,:),'Color',[0,0.4,1],'LineWidth',1);
plot(1:T_cl,yZ1(2,:),'Color',[1,0.4,0],'LineWidth',1);
plot(1:T_cl,yZ2(2,:),'Color',[0,0.7,0.4],'LineWidth',1);
yline(y_max,'g--','LineWidth',1.2);
for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
legend('ref',sprintf('C(%.3f,%d)',eC,violC),sprintf('Z1(%.3f,%d)',eZ1,violZ1),sprintf('Z2(%.3f,%d)',eZ2,violZ2),'y_{max}','Location','best');
ylabel('y_2'); title(sprintf('y tracking  Rw=%d  y_{max}=%.2f',2,y_max));

% 行2: 噪声eps
subplot(5,1,2); hold on; grid on;
plot(500:T_cl,sC_eps(500:end),'Color',[0,0.4,1],'LineWidth',1);
plot(500:T_cl,sZ1_eps(500:end),'Color',[1,0.4,0],'LineWidth',1);
plot(500:T_cl,sZ2_eps(500:end),'Color',[0,0.7,0.4],'LineWidth',1);
plot(500:T_cl,act_sw(500:end),'k-','LineWidth',1);
for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
legend('C','Z1','Z2','true sw','Location','best');
ylabel('sigma_eps');

% 行3: 静态噪声 ebar
subplot(5,1,3); hold on; grid on;
plot(500:T_cl,sC_ebar(500:end),'Color',[0,0.4,1],'LineWidth',1);
plot(500:T_cl,sZ1_ebar(500:end),'Color',[1,0.4,0],'LineWidth',1);
plot(500:T_cl,sZ2_ebar(500:end),'Color',[0,0.7,0.4],'LineWidth',1);
plot(500:T_cl,zeros(1,T_cl-499)+se0,'k-','LineWidth',0.8);
for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
legend('C','Z1','Z2','true se','Location','best');
ylabel('sigma_ebar');

% 行4: 控制 u
subplot(5,1,4); hold on; grid on;
plot(1:T_cl,uC(1,:),'Color',[0,0.4,1]); plot(1:T_cl,uC(2,:),'Color',[0.8,0,0.7]); plot(1:T_cl,uC(3,:),'Color',[0,0.7,0.4]);
for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
ylabel('u'); title('C: u_1 u_2 u_3'); ylim([-6 6]);

% 行4: 代价 J (log)
subplot(5,1,5); hold on; grid on;
semilogy(500:T_cl,max(costC(500:end),1e-10),'Color',[0,0.4,1],'LineWidth',1);
semilogy(500:T_cl,max(costZ1(500:end),1e-10),'Color',[1,0.4,0],'LineWidth',1);
semilogy(500:T_cl,max(costZ2(500:end),1e-10),'Color',[0,0.7,0.4],'LineWidth',1);
for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
legend('C','Z1','Z2','Location','best'); ylabel('J_k (log)'); xlabel('k');

sgtitle(sprintf('副本Z: C(%.3f,%d) Z1(%.3f,%d) Z2(%.3f,%d)  Rw=%d',eC,violC,eZ1,violZ1,eZ2,violZ2,2));
saveas(gcf,'copyZ_4row.png');
fprintf('图: copyZ_4row.png\n');