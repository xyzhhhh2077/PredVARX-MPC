%% 副本 S_nosoft — 滑动窗口重辨识 + Pbar同步更新 
clear; clc; addpath(fileparts(mfilename('fullpath')));

n=6; m=3; p=15; ell=2; N_pred=30; T_cl=2000; T_off=200;
beta_u=0.5; sw0=0.1; se0=0.3; u_min=-5; u_max=5;
Nr=30; Nw=30; gb=0.02; Nv=12;
ni=200; nf=[1,1.5,0.8,2,0.5,1.2,3,0.6,1,1.8];
rv=[1,3,0.2,3,0.5]; y_max=3.10;


rng(42,'twister');
A=diag([.85,.70,.55,.40,.30,.20]); A(1,2)=.1; A(3,4)=.01;
B=randn(n,m); C=randn(p,n); [C,~]=qr(C,0);

Rf=zeros(p,T_cl);
for s=1:5, ks=(s-1)*400+50; ke=min(s*400+49,T_cl);
    if ks<=T_cl, Rf(1:ell,ks:ke)=rv(s); end; end

oq=optimset('Display','off','LargeScale','off');
Q=zeros(p); Q(1:ell,1:ell)=eye(ell); Rw=2*eye(m); nz=ell; Ex=eye(ell);

xo=zeros(n,T_off+1); yo=zeros(p,T_off); uo=10*randn(m,T_off);
for k=1:T_off, yo(:,k)=C*xo(:,k)+se0*randn(p,1); xo(:,k+1)=A*xo(:,k)+B*uo(:,k)+sw0*randn(n,1); end
[A0,B0,P0,~,R0,~,G0,Se0,~,F0,H0]=predvarx_identify(yo,uo,ell,beta_u,2,A,B,C,n,m,p);
Ac=F0; Pbar0=null(P0');
% 为后续离线分析保留每次辨识得到的模型快照。
modelC_time=0; modelC_A={A0}; modelC_B={B0}; modelC_P={P0}; modelC_R={R0}; modelC_G={G0};
modelZ_time=0; modelZ_A={A0}; modelZ_B={B0}; modelZ_P={P0}; modelZ_R={R0}; modelZ_G={G0};
Mc=cell(1,N_pred); Nc=cell(1,N_pred);
for j=1:N_pred, Mc{j}=P0*Ex'*Ac^j; Nc{j}=zeros(p,N_pred*m);
    for i=0:j-1, Nc{j}(:,i*m+1:(i+1)*m)=P0*Ex'*Ac^(j-1-i)*H0; end; end

% ===== C (固定:累积重辨识, Pbar不更新) =====
fprintf('C ...\n'); yC=zeros(p,T_cl); uC=zeros(m,T_cl);
sC_eps=zeros(1,T_cl); sC_ebar=zeros(1,T_cl); costC=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0; PbrS=Pbar0;
cp=zeros(Nv*m,1); eb=zeros(ell,Nw); ec=0; ebb=zeros(p-ell,Nw); ecb=0; SebE=zeros(p-ell,p-ell);
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc;
for k=1:T_cl
    sw=sw0; se=se0;
    yk=C*xt+se*randn(p,1); yC(:,k)=yk; xk=Rk'*yk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
    if k>=2, ek=xk-Akh*xk_prev-Bk*uC(:,k-1);
        bi=mod(k-1,Nw)+1; eb(:,bi)=ek; ec=min(ec+1,Nw);
        if ec>=Nw, Sk=(eb*eb')/Nw; Sk=(Sk+Sk')/2; end; end
    sC_eps(k)=sqrt(trace(Sk)/ell);
    ebar=PbrS'*(yk-Pk*xk); ecb=ecb+1;
    SebE=SebE+(ebar*ebar'-SebE)/ecb; SebE=(SebE+SebE')/2;
    sC_ebar(k)=sqrt(trace(SebE)/(p-ell+1e-30));
    zk=xk; Nv=12;
    H=zeros(Nv*m); f=zeros(Nv*m,1); ya=Pk*zk;
    for j=1:Nv
        rj=Rf(:,min(k+j,T_cl))-bk; rp=Rf(:,min(k+j-1,T_cl))-bk;
        if j==1, dr=rj-ya; Md=Mk{j}; Nd=Nk{j}(:,1:Nv*m); fc=-ya;
        else, dr=rj-rp; Md=Mk{j}-Mk{j-1}; Nd=Nk{j}(:,1:Nv*m)-Nk{j-1}(:,1:Nv*m); fc=zeros(p,1); end
        H=H+Nd'*Q*Nd; f=f+Nd'*Q*(Md*zk+fc-dr);
    end
    H=H+kron(eye(Nv),Rw); H=(H+H')/2+1e-8*eye(size(H));
    lb=u_min*ones(Nv*m,1); ub=u_max*ones(Nv*m,1);
    try [co,~,~]=quadprog(H,f,[],[],[],[],lb,ub,cp(1:Nv*m),oq); catch, co=max(min(-H\f,ub),lb); end
    uk=co(1:m); uC(:,k)=uk;
    % 完整代价 = QP部分 + 常数项 Σ_j ||M_j*zk - r_j||_Q^2
    J_const = 0;
    for jj=1:Nv
        rjj=Rf(:,min(k+jj,T_cl))-bk; rpj=Rf(:,min(k+jj-1,T_cl))-bk;
        if jj==1, drj=rjj-ya; else, drj=rjj-rpj; end
        J_const = J_const + (Mk{jj}*zk - drj)'*Q*(Mk{jj}*zk - drj);
    end
    costC(k)=0.5*co'*H*co+f'*co+J_const;
    xt=A*xt+B*uk+sw*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0; cp(1:Nv*m)=max(min(cp(1:Nv*m),ub),lb);
    xk_prev=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && k>0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,beta_u,2,A,B,C,n,m,p);
        modelC_time(end+1)=k; modelC_A{end+1}=Akh; modelC_B{end+1}=Bk;
        modelC_P{end+1}=Pk; modelC_R{end+1}=Rk; modelC_G{end+1}=G0;
        Ak=Akh; for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb=zeros(ell,Nw); ec=0; ebb=zeros(p-ell,Nw); ecb=0; SebE=zeros(p-ell,p-ell); tl=k;
        % C: 累积(对标旧版)
    end
end
eC=mean(abs(yC(2,500:end)-Rf(2,500:end))); vC=sum(yC(2,:)>y_max);

% ===== SMPC-only (无软约束, 只用最近200步) =====
Sigma_v = se0^2 * eye(p);
alpha_cc = 0.45; z_score = norminv(alpha_cc);
fprintf('SMPC-only ...\n'); yZ=zeros(p,T_cl); uZ=zeros(m,T_cl);
sZ_eps=zeros(1,T_cl); sZ_ebar=zeros(1,T_cl); costZ=zeros(1,T_cl);
Sz=zeros(ell,ell,N_pred);
xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0; PbrS=Pbar0;
cp=zeros(Nv*m,1); eb=zeros(ell,Nw); ec=0; ebb=zeros(p-ell,Nw); ecb=0; SebE=zeros(p-ell,p-ell);
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc;
for k=1:T_cl
    sw=sw0; se=se0;
    yk=C*xt+se*randn(p,1); yZ(:,k)=yk; xk=Rk'*yk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
    if k>=2, ek=xk-Akh*xk_prev-Bk*uZ(:,k-1);
        bi=mod(k-1,Nw)+1; eb(:,bi)=ek; ec=min(ec+1,Nw);
        if ec>=Nw, Sk=(eb*eb')/Nw; Sk=(Sk+Sk')/2; end; end
    sig=sqrt(trace(Sk)/ell); sZ_eps(k)=sig;
    ebar=PbrS'*(yk-Pk*xk); ecb=ecb+1;
    SebE=SebE+(ebar*ebar'-SebE)/ecb; SebE=(SebE+SebE')/2;
    sZ_ebar(k)=sqrt(trace(SebE)/(p-ell+1e-30));
    zk=xk;
    % ★ 协方差传播
    Qa=G0*Sk*G0';
    if k>1, Sp=Ak*Sz(:,:,N_pred)*Ak'+Qa; Sj=Sp; else Sp=zeros(ell); Sj=zeros(ell); end
    Sz(:,:,1)=Sj; for jj=2:N_pred, Sj=Ak*Sp*Ak'+Qa; Sp=Sj; Sz(:,:,jj)=Sj; end

    Nv=12;
    H=zeros(Nv*m); f=zeros(Nv*m,1); ya=Pk*zk;
    for j=1:Nv
        rj=Rf(:,min(k+j,T_cl))-bk; rp=Rf(:,min(k+j-1,T_cl))-bk;
        if j==1, dr=rj-ya; Md=Mk{j}; Nd=Nk{j}(:,1:Nv*m); fc=-ya;
        else, dr=rj-rp; Md=Mk{j}-Mk{j-1}; Nd=Nk{j}(:,1:Nv*m)-Nk{j-1}(:,1:Nv*m); fc=zeros(p,1); end
        H=H+Nd'*Q*Nd; f=f+Nd'*Q*(Md*zk+fc-dr);
    end
    H=H+kron(eye(Nv),Rw); H=(H+H')/2+1e-8*eye(size(H));

    % ★ SMPC 机会约束: μ_y_i(j) + z_α*sqrt(Σ_y_ii(j)) ≤ y_max
    n_cc=2*ell*Nv; A_ch=zeros(n_cc,Nv*m); b_ch=zeros(n_cc,1); ri=0;
    for jj=1:min(Nv,N_pred)
        N_fj=Nk{jj}(:,1:Nv*m); mu0j=Mk{jj}*zk;
        Sig_j=Pk*Sz(:,:,jj)*Pk'+Sigma_v;
        for ch=1:ell
            tv=z_score*sqrt(max(Sig_j(ch,ch),1e-6));
            ri=ri+1; A_ch(ri,:)=N_fj(ch,:); b_ch(ri)=y_max-mu0j(ch)-tv;
            ri=ri+1; A_ch(ri,:)=-N_fj(ch,:); b_ch(ri)=mu0j(ch)-tv+y_max;
        end
    end
    A_ch=A_ch(1:ri,:); b_ch=b_ch(1:ri);

    lb=u_min*ones(Nv*m,1); ub=u_max*ones(Nv*m,1);
    try [co0,~,~]=quadprog(H,f,[],[],[],[],lb,ub,cp(1:Nv*m),oq); catch, co0=cp(1:Nv*m); end
    try [co,~,ef]=quadprog(H,f,A_ch,b_ch,[],[],lb,ub,co0,oq);
        if ef<0, co=co0; end
    catch, co=co0; end
    uk=co(1:m); uZ(:,k)=uk;
    % 完整代价 = QP部分 + 常数项 Σ_j ||state_pred_j - target_j||_Q^2
    J_const = 0;
    for jj=1:Nv
        rjj=Rf(:,min(k+jj,T_cl))-bk; rpj=Rf(:,min(k+jj-1,T_cl))-bk;
        if jj==1
            state_j = Mk{jj}*zk;
            target_j = rjj - ya;
        else
            state_j = Mk{jj}*zk - Mk{jj-1}*zk;
            target_j = rjj - rpj;
        end
        J_const = J_const + (state_j - target_j)'*Q*(state_j - target_j);
    end
    costZ(k)=0.5*co'*H*co+f'*co+J_const;
    xt=A*xt+B*uk+sw*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0; cp(1:Nv*m)=max(min(cp(1:Nv*m),ub),lb);
    xk_prev=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && k>0 && (k-tl)>=Nr
        % ★ 滑动窗口: 只用最近200步数据
        % ★ 离线200 + 最近100步
        n_recent = min(100, size(yc,2) - T_off);
        yc_in = [yc(:,1:T_off), yc(:,end-n_recent+1:end)];
        uc_in = [uc(:,1:T_off), uc(:,end-n_recent+1:end)];
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc_in,uc_in,ell,beta_u,2,A,B,C,n,m,p);
        PbrS=null(Pk');  % ★ 同步更新Pbar
        modelZ_time(end+1)=k; modelZ_A{end+1}=Akh; modelZ_B{end+1}=Bk;
        modelZ_P{end+1}=Pk; modelZ_R{end+1}=Rk; modelZ_G{end+1}=G0;
        Ak=Akh; for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb=zeros(ell,Nw); ec=0; ebb=zeros(p-ell,Nw); ecb=0; SebE=zeros(p-ell,p-ell); tl=k;
    end
end
eZ=mean(abs(yZ(2,500:end)-Rf(2,500:end))); vZ=sum(yZ(2,:)>y_max);

fprintf('C=%.3f(%d) Zfix=%.3f(%d) imp=%+.0f%% viol_red=%.0f%%\n',...
    eC,vC,eZ,vZ,(eC-eZ)/eC*100,(vC-vZ)/max(vC,1)*100);
fprintf('ebar: C=%.4f Zfix=%.4f\\n', mean(sC_ebar(500:end)), mean(sZ_ebar(500:end)));

% 保存：将轨迹、噪声统计、真实系统、离线数据和在线重辨识快照统一写入 MAT。
results_dir=fullfile(fileparts(mfilename('fullpath')),'results'); if ~exist(results_dir,'dir'), mkdir(results_dir); end;
sv.yC=yC; sv.yZ=yZ; sv.uC=uC; sv.uZ=uZ;
sv.sC_eps=sC_eps; sv.sC_ebar=sC_ebar;
sv.sZ_eps=sZ_eps; sv.sZ_ebar=sZ_ebar;
sv.true_sigma_w=sw0*ones(1,T_cl);
sv.true_sigma_e=se0*ones(1,T_cl);
sv.costC=costC; sv.costZ=costZ;
sv.Rf=Rf; sv.y_max=y_max; sv.T_cl=T_cl;
sv.eC=eC; sv.eZ=eZ; sv.vC=vC; sv.vZ=vZ;
sv.nf=nf; sv.ni=ni; sv.sw0=sw0; sv.se0=se0; sv.rv=rv;
sv.n=n; sv.m=m; sv.p=p; sv.ell=ell; sv.N_pred=N_pred; sv.Nv=Nv;
sv.T_off=T_off; sv.Nr=Nr; sv.Nw=Nw; sv.beta_u=beta_u; sv.gb=gb;
sv.u_min=u_min; sv.u_max=u_max; sv.alpha_cc=alpha_cc; sv.z_score=z_score;
sv.A_true=A; sv.B_true=B; sv.C_true=C;
sv.x_off=xo; sv.y_off=yo; sv.u_off=uo;
sv.modelC_time=modelC_time; sv.modelC_A=modelC_A; sv.modelC_B=modelC_B;
sv.modelC_P=modelC_P; sv.modelC_R=modelC_R; sv.modelC_G=modelC_G;
sv.modelZ_time=modelZ_time; sv.modelZ_A=modelZ_A; sv.modelZ_B=modelZ_B;
sv.modelZ_P=modelZ_P; sv.modelZ_R=modelZ_R; sv.modelZ_G=modelZ_G;
sv.schema_version='copyS_nosoft_analysis_v2';
save(fullfile(results_dir,'copyZfix_data.mat'),'-struct','sv','-v7.3');
fprintf('done: %s\n',fullfile(results_dir,'copyZfix_data.mat'));

%% ═══════════════════════════════════════════════════════════════
%% 画图: 4行1列纵向, 副本C蓝紫 vs Zfix红橙
%% ═══════════════════════════════════════════════════════════════
cC = [0.2 0.4 0.8];    % C 蓝
cC2 = [0.6 0.2 0.7];   % C 紫
cZ = [0.9 0.3 0.1];    % Zfix 红
cZ2 = [0.9 0.6 0.1];   % Zfix 橙
t_axis = 1:T_cl;

figure('Position',[50 50 2000 1200],'Color','w');

% --- 子图1: y1 & y2 跟踪 ---
subplot(5,1,1);
plot(t_axis, Rf(1,:), 'k:', 'LineWidth',1.2, 'DisplayName','ref_1'); hold on;
plot(t_axis, Rf(2,:), 'k--', 'LineWidth',1.5, 'DisplayName','ref_2');
plot(t_axis, yC(1,:), ':', 'Color',cC, 'LineWidth',0.6, 'DisplayName','C y_1');
plot(t_axis, yC(2,:), 'Color',cC, 'LineWidth',0.8, 'DisplayName',sprintf('C y_2 (e=%.3f)',eC));
plot(t_axis, yZ(1,:), ':', 'Color',cZ, 'LineWidth',0.6, 'DisplayName','Zfix y_1');
plot(t_axis, yZ(2,:), 'Color',cZ, 'LineWidth',0.8, 'DisplayName',sprintf('Zfix y_2 (e=%.3f)',eZ));
yline(y_max,'r--','y_max','LabelVerticalAlignment','bottom','FontSize',10);
for s=1:5, ks=(s-1)*400; xline(ks,'Color',[.7 .7 .7],'HandleVisibility','off'); end
ylabel('y'); title('y_1 (dashed, not tracked) & y_2 (solid, tracked by MPC)','FontSize',13);
legend('Location','best','FontSize',9,'NumColumns',2); grid on; hold off;

% --- 子图2: 噪声估计 vs 真实噪声（四根线） ---
subplot(5,1,2);
plot(t_axis, sC_eps, 'Color',cC, 'LineWidth',0.8, 'DisplayName',sprintf('C estimated sigma_{eps} (avg=%.3f)',mean(sC_eps(500:end)))); hold on;
plot(t_axis, sZ_eps, 'Color',cZ, 'LineWidth',0.8, 'DisplayName',sprintf('Zfix estimated sigma_{eps} (avg=%.3f)',mean(sZ_eps(500:end))));
plot(t_axis, sw0*ones(1,T_cl), 'k-', 'LineWidth',1.4, 'DisplayName',sprintf('true process noise sigma_w=%.2f',sw0));
plot(t_axis, se0*ones(1,T_cl), 'm-', 'LineWidth',1.4, 'DisplayName',sprintf('true sensor noise sigma_e=%.2f',se0));
for s=1:5, ks=(s-1)*400; xline(ks,'Color',[.7 .7 .7],'HandleVisibility','off'); end
ylabel('sigma'); title('Noise diagnosis: estimated latent process noise vs true process/sensor noise','FontSize',13);
legend('Location','best','FontSize',8,'NumColumns',2); grid on; hold off;

% --- 子图3: 控制输入 ---
subplot(5,1,3);
plot(t_axis, uC(1,:), 'Color',cC, 'LineWidth',0.6, 'DisplayName','C u_1'); hold on;
plot(t_axis, uC(2,:), '--', 'Color',cC2, 'LineWidth',0.6, 'DisplayName','C u_2');
plot(t_axis, uZ(1,:), 'Color',cZ, 'LineWidth',0.6, 'DisplayName','Zfix u_1');
plot(t_axis, uZ(2,:), '--', 'Color',cZ2, 'LineWidth',0.6, 'DisplayName','Zfix u_2');
yline(u_min,'k--','HandleVisibility','off'); yline(u_max,'k--','HandleVisibility','off');
for s=1:5, ks=(s-1)*400; xline(ks,'Color',[.7 .7 .7],'HandleVisibility','off'); end
ylabel('u'); title('控制输入','FontSize',14);
legend('Location','best','FontSize',9); grid on; hold off;

% --- 子图4: MPC 代价 ---
subplot(5,1,4);
plot(t_axis, costC, 'Color',cC, 'LineWidth',0.8, 'DisplayName',sprintf('C (avg=%.2f)',mean(costC(500:end)))); hold on;
plot(t_axis, costZ, 'Color',cZ, 'LineWidth',0.8, 'DisplayName',sprintf('Zfix (avg=%.2f)',mean(costZ(500:end))));
for s=1:5, ks=(s-1)*400; xline(ks,'Color',[.7 .7 .7],'HandleVisibility','off'); end
ylabel('J'); title('MPC full cost (Q-tracking + R-regularization)','FontSize',13);
legend('Location','best','FontSize',10); grid on; hold off;

% --- 子图5: 跟踪误差 ---
subplot(5,1,5);
plot(t_axis, abs(yC(2,:)-Rf(2,:)), 'Color',cC, 'LineWidth',0.8, 'DisplayName',sprintf('C (MAE=%.3f)',eC)); hold on;
plot(t_axis, abs(yZ(2,:)-Rf(2,:)), 'Color',cZ, 'LineWidth',0.8, 'DisplayName',sprintf('Zfix (MAE=%.3f)',eZ));
for s=1:5, ks=(s-1)*400; xline(ks,'Color',[.7 .7 .7],'HandleVisibility','off'); end
ylabel('|y_2 - ref|'); xlabel('Time step');
title('Tracking error |y_2 - ref|','FontSize',14);
legend('Location','best','FontSize',10); grid on; hold off;

sgtitle(sprintf('copyZfix: C(e=%.3f,viol=%d) vs Zfix(e=%.3f,viol=%d)',eC,vC,eZ,vZ),'FontSize',15,'FontWeight','bold');
print(fullfile(results_dir,'copySfix_fig'),'-dpng','-r150');
fprintf('figure saved: %s\n',fullfile(results_dir,'copySfix_fig.png'));