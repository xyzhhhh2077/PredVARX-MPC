%% 副本 D — 增量MPC (Δ-MPC) + SMPC
clear; clc; addpath(fileparts(mfilename('fullpath')));

n=6; m=3; p=15; ell=2; N_pred=30; T_cl=2000; T_off=200;
beta_u=0.5; sw0=0.1; se0=0.1; u_min=-5; u_max=5;
du_min=-2; du_max=2;  % Δu 增量约束
Nr=50; Nw=50; gb=0.02;
rv=[1,3,0.2,2.5,0.5]; y_max=3.10;
Sigma_v = se0^2 * eye(p);
alpha_cc = 0.84; z_score = norminv(alpha_cc);
Nv=12;  % 固定预测步长

rng(42,'twister');
A=diag([.85,.70,.55,.40,.30,.20]); A(1,2)=.1; A(3,4)=.01;
B=randn(n,m); C=randn(p,n); [C,~]=qr(C,0);

Rf=zeros(p,T_cl);
for s=1:5, ks=(s-1)*400+50; ke=min(s*400+49,T_cl);
    if ks<=T_cl, Rf(1:ell,ks:ke)=rv(s); end; end

oq=optimset('Display','off','LargeScale','off','Algorithm','interior-point-convex');
Q=zeros(p); Q(1:ell,1:ell)=eye(ell); Rw=2*eye(m);
nz_aug = ell + m;  % 增广状态维度: x̂(ℓ) + u_prev(m)
Ex_aug = [eye(ell), zeros(ell,m)];  % 从增广状态提取 x̂

xo=zeros(n,T_off+1); yo=zeros(p,T_off); uo=10*randn(m,T_off);
for k=1:T_off, yo(:,k)=C*xo(:,k)+se0*randn(p,1); xo(:,k+1)=A*xo(:,k)+B*uo(:,k)+sw0*randn(n,1); end
[A0,B0,P0,~,R0,~,G0,Se0,~,F0,H0]=predvarx_identify(yo,uo,ell,beta_u,2,A,B,C,n,m,p);

% ★ 构建增广动力学矩阵
% z = [x̂; u_prev], z(k+1) = A_aug*z(k) + B_aug*Δu(k)
% A_aug = [Â, B̂; 0, I], B_aug = [B̂; I]
A_aug = [A0, B0; zeros(m,ell), eye(m)];
B_aug = [B0; eye(m)];
G_aug = [G0; zeros(m,ell)];  % 噪声映射

% ★ 增广预测矩阵: μ_y(j) = P̂*Ex_aug*A_aug^j*z_k + Σ P̂*Ex_aug*A_aug^(j-1-i)*B_aug*Δu_i
Mc=cell(1,N_pred); Nc=cell(1,N_pred);
for j=1:N_pred
    Mc{j}=P0*Ex_aug*A_aug^j;
    Nc{j}=zeros(p,N_pred*m);
    for i=0:j-1
        Nc{j}(:,i*m+1:(i+1)*m)=P0*Ex_aug*A_aug^(j-1-i)*B_aug;
    end
end

%% ===== Δ-MPC =====
fprintf('Delta-MPC ...\n');
yD=zeros(p,T_cl); uD=zeros(m,T_cl); sD_eps=zeros(1,T_cl); sD_ebar=zeros(1,T_cl);
costD=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; Ak=A_aug; Pk=P0; Rk=R0;
Akh=A0; Bk=B0; PbrS=null(P0');
cp=zeros(Nv*m,1); eb=zeros(ell,Nw); ec=0;
ecb=0; SebE=zeros(p-ell,p-ell);
bk=zeros(p,1); yc=zeros(p,0); uc=zeros(m,0); tl=-Nr;
Mk=Mc; Nk=Nc; u_prev=zeros(m,1);  % ★ 初始 u(k-1)=0

for k=1:T_cl
    yk=C*xt+se0*randn(p,1); yD(:,k)=yk; xk=Rk'*yk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
    if k>=2, ek=xk-Akh*xk_prev-Bk*uD(:,k-1);
        bi=mod(k-1,Nw)+1; eb(:,bi)=ek; ec=min(ec+1,Nw);
        if ec>=Nw, Sk=(eb*eb')/Nw; Sk=(Sk+Sk')/2; end; end
    sD_eps(k)=sqrt(trace(Sk)/ell);
    ebar=PbrS'*(yk-Pk*xk); ecb=ecb+1;
    SebE=SebE+(ebar*ebar'-SebE)/ecb; SebE=(SebE+SebE')/2;
    sD_ebar(k)=sqrt(trace(SebE)/(p-ell+1e-30));

    % ★ 增广状态: z = [x̂; u_prev]
    zk = [xk; u_prev];
    Qa=G_aug*Sk*G_aug';
    if k>1, Sp=Ak*Sz_aug(:,:,N_pred)*Ak'+Qa; Sj=Sp; else Sp=zeros(nz_aug); Sj=zeros(nz_aug); end
    Sz_aug(:,:,1)=Sj; for j=2:N_pred, Sj=Ak*Sp*Ak'+Qa; Sp=Sj; Sz_aug(:,:,j)=Sj; end

    H=zeros(Nv*m); f=zeros(Nv*m,1); ya=Pk*Ex_aug*zk;
    for j=1:Nv
        rj=Rf(:,min(k+j,T_cl))-bk; rp=Rf(:,min(k+j-1,T_cl))-bk;
        if j==1, dr=rj-ya; Md=Mk{j}; Nd=Nk{j}(:,1:Nv*m); fc=-ya;
        else, dr=rj-rp; Md=Mk{j}-Mk{j-1}; Nd=Nk{j}(:,1:Nv*m)-Nk{j-1}(:,1:Nv*m); fc=zeros(p,1); end
        H=H+Nd'*Q*Nd; f=f+Nd'*Q*(Md*zk+fc-dr);
    end
    H=H+kron(eye(Nv),Rw); H=(H+H')/2+1e-8*eye(size(H));

    % 只用 Δu box 约束，不加 u 绝对约束和 SMPC
    lb_du=du_min*ones(Nv*m,1); ub_du=du_max*ones(Nv*m,1);

    % 求解 QP
    try [co,~,~]=quadprog(H,f,[],[],[],[],lb_du,ub_du,cp(1:Nv*m),oq); catch, co=max(min(-H\f,ub_du),lb_du); end

    du_k=co(1:m);  % ★ 当前增量
    uk=u_prev+du_k;  % ★ u(k) = u(k-1) + Δu(k)
    uk=max(min(uk,u_max),u_min);  % 安全截断
    uD(:,k)=uk; u_prev=uk;  % ★ 更新 u_prev

    J_const=0;
    for jj=1:Nv
        rjj=Rf(:,min(k+jj,T_cl))-bk; rpj=Rf(:,min(k+jj-1,T_cl))-bk;
        if jj==1, drj=rjj-ya; else, drj=rjj-rpj; end
        J_const=J_const+(Mk{jj}*zk-drj)'*Q*(Mk{jj}*zk-drj);
    end
    costD(k)=0.5*co'*H*co+f'*co+J_const;

    xt=A*xt+B*uk+sw0*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0;
    cp(1:Nv*m)=max(min(cp(1:Nv*m),ub_du),lb_du);
    xk_prev=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && k>0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,beta_u,2,A,B,C,n,m,p);
        PbrS=null(Pk');
        % 重建增广矩阵
        A_aug_new = [Akh, Bk; zeros(m,ell), eye(m)];
        B_aug_new = [Bk; eye(m)];
        G_aug_new = [G0; zeros(m,ell)];
        Ak=A_aug_new; G_aug=G_aug_new;
        for j=1:N_pred, Mk{j}=Pk*Ex_aug*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex_aug*Ak^(j-1-i)*B_aug_new; end; end
        eb=zeros(ell,Nw); ec=0; ecb=0; SebE=zeros(p-ell,p-ell); tl=k;
    end
end
eD=mean(abs(yD(2,500:end)-Rf(2,500:end))); vD=sum(yD(2,:)>y_max);

%% 结果
fprintf('\n=== Delta-MPC 结果 ===\n');
fprintf('e=%.3f viol=%d ebar=%.4f\n', eD, vD, mean(sD_ebar(500:end)));
fprintf('du range: [%.3f, %.3f]\n', min(diff([zeros(m,1), uD],1,2),[],'all'), max(diff([zeros(m,1), uD],1,2),[],'all'));

%% 画图
cD=[0.2 0.6 0.2]; t_axis=1:T_cl;
figure('Position',[50 50 2000 1200],'Color','w');

subplot(5,1,1);
plot(t_axis,Rf(2,:),'k--','LineWidth',1.5,'DisplayName','ref_2'); hold on;
plot(t_axis,yD(2,:),'Color',cD,'LineWidth',0.8,'DisplayName',sprintf('Delta-MPC (e=%.3f,v=%d)',eD,vD));
yline(y_max,'r--','y_max','FontSize',10,'LabelVerticalAlignment','bottom');
for s=1:5, xline((s-1)*400,'Color',[.7 .7 .7],'HandleVisibility','off'); end
ylabel('y_2'); title('y_2 tracking: Delta-MPC','FontSize',13);
legend('Location','best','FontSize',9); grid on; hold off;

subplot(5,1,2);
plot(t_axis,sD_eps,'Color',cD,'LineWidth',0.8,'DisplayName','sigma_eps'); hold on;
plot(t_axis,sD_ebar,'--','Color',cD,'LineWidth',0.8,'DisplayName','sigma_ebar');
yline(sw0,'g-','true sw','LineWidth',1.5); yline(se0,'m-','true se','LineWidth',1.5);
for s=1:5, xline((s-1)*400,'Color',[.7 .7 .7],'HandleVisibility','off'); end
ylabel('sigma'); title('Noise estimation','FontSize',13);
legend('Location','best','FontSize',9); grid on; hold off;

subplot(5,1,3);
plot(t_axis,uD(1,:),'Color',cD,'LineWidth',0.6,'DisplayName','u_1'); hold on;
plot(t_axis,uD(2,:),'--','Color',cD,'LineWidth',0.6,'DisplayName','u_2');
yline(u_min,'k--'); yline(u_max,'k--');
for s=1:5, xline((s-1)*400,'Color',[.7 .7 .7],'HandleVisibility','off'); end
ylabel('u'); title('Control input u (absolute)','FontSize',13);
legend('Location','best','FontSize',9); grid on; hold off;

subplot(5,1,4);
du_all = diff([zeros(m,1), uD], 1, 2);
plot(t_axis,du_all(1,:),'Color',cD,'LineWidth',0.6,'DisplayName','du_1'); hold on;
plot(t_axis,du_all(2,:),'--','Color',cD,'LineWidth',0.6,'DisplayName','du_2');
yline(du_min,'k--'); yline(du_max,'k--');
for s=1:5, xline((s-1)*400,'Color',[.7 .7 .7],'HandleVisibility','off'); end
ylabel('du'); title('Control increment Delta-u','FontSize',13);
legend('Location','best','FontSize',9); grid on; hold off;

subplot(5,1,5);
plot(t_axis,abs(yD(2,:)-Rf(2,:)),'Color',cD,'LineWidth',0.8,'DisplayName','|y_2-ref|');
for s=1:5, xline((s-1)*400,'Color',[.7 .7 .7],'HandleVisibility','off'); end
ylabel('|y_2-ref|'); xlabel('Time step'); title('Tracking error','FontSize',13);
legend('Location','best','FontSize',9); grid on; hold off;

sgtitle(sprintf('copyD: Delta-MPC (nz=%d, du=[%.1f,%.1f], alpha=%.2f)',nz_aug,du_min,du_max,alpha_cc),'FontSize',14,'FontWeight','bold');
print('copyD_fig','-dpng','-r150');
fprintf('figure saved: copyD_fig.png\n');
