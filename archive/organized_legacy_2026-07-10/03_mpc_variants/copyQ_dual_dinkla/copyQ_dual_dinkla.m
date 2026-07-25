%% 副本 Q — 50步重辨识 + 双 Dinkla Σ_ε, Σ_ē (都从数据估)
clear; clc; addpath(fileparts(mfilename('fullpath')));

n=6; m=3; p=15; ell=2; N_pred=30; T_cl=10000; T_off=200;
beta_u=0.5; sw0=0.1; se0=0.1; u_min=-5; u_max=5;
Nr=50; Nw=50; Ns=5; Nm=12; Nx=30; gb=0.02;
ni=1000; nf=[1,1.5,0.8,2,0.5,1.2,3,0.6,1,1.8];
rv=[1,3,0.2,2.5,0.5];

rng(42,'twister');
A=diag([.85,.70,.55,.40,.30,.20]); A(1,2)=.1; A(3,4)=.01;
B=randn(n,m); C=randn(p,n); [C,~]=qr(C,0);

Rf=zeros(p,T_cl);
for s=1:5, ks=(s-1)*2000+500; ke=min(s*2000+499,T_cl);
    if ks<=T_cl, Rf(1:ell,ks:ke)=rv(s); end; end

oq=optimset('Display','off','LargeScale','off');
Q=zeros(p); Q(1:ell,1:ell)=eye(ell); Rw=eye(m);
nz=ell; Ex=eye(ell);

xo=zeros(n,T_off+1); yo=zeros(p,T_off); uo=10*randn(m,T_off);
for k=1:T_off, yo(:,k)=C*xo(:,k)+se0*randn(p,1); xo(:,k+1)=A*xo(:,k)+B*uo(:,k)+sw0*randn(n,1); end

[A0,B0,P0,~,R0,~,G0,Se0,Seb0,F0,H0]=predvarx_identify(yo,uo,ell,beta_u,2,A,B,C,n,m,p);
Ac=F0; nP=norm(P0); Pbar0=null(P0');
Mc=cell(1,N_pred); Nc=cell(1,N_pred);
for j=1:N_pred, Mc{j}=P0*Ex'*Ac^j; Nc{j}=zeros(p,N_pred*m);
    for i=0:j-1, Nc{j}(:,i*m+1:(i+1)*m)=P0*Ex'*Ac^(j-1-i)*H0; end; end

% ===== 副本 C =====
yC=zeros(p,T_cl); uC=zeros(m,T_cl); xh=zeros(ell,T_cl);
costC=zeros(1,T_cl); sCe=zeros(1,T_cl); swT=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; cp=zeros(Nx*m,1); Sz=zeros(nz,nz,N_pred);
for k=1:T_cl
    sg=min(floor((k-1)/ni)+1,length(nf));
    sw=sw0*nf(sg); se=se0*nf(sg); swT(k)=sqrt(sw^2+se^2);
    yk=C*xt+se*randn(p,1); yC(:,k)=yk; xk=R0'*yk; xh(:,k)=xk;
    if k>=2, ek=xk-A0*xh(:,k-1)-B0*uC(:,k-1); Sk=.95*Sk+.05*(ek*ek'); end
    if k>=2, sCe(k)=sqrt(trace(Sk)/ell); end
    zk=xk; Qa=G0*Sk*G0';
    if k>1, Sp=Ac*Sz(:,:,N_pred)*Ac'+Qa; Sj=Sp; else Sp=zeros(nz); Sj=zeros(nz); end
    Sz(:,:,1)=Sj; for j=2:N_pred, Sj=Ac*Sp*Ac'+Qa; Sp=Sj; Sz(:,:,j)=Sj; end
    if k<500, Nk=Ns; elseif k<=8000, Nk=Nm; else Nk=Nx; end
    H=zeros(Nk*m); f=zeros(Nk*m,1); ya=P0*xk;
    for j=1:Nk
        rj=Rf(:,min(k+j,T_cl));
        if j==1, dr=rj-ya; Md=Mc{j}; Nd=Nc{j}(:,1:Nk*m); fc=-ya;
        else, rp=Rf(:,min(k+j-1,T_cl)); dr=rj-rp; Md=Mc{j}-Mc{j-1}; Nd=Nc{j}(:,1:Nk*m)-Nc{j-1}(:,1:Nk*m); fc=zeros(p,1); end
        H=H+Nd'*Q*Nd; f=f+Nd'*Q*(Md*zk+fc-dr);
    end
    H=H+kron(eye(Nk),Rw); H=(H+H')/2+1e-8*eye(size(H));
    lb=u_min*ones(Nk*m,1); ub=u_max*ones(Nk*m,1);
    try [co,~,~]=quadprog(H,f,[],[],[],[],lb,ub,cp(1:Nk*m),oq); catch, co=max(min(-H\f,ub),lb); end
    uk=co(1:m); uC(:,k)=uk;
    cst=0; for jj=1:Nk
        rjj=Rf(:,min(k+jj,T_cl)); rpp=Rf(:,min(k+jj-1,T_cl));
        if jj==1, drr=rjj-ya; ress=Mc{jj}*zk+Nc{jj}(:,1:Nk*m)*co-ya-drr;
        else, drr=rjj-rpp; ress=(Mc{jj}-Mc{jj-1})*zk+(Nc{jj}(:,1:Nk*m)-Nc{jj-1}(:,1:Nk*m))*co-drr; end
        cst=cst+ress'*Q*ress;
    end
    costC(k)=0.5*co'*H*co+f'*co+cst;
    xt=A*xt+B*uk+sw*randn(n,1);
    cp=[co(2:end);co(end)]; cp=max(min(cp,ub),lb);
end
eC=mean(abs(yC(2,500:end)-Rf(2,500:end)));

% ===== 副本 Q: 50步 IVR + Dinkla Σ_ε + Dinkla Σ_ē =====
fprintf('副本 Q (Nr=50, 双Dinkla) ...\n');
yQ=zeros(p,T_cl); uQ=zeros(m,T_cl); xhQ=zeros(ell,T_cl);
costQ=zeros(1,T_cl); sQe=zeros(1,T_cl); sQeB=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Akh=A0; Bk=B0; Rk=R0; Pbr=Pbar0;
cp=zeros(Nx*m,1); Sz=zeros(nz,nz,N_pred);
eb_eps=zeros(ell,Nw); ec_eps=0;           % ★ Dinkla 窗口 1: Σ_ε
eb_ebar=zeros(p-ell,Nw); ec_ebar=0;       % ★ Dinkla 窗口 2: Σ_ē
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr;
Mk=Mc; Nk=Nc; SebE=zeros(p-ell,p-ell);
for k=1:T_cl
    sg=min(floor((k-1)/ni)+1,length(nf));
    sw=sw0*nf(sg); se=se0*nf(sg);
    yk=C*xt+se*randn(p,1); yQ(:,k)=yk; xk=Rk'*yk; xhQ(:,k)=xk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);

    % ★ Dinkla Σ_ε (潜空间新息)
    if k>=2, ek=xk-Akh*xhQ(:,k-1)-Bk*uQ(:,k-1);
        bi_eps=mod(k-1,Nw)+1; eb_eps(:,bi_eps)=ek; ec_eps=min(ec_eps+1,Nw);
        if ec_eps>=Nw, Sk=(eb_eps*eb_eps')/Nw; Sk=(Sk+Sk')/2; end; end
    if k>=2, sQe(k)=sqrt(trace(Sk)/(ell+1e-30)); end

    % ★ Dinkla Σ_ē (观测空间残差, 投影到 P̄ 后)
    ebar_k = Pbr'*(yk-Pk*xk);  % (p-ell)×1 = 13×1
    bi_ebar=mod(k-1,Nw)+1; eb_ebar(:,bi_ebar)=ebar_k; ec_ebar=min(ec_ebar+1,Nw);
    if ec_ebar>=Nw, SebE=(eb_ebar*eb_ebar')/Nw; SebE=(SebE+SebE')/2; end
    sQeB(k)=sqrt(trace(SebE)/(p-ell+1e-30));

    zk=xk; Qa=G0*Sk*G0';
    if k>1, Sp=Ak*Sz(:,:,N_pred)*Ak'+Qa; Sj=Sp; else Sp=zeros(nz); Sj=zeros(nz); end
    Sz(:,:,1)=Sj; for j=2:N_pred, Sj=Ak*Sp*Ak'+Qa; Sp=Sj; Sz(:,:,j)=Sj; end
    if k<500, Nv=Ns; elseif k<=8000, Nv=Nm; else Nv=Nx; end
    H=zeros(Nv*m); f=zeros(Nv*m,1); ya=Pk*xk;
    for j=1:Nv
        rj=Rf(:,min(k+j,T_cl))-bk; rp=Rf(:,min(k+j-1,T_cl))-bk;
        if j==1, dr=rj-ya; Md=Mk{j}; Nd=Nk{j}(:,1:Nv*m); fc=-ya;
        else, dr=rj-rp; Md=Mk{j}-Mk{j-1}; Nd=Nk{j}(:,1:Nv*m)-Nk{j-1}(:,1:Nv*m); fc=zeros(p,1); end
        H=H+Nd'*Q*Nd; f=f+Nd'*Q*(Md*zk+fc-dr);
    end
    % 构建 Li 2026 Eq.17a chance constraints (输出通道 1:ell)
    %  Σ_y^cl(j) = P̂·Sz·P̂^T + P̄·Σ̂_ē·P̄^T
    p_fail=0.16; beta=norminv(1-p_fail); y_max_cc=4.0; y_min_cc=-1.0;
    n_cc=2*ell*Nv; A_ch=zeros(n_cc,Nv*m); b_ch=zeros(n_cc,1); ri=0;
    for jj=1:min(Nv,N_pred)
        N_full_j=Nk{jj}(:,1:Nv*m); mu0_j=Mk{jj}*zk;
        Sig_y_j=Pk*Sz(:,:,jj)*Pk'+Pbr*SebE*Pbr';
        for ch=1:ell
            tight=beta*sqrt(max(Sig_y_j(ch,ch),1e-6));
            % Upper: e_i^T μ_y + tight ≤ y_max → N*c ≤ y_max-mu0-tight
            ri=ri+1; A_ch(ri,:)=N_full_j(ch,:); b_ch(ri)=y_max_cc-mu0_j(ch)-tight;
            % Lower: e_i^T μ_y - tight ≥ y_min → -N*c ≤ mu0-tight-y_min
            ri=ri+1; A_ch(ri,:)=-N_full_j(ch,:); b_ch(ri)=mu0_j(ch)-tight-y_min_cc;
        end
    end
    A_ch=A_ch(1:ri,:); b_ch=b_ch(1:ri);
    H=H+kron(eye(Nv),Rw); H=(H+H')/2+1e-8*eye(size(H));
    lb=u_min*ones(Nv*m,1); ub=u_max*ones(Nv*m,1);
    try [co,~,~]=quadprog(H,f,A_ch,b_ch,[],[],lb,ub,cp(1:Nv*m),oq); catch, co=max(min(-H\f,ub),lb); end
    uk=co(1:m); uQ(:,k)=uk;
    cst=0; for jj=1:Nv
        rjj=Rf(:,min(k+jj,T_cl))-bk; rpp=Rf(:,min(k+jj-1,T_cl))-bk;
        if jj==1, drr=rjj-ya; ress=Mk{jj}*zk+Nk{jj}(:,1:Nv*m)*co-ya-drr;
        else, drr=rjj-rpp; ress=(Mk{jj}-Mk{jj-1})*zk+(Nk{jj}(:,1:Nv*m)-Nk{jj-1}(:,1:Nv*m))*co-drr; end
        cst=cst+ress'*Q*ress;
    end
    costQ(k)=0.5*co'*H*co+f'*co+cst;
    xt=A*xt+B*uk+sw*randn(n,1);
    cp=[co(2:end);co(end)]; cp=max(min(cp,ub),lb);
    yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && k>0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,beta_u,2,A,B,C,n,m,p);
        Pbr=null(Pk'); Ak=Akh;
        for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb_eps=zeros(ell,Nw); ec_eps=0;
        eb_ebar=zeros(p-ell,Nw); ec_ebar=0;
        tl=k;
    end
end
eQ=mean(abs(yQ(2,500:end)-Rf(2,500:end)));
fprintf('C:%.3f Q:%.3f imp=+%.0f%%\n',eC,eQ,(eC-eQ)/eC*100);

% 图
figure('Position',[50,50,2000,1000],'Name','Copy Q');
subplot(4,1,1); hold on; grid on;
plot(1:T_cl,Rf(1,:),'k--','LineWidth',1); plot(1:T_cl,Rf(2,:),'k:','LineWidth',1);
plot(1:T_cl,yQ(1,:),'Color',[1 0.4 0],'LineWidth',1); plot(1:T_cl,yQ(2,:),'Color',[0 0.3 0.8],'LineWidth',1);
for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
ylabel('y'); title(sprintf('副本 Q: Nr=50, 双Dinkla 估算噪音 (误差=%.3f, +%.0f%%)',eQ,(eC-eQ)/eC*100));
legend('ref y_1','ref y_2','Q y_1','Q y_2','Location','best','FontSize',8);

subplot(4,1,2); hold on; grid on;
plot(1:T_cl,swT,'k-','LineWidth',1.3,'DisplayName','真实 ε (w+v)');
plot(1:T_cl,sQe,'Color',[1 0.4 0],'LineWidth',1.0,'DisplayName','Q σ_ε (Dinkla,动态)');
plot(1:T_cl,sQeB,'Color',[0 0.3 0.8],'LineWidth',1.0,'DisplayName','Q σ_ē (Dinkla,静态)');
for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
ylabel('sigma'); title('Q 噪声: 双 Dinkla 50步 (σ_ε + σ_ē)'); legend('Location','best','FontSize',8);

subplot(4,1,3); hold on; grid on;
u_Q=[0.6 0.2 0.6; 0.9 0.5 0.1; 0.3 0.6 0.6];
for ch=1:m, plot(1:T_cl,uQ(ch,:),'Color',[u_Q(ch,:) 0.6],'LineWidth',0.7); end
yline(5,'r--'); yline(-5,'r--');
for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
ylabel('u'); title('Q 控制输入');

subplot(4,1,4); hold on; grid on;
semilogy(1:T_cl,max(costQ,1e-10),'Color',[1 0.4 0],'LineWidth',1);
for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
xlabel('k'); ylabel('J_k'); title('Q 代价');
saveas(gcf,'copyQ_result.png');
fprintf('图已保存: copyQ_result.png\n');
% 保存数据 + 越界检查
y_max_val=3.2; viol=sum(yQ(2,:)>y_max_val);
fprintf('\n=== 越界: ymax=%.1f max_y=%.3f viol=%d/%d ===\n', y_max_val, max(yQ(2,:)), viol, T_cl);
sv.yQ=yQ; sv.uQ=uQ; sv.sQe=sQe; sv.sQeb=sQeb; sv.Rf=Rf; sv.eQ=eQ; sv.T_cl=T_cl; sv.y_max=y_max_val;
save('copyQ_data.mat','-struct','sv');