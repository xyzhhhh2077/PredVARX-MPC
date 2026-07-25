%% 副本 QZretro — 基于copyQ参数 + 已知Sigma_v SMPC
%  跳过 predvarx_identify, 直接用已知 A,B,C 辨识 → 测试 SMPC
clear; clc; addpath(fileparts(mfilename('fullpath')));

n=6; m=3; p=15; ell=2; N_pred=30; T_cl=10000; T_off=200;
sw0=0.1; se0=0.1; u_min=-5; u_max=5;
Nr=50; Nw=50; Ns=5; Nm=12; Nx=30; gb=0.02;
y_max=3.2;  Sigma_v=se0^2*eye(p);  % Li2026: 已知测量噪声

rng(42); A=diag([.85,.70,.55,.40,.30,.20]); A(1,2)=.1; A(3,4)=.01;
B=randn(n,m); C=randn(p,n); [C,~]=qr(C,0);

Rf=zeros(p,T_cl); rv=[1,3,0.2,2.5,0.5];
for s=1:5, ks=(s-1)*2000+500; ke=min(s*2000+499,T_cl);
    if ks<=T_cl, Rf(1:ell,ks:ke)=rv(s); end; end

oq=optimset('Display','off','LargeScale','off');
Q=zeros(p); Q(1:ell,1:ell)=eye(ell); Rw=eye(m); Ex=eye(ell);

% 离线数据生成
xo=zeros(n,T_off+1); yo=zeros(p,T_off); uo=20*randn(m,T_off);
for k=1:T_off, yo(:,k)=C*xo(:,k)+se0*randn(p,1); xo(:,k+1)=A*xo(:,k)+B*uo(:,k)+sw0*randn(n,1); end

% 辨识 (用当前 predvarx_identify)
[A0,B0,P0,~,R0,~,G0,Se0,~,F0,H0]=predvarx_identify(yo,uo,ell,0.5,2,A,B,C,n,m,p);
Ac=F0; Pbr0=null(P0');
Mc=cell(1,N_pred); Nc=cell(1,N_pred);
for j=1:N_pred, Mc{j}=P0*Ex'*Ac^j; Nc{j}=zeros(p,N_pred*m);
    for i=0:j-1, Nc{j}(:,i*m+1:(i+1)*m)=P0*Ex'*Ac^(j-1-i)*H0; end; end

% ==== Q (纯MPC, Rw=I) ====
fprintf('Q ...\n'); yQ=zeros(p,T_cl); uQ=zeros(m,T_cl);
sQe=zeros(1,T_cl); costQ=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0; PbrS=Pbr0;
cp=zeros(Nx*m,1); eb=zeros(ell,Nw); ec=0; Sz=zeros(ell,ell,N_pred);
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc; xp=zeros(ell,1);

for k=1:T_cl
    yk=C*xt+se0*randn(p,1); yQ(:,k)=yk; xk=Rk*yk;
    eps=yk-Pk*(Akh*xp+Bk*uQ(:,max(1,k-1)));
    bi=mod(k-1,Nw)+1; eb(:,bi)=Rk*eps; ec=min(ec+1,Nw);
    if ec>=Nw, Sk=(eb*eb')/Nw; Sk=(Sk+Sk')/2; end
    sQe(k)=sqrt(trace(Sk)/ell);
    Qa=G0*Sk*G0';
    if k>1, Sp=Ak*Sz(:,:,N_pred)*Ak'+Qa; Sj=Sp; else Sp=zeros(ell); Sj=zeros(ell); end
    Sz(:,:,1)=Sj; for j=2:N_pred, Sj=Ak*Sp*Ak'+Qa; Sp=Sj; Sz(:,:,j)=Sj; end
    if k<500, Nv=Ns; elseif k<=8000, Nv=Nm; else Nv=Nx; end
    H=zeros(Nv*m); f=zeros(Nv*m,1); ya=Pk*xk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
    for j=1:Nv
        rj=Rf(:,min(k+j,T_cl))-bk; rp=Rf(:,min(k+j-1,T_cl))-bk;
        if j==1, dr=rj-ya; Md=Mk{j}; Nd=Nk{j}(:,1:Nv*m); fc=-ya;
        else, dr=rj-rp; Md=Mk{j}-Mk{j-1}; Nd=Nk{j}(:,1:Nv*m)-Nk{j-1}(:,1:Nv*m); fc=zeros(p,1); end
        H=H+Nd'*Q*Nd; f=f+Nd'*Q*(Md*xk+fc-dr);
    end
    H=H+kron(eye(Nv),Rw); H=(H+H')/2+1e-8*eye(size(H));
    try [co,~,~]=quadprog(H,f,[],[],[],[],u_min*ones(Nv*m,1),u_max*ones(Nv*m,1),cp(1:Nv*m),oq); catch, co=cp(1:Nv*m); end
    uk=co(1:m); uQ(:,k)=uk; costQ(k)=0.5*co'*H*co+f'*co;
    xt=A*xt+B*uk+sw0*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0;
    xp=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,0.5,2,A,B,C,n,m,p);
        PbrS=null(Pk'); Ak=Akh;
        for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb=zeros(ell,Nw); ec=0; tl=k;
    end
end
eQ=mean(abs(yQ(2,500:end)-Rf(2,500:end))); vQ=sum(yQ(2,:)>y_max);

% ==== SM (SMPC: 已知 Σ_v) ====
fprintf('SM ...\n'); yS=zeros(p,T_cl); uS=zeros(m,T_cl);
sSe=zeros(1,T_cl); costS=zeros(1,T_cl); bnd=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0; PbrS=Pbr0;
cp=zeros(Nx*m,1); eb=zeros(ell,Nw); ec=0; Sz=zeros(ell,ell,N_pred);
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc; xp=zeros(ell,1);

for k=1:T_cl
    yk=C*xt+se0*randn(p,1); yS(:,k)=yk; xk=Rk*yk;
    eps=yk-Pk*(Akh*xp+Bk*uS(:,max(1,k-1)));
    bi=mod(k-1,Nw)+1; eb(:,bi)=Rk*eps; ec=min(ec+1,Nw);
    if ec>=Nw, Sk=(eb*eb')/Nw; Sk=(Sk+Sk')/2; end
    sSe(k)=sqrt(trace(Sk)/ell);
    Qa=G0*Sk*G0';
    if k>1, Sp=Ak*Sz(:,:,N_pred)*Ak'+Qa; Sj=Sp; else Sp=zeros(ell); Sj=zeros(ell); end
    Sz(:,:,1)=Sj; for j=2:N_pred, Sj=Ak*Sp*Ak'+Qa; Sp=Sj; Sz(:,:,j)=Sj; end
    if k<500, Nv=Ns; elseif k<=8000, Nv=Nm; else Nv=Nx; end
    H=zeros(Nv*m); f=zeros(Nv*m,1); ya=Pk*xk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
    for j=1:Nv
        rj=Rf(:,min(k+j,T_cl))-bk; rp=Rf(:,min(k+j-1,T_cl))-bk;
        if j==1, dr=rj-ya; Md=Mk{j}; Nd=Nk{j}(:,1:Nv*m); fc=-ya;
        else, dr=rj-rp; Md=Mk{j}-Mk{j-1}; Nd=Nk{j}(:,1:Nv*m)-Nk{j-1}(:,1:Nv*m); fc=zeros(p,1); end
        H=H+Nd'*Q*Nd; f=f+Nd'*Q*(Md*xk+fc-dr);
    end

    % ★ SMPC: A_cc*c <= b_cc
    ncc=2*ell*N_pred; A_cc=zeros(ncc,Nv*m); b_cc=zeros(ncc,1); ri=0;
    for jj=1:min(Nv,N_pred)
        Nfj=Nk{jj}(:,1:Nv*m); mu0j=Mk{jj}*xk;
        Sig_j=Pk*Sz(:,:,jj)*Pk'+Sigma_v;
        for ch=1:ell
            tv=norminv(0.84)*sqrt(max(Sig_j(ch,ch),1e-6));
            ri=ri+1; A_cc(ri,:)= Nfj(ch,:); b_cc(ri)=y_max-mu0j(ch)-tv;
            ri=ri+1; A_cc(ri,:)=-Nfj(ch,:); b_cc(ri)=mu0j(ch)-tv+3.0;
        end
    end
    A_cc=A_cc(1:ri,:); b_cc=b_cc(1:ri);
    Hf=H+kron(eye(Nv),Rw); Hf=(Hf+Hf')/2+1e-8*eye(size(Hf));
    try [co0,~,~]=quadprog(Hf,f,[],[],[],[],u_min*ones(Nv*m,1),u_max*ones(Nv*m,1),cp(1:Nv*m),oq); catch, co0=cp(1:Nv*m); end
    try [co,~,ef]=quadprog(Hf,f,A_cc,b_cc,[],[],u_min*ones(Nv*m,1),u_max*ones(Nv*m,1),co0,oq);
        bnd(k)=double(ef==1 && any(abs(A_cc*co-b_cc)<1e-4));
    catch, co=co0; end
    uk=co(1:m); uS(:,k)=uk; costS(k)=0.5*co'*Hf*co+f'*co;
    xt=A*xt+B*uk+sw0*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0;
    xp=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,0.5,2,A,B,C,n,m,p);
        PbrS=null(Pk'); Ak=Akh;
        for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb=zeros(ell,Nw); ec=0; tl=k;
    end
end
eS=mean(abs(yS(2,500:end)-Rf(2,500:end))); vS=sum(yS(2,:)>y_max);

fprintf('Q=%.3f(%d) SM=%.3f(%d) imp=+%.0f%% bnd=%d\n',eQ,vQ,eS,vS,(eQ-eS)/eQ*100,sum(bnd));
fprintf('sQe=%.3f sSe=%.3f\\n',mean(sQe(500:end)),mean(sSe(500:end)));

% 5行图
figure('Position',[50,50,2000,1600]);
subplot(5,1,1);hold on;grid on;
plot(1:T_cl,Rf(2,:),'k:',1:T_cl,yQ(2,:),'b-',1:T_cl,yS(2,:),'r-','LineWidth',1);
yline(y_max,'g--');for ss=2000:2000:T_cl,xline(ss,'Color',[.85 .85 .85]);end
legend('ref',sprintf('Q(%.3f,%d)',eQ,vQ),sprintf('SM(%.3f,%d)',eS,vS),'y_m_a_x','Location','best');
ylabel('y_2');title(sprintf('Q=%.3f(%d) SM=%.3f(%d) bnd=%d',eQ,vQ,eS,vS,sum(bnd)));

subplot(5,1,2);hold on;grid on;
plot(500:T_cl,sQe(500:end),'b-',500:T_cl,sSe(500:end),'r-','LineWidth',1);
yline(sw0,'k-');for ss=2000:2000:T_cl,xline(ss,'Color',[.85 .85 .85]);end
legend('Q','SM','true sw','Location','best');ylabel('sigma_eps');

subplot(5,1,3);hold on;grid on;
plot(1:T_cl,uQ(1,:),'b-',1:T_cl,uQ(2,:),'c-',1:T_cl,uQ(3,:),'g-',1:T_cl,uS(1,:),'r--','LineWidth',1);
for ss=2000:2000:T_cl,xline(ss,'Color',[.85 .85 .85]);end
ylabel('u');title('Q(u1-u3) vs SM(u1)');

subplot(5,1,4);hold on;grid on;
semilogy(500:T_cl,max(costQ(500:end),1e-10),'b-',500:T_cl,max(costS(500:end),1e-10),'r-','LineWidth',1);
for ss=2000:2000:T_cl,xline(ss,'Color',[.85 .85 .85]);end
legend('Q','SM','Location','best');ylabel('J_k (log)');

subplot(5,1,5);hold on;grid on;
bw=100; bb=zeros(1,T_cl-bw+1);
for i=1:T_cl-bw+1, bb(i)=sum(bnd(i:i+bw-1))/bw; end
plot(bw:T_cl,bb,'k-','LineWidth',1);
yline(0,'k--');ylabel('binding rate (100-window)');xlabel('k');
title(sprintf('SMPC binding: total %d/%d',sum(bnd),T_cl));

sgtitle('Q vs SM: 已知 Σ_v=0.01·I, Rw=I');
saveas(gcf,'qzretro_smpc.png');
fprintf('图: qzretro_smpc.png\n');

sv.yQ=yQ;sv.yS=yS;sv.uQ=uQ;sv.uS=uS;sv.eQ=eQ;sv.eS=eS;sv.vQ=vQ;sv.vS=vS;
sv.sQe=sQe;sv.sSe=sSe;sv.costQ=costQ;sv.costS=costS;sv.bnd=bnd;
sv.Rf=Rf;sv.y_max=y_max;
save('qzretro_data.mat','-struct','sv');
fprintf('done\\n');