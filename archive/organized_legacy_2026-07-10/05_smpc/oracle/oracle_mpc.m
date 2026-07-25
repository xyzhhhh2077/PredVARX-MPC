%% 副本 Oracle — 真系统参数(A,B,C已知) → MPC vs SMPC
%  跳过predvarx_identify, 直接用地真模型
clear; clc;

n=6; m=3; p=15; ell=2; N_pred=30; T_cl=10000;
sw0=0.1; se0=0.1; u_min=-5; u_max=5;
Ns=5; Nm=12; Nx=30; gb=0.02;
y_max=3.2; rv=[1,3,0.2,2.5,0.5];

rng(42); A=diag([.85,.70,.55,.40,.30,.20]); A(1,2)=.1; A(3,4)=.01;
B=randn(n,m); C=randn(p,n); [C,~]=qr(C,0);

Rf=zeros(p,T_cl);
for s=1:5, ks=(s-1)*2000+500; ke=min(s*2000+499,T_cl);
    if ks<=T_cl, Rf(1:ell,ks:ke)=rv(s); end; end

oq=optimset('Display','off','LargeScale','off');
Q=zeros(p); Q(1:ell,1:ell)=eye(ell); Rw=eye(m);
Ex=eye(ell); Pi=eye(ell);

% ★ 用真C的前ℓ列作为P̂ (简化: 选前2个输出方向)
[Uc,Sc,~]=svd(C*(A*C'),'econ');
P0=Uc(:,1:ell);  % p×ℓ, 最可预测方向
R0=(P0'*P0)\P0';  % ℓ×p
A0=R0*C*A*pinv(R0*C);  % ℓ×ℓ (投影后的A)
B0=R0*C*B;             % ℓ×m (投影后的B)
G0=eye(ell);
Pbr0=null(P0');
Se0=sw0^2*eye(ell);  % 已知σ_w²

Ac=A0;
Mc=cell(1,N_pred); Nc=cell(1,N_pred);
for j=1:N_pred, Mc{j}=P0*Ex'*Ac^j; Nc{j}=zeros(p,N_pred*m);
    for i=0:j-1, Nc{j}(:,i*m+1:(i+1)*m)=P0*Ex'*Ac^(j-1-i)*B0; end; end

% ==== C (纯MPC) ====
fprintf('C ...\n'); yC=zeros(p,T_cl); uC=zeros(m,T_cl);
costC=zeros(1,T_cl); sCe=zeros(1,T_cl);
xt=zeros(n,1); cp=zeros(Nx*m,1); Sz=zeros(ell,ell,N_pred);
bk=zeros(p,1); xp=zeros(ell,1); Nw=50; eb=zeros(ell,Nw); ec=0;
Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0; Mk=Mc; Nk=Nc;

for k=1:T_cl
    yk=C*xt+se0*randn(p,1); yC(:,k)=yk; xk=Rk*yk;
    % Dinkla Σ_ε
    eps_pred=Rk*(yk-Pk*(Akh*xp+Bk*uC(:,max(1,k-1))));
    bi=mod(k-1,Nw)+1; eb(:,bi)=eps_pred; ec=min(ec+1,Nw);
    Sk=Se0; if ec>=Nw, Sk=(eb*eb')/Nw; end
    sCe(k)=sqrt(trace(Sk)/ell);
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
    try co=quadprog(H,f,[],[],[],[],u_min*ones(Nv*m,1),u_max*ones(Nv*m,1),cp(1:Nv*m),oq); catch, co=cp(1:Nv*m); end
    uk=co(1:m); uC(:,k)=uk; costC(k)=0.5*co'*H*co+f'*co;
    xt=A*xt+B*uk+sw0*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0;
    xp=xk;
end
eC=mean(abs(yC(2,500:end)-Rf(2,500:end))); vC=sum(yC(2,:)>y_max);

% ==== SM (SMPC: 真Σ_e=se0²·I) ====
fprintf('SM ...\n'); yS=zeros(p,T_cl); uS=zeros(m,T_cl);
costS=zeros(1,T_cl); sSe=zeros(1,T_cl); bnd=zeros(1,T_cl);
xt=zeros(n,1); cp=zeros(Nx*m,1); Sz=zeros(ell,ell,N_pred);
bk=zeros(p,1); xp=zeros(ell,1); eb=zeros(ell,Nw); ec=0;
Sigma_v=se0^2*eye(p);
Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0; Mk=Mc; Nk=Nc;

for k=1:T_cl
    yk=C*xt+se0*randn(p,1); yS(:,k)=yk; xk=Rk*yk;
    eps_pred=Rk*(yk-Pk*(Akh*xp+Bk*uS(:,max(1,k-1))));
    bi=mod(k-1,Nw)+1; eb(:,bi)=eps_pred; ec=min(ec+1,Nw);
    Sk=Se0; if ec>=Nw, Sk=(eb*eb')/Nw; end
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

    % SMPC: Σ_y = P Σ_z P^T + Σ_v
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
    lb=u_min*ones(Nv*m,1); ub=u_max*ones(Nv*m,1);
    try co0=quadprog(Hf,f,[],[],[],[],lb,ub,cp(1:Nv*m),oq); catch, co0=cp(1:Nv*m); end
    try [co,~,ef]=quadprog(Hf,f,A_cc,b_cc,[],[],lb,ub,co0,oq);
        bnd(k)=double(ef==1 && any(abs(A_cc*co-b_cc)<1e-4));
        if ef<0, co=co0; end
    catch, co=co0; end
    uk=co(1:m); uS(:,k)=uk; costS(k)=0.5*co'*Hf*co+f'*co;
    xt=A*xt+B*uk+sw0*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0;
    xp=xk;
end
eS=mean(abs(yS(2,500:end)-Rf(2,500:end))); vS=sum(yS(2,:)>y_max);

fprintf('C=%.3f(%d) SM=%.3f(%d) imp=+%.0f%% bnd=%d\n',eC,vC,eS,vS,(eC-eS)/eC*100,sum(bnd));

% 图
figure('Position',[50,50,2000,1600]);
subplot(5,1,1);hold on;grid on;
plot(1:T_cl,Rf(2,:),'k:',1:T_cl,yC(2,:),'b-',1:T_cl,yS(2,:),'r-','LineWidth',1);
yline(y_max,'g--');for ss=2000:2000:T_cl,xline(ss,'Color',[.85 .85 .85]);end
legend('ref',sprintf('C(%.3f,%d)',eC,vC),sprintf('SM(%.3f,%d)',eS,vS),'y_m_a_x','Location','best');
ylabel('y_2');title(sprintf('Oracle: C=%.3f(%d) SM=%.3f(%d) bnd=%d',eC,vC,eS,vS,sum(bnd)));

subplot(5,1,2);hold on;grid on;
plot(500:T_cl,sCe(500:end),'b-',500:T_cl,sSe(500:end),'r-','LineWidth',1);
yline(sw0,'k-');for ss=2000:2000:T_cl,xline(ss,'Color',[.85 .85 .85]);end
legend('C','SM','true','Location','best');ylabel('sigma_eps');

subplot(5,1,3);hold on;grid on;
plot(1:T_cl,uC(1,:),'b-',1:T_cl,uS(1,:),'r--','LineWidth',1);
for ss=2000:2000:T_cl,xline(ss,'Color',[.85 .85 .85]);end
ylabel('u_1');title('Control');

subplot(5,1,4);hold on;grid on;
semilogy(500:T_cl,max(costC(500:end),1e-10),'b-',500:T_cl,max(costS(500:end),1e-10),'r-','LineWidth',1);
for ss=2000:2000:T_cl,xline(ss,'Color',[.85 .85 .85]);end
legend('C','SM','Location','best');ylabel('J_k (log)');

subplot(5,1,5);hold on;grid on;
bw=100; bb=zeros(1,T_cl-bw+1);
for i=1:T_cl-bw+1, bb(i)=sum(bnd(i:i+bw-1))/bw; end
plot(bw:T_cl,bb,'k-','LineWidth',1);
ylabel('binding rate');xlabel('k');
title(sprintf('SMPC binding: %d/%d',sum(bnd),T_cl));

sgtitle('Oracle: 真A,B,C已知 → MPC vs SMPC');
saveas(gcf,'oracle_smpc.png');
fprintf('图: oracle_smpc.png\\n');

sv.yC=yC;sv.yS=yS;sv.uC=uC;sv.uS=uS;sv.eC=eC;sv.eS=eS;sv.vC=vC;sv.vS=vS;
sv.sCe=sCe;sv.sSe=sSe;sv.costC=costC;sv.costS=costS;sv.bnd=bnd;
sv.Rf=Rf;sv.y_max=y_max;
save('oracle_data.mat','-struct','sv');
fprintf('done\\n');