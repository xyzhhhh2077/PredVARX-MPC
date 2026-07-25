%% 副本 FINAL2 — C(纯MPC) vs SM(SMPC) + 准确噪声 + 图
clear; clc; addpath(fileparts(mfilename('fullpath')));

n=6; m=3; p=15; ell=2; N_pred=30; T_cl=10000; T_off=200;
beta_u=0.5; sw0=0.1; se0=0.1; u_min=-5; u_max=5;
Nr=50; Nw=50; Ns=5; Nm=12; Nx=30; gb=0.02;
y_max=3.2;  Sigma_v=se0^2*eye(p);  % 已知测量噪声

rng(42,'twister');
A=diag([.85,.70,.55,.40,.30,.20]); A(1,2)=.1; A(3,4)=.01;
B=randn(n,m); C=randn(p,n); [C,~]=qr(C,0);

Rf=zeros(p,T_cl); rv=[1,3,0.2,2.5,0.5];
for s=1:5, ks=(s-1)*2000+500; ke=min(s*2000+499,T_cl);
    if ks<=T_cl, Rf(1:ell,ks:ke)=rv(s); end; end

oq=optimset('Display','off','LargeScale','off','Algorithm','interior-point-convex');
Q=zeros(p); Q(1:ell,1:ell)=eye(ell); Rw=eye(m); nz=ell; Ex=eye(ell);

xo=zeros(n,T_off+1); yo=zeros(p,T_off); uo=10*randn(m,T_off);
for k=1:T_off, yo(:,k)=C*xo(:,k)+se0*randn(p,1); xo(:,k+1)=A*xo(:,k)+B*uo(:,k)+sw0*randn(n,1); end
[A0,B0,P0,~,R0,~,G0,Se0,~,F0,H0]=predvarx_identify(yo,uo,ell,beta_u,2,A,B,C,n,m,p);
Ac=F0; Pbar0=null(P0');
Mc=cell(1,N_pred); Nc=cell(1,N_pred);
for j=1:N_pred, Mc{j}=P0*Ex'*Ac^j; Nc{j}=zeros(p,N_pred*m);
    for i=0:j-1, Nc{j}(:,i*m+1:(i+1)*m)=P0*Ex'*Ac^(j-1-i)*H0; end; end

% ==== C (纯MPC) ====
fprintf('C ...\n'); yC=zeros(p,T_cl); uC=zeros(m,T_cl);
sC_eps=zeros(1,T_cl); sC_xi=zeros(1,T_cl); costC=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0; PbrS=Pbar0;
cp=zeros(Nx*m,1); eb_eps=zeros(ell,Nw); ec_eps=0; eb_xi=zeros(p-ell,Nw); ec_xi=0;
Sxi=zeros(p-ell,p-ell); Sz=zeros(nz,nz,N_pred);
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc; xk_prev=zeros(ell,1);

for k=1:T_cl
    yk=C*xt+se0*randn(p,1); yC(:,k)=yk; xk=Rk'*yk;
    v_hat=Akh*xk_prev+Bk*uC(:,max(1,k-1));
    e_full=yk-Pk*v_hat; eps_k=Rk'*e_full; xi_k=PbrS'*e_full;
    bi=mod(k-1,Nw)+1; eb_eps(:,bi)=eps_k; ec_eps=min(ec_eps+1,Nw);
    if ec_eps>=Nw, Sk=(eb_eps*eb_eps')/Nw; Sk=(Sk+Sk')/2; end
    sC_eps(k)=sqrt(trace(Sk)/ell); sC_xi(k)=sqrt(trace(Sxi)/max(p-ell,1));
    zk=xk; Qa=G0*Sk*G0';
    if k>1, Sp=Ak*Sz(:,:,N_pred)*Ak'+Qa; Sj=Sp; else Sp=zeros(nz); Sj=zeros(nz); end
    Sz(:,:,1)=Sj; for j=2:N_pred, Sj=Ak*Sp*Ak'+Qa; Sp=Sj; Sz(:,:,j)=Sj; end
    if k<500, Nv=Ns; elseif k<=8000, Nv=Nm; else Nv=Nx; end
    H=zeros(Nv*m); f=zeros(Nv*m,1); ya=Pk*zk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
    for j=1:Nv
        rj=Rf(:,min(k+j,T_cl))-bk; rp=Rf(:,min(k+j-1,T_cl))-bk;
        if j==1, dr=rj-ya; Md=Mk{j}; Nd=Nk{j}(:,1:Nv*m); fc=-ya;
        else, dr=rj-rp; Md=Mk{j}-Mk{j-1}; Nd=Nk{j}(:,1:Nv*m)-Nk{j-1}(:,1:Nv*m); fc=zeros(p,1); end
        H=H+Nd'*Q*Nd; f=f+Nd'*Q*(Md*zk+fc-dr);
    end
    H=H+kron(eye(Nv),Rw); H=(H+H')/2+1e-8*eye(size(H));
    lb=u_min*ones(Nv*m,1); ub=u_max*ones(Nv*m,1);
    try [co,~,~]=quadprog(H,f,[],[],[],[],lb,ub,cp(1:Nv*m),oq); catch, co=max(min(-H\f,ub),lb); end
    uk=co(1:m); uC(:,k)=uk; costC(k)=0.5*co'*H*co+f'*co;
    xt=A*xt+B*uk+sw0*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0; cp(1:Nv*m)=max(min(cp(1:Nv*m),ub),lb);
    xk_prev=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && k>0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,beta_u,2,A,B,C,n,m,p);
        PbrS=null(Pk'); Ak=Akh;
        for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb_eps=zeros(ell,Nw); ec_eps=0; tl=k;
    end
end
eC=mean(abs(yC(2,500:end)-Rf(2,500:end))); vC=sum(yC(2,:)>y_max);

% ==== SM (SMPC: 已知Σ_ξ + 准确Σ_ε) ====
fprintf('SM ...\n'); yS=zeros(p,T_cl); uS=zeros(m,T_cl);
sS_eps=zeros(1,T_cl); sS_xi=zeros(1,T_cl); costS=zeros(1,T_cl);
tight_rec=zeros(ell,T_cl); bnd_rec=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0; PbrS=Pbar0;
cp=zeros(Nx*m,1); eb_eps=zeros(ell,Nw); ec_eps=0; eb_xi=zeros(p-ell,Nw); ec_xi=0;
Sxi=zeros(p-ell,p-ell); Sz=zeros(nz,nz,N_pred);
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc; xk_prev=zeros(ell,1);

for k=1:T_cl
    yk=C*xt+se0*randn(p,1); yS(:,k)=yk; xk=Rk'*yk;
    v_hat=Akh*xk_prev+Bk*uS(:,max(1,k-1));
    e_full=yk-Pk*v_hat; eps_k=Rk'*e_full; xi_k=PbrS'*e_full;
    bi=mod(k-1,Nw)+1; eb_eps(:,bi)=eps_k; ec_eps=min(ec_eps+1,Nw);
    if ec_eps>=Nw, Sk=(eb_eps*eb_eps')/Nw; Sk=(Sk+Sk')/2; end
    sS_eps(k)=sqrt(trace(Sk)/ell); sS_xi(k)=sqrt(trace(Sxi)/max(p-ell,1));
    zk=xk; Qa=G0*Sk*G0';
    if k>1, Sp=Ak*Sz(:,:,N_pred)*Ak'+Qa; Sj=Sp; else Sp=zeros(nz); Sj=zeros(nz); end
    Sz(:,:,1)=Sj; for j=2:N_pred, Sj=Ak*Sp*Ak'+Qa; Sp=Sj; Sz(:,:,j)=Sj; end
    if k<500, Nv=Ns; elseif k<=8000, Nv=Nm; else Nv=Nx; end
    H=zeros(Nv*m); f=zeros(Nv*m,1); ya=Pk*zk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
    for j=1:Nv
        rj=Rf(:,min(k+j,T_cl))-bk; rp=Rf(:,min(k+j-1,T_cl))-bk;
        if j==1, dr=rj-ya; Md=Mk{j}; Nd=Nk{j}(:,1:Nv*m); fc=-ya;
        else, dr=rj-rp; Md=Mk{j}-Mk{j-1}; Nd=Nk{j}(:,1:Nv*m)-Nk{j-1}(:,1:Nv*m); fc=zeros(p,1); end
        H=H+Nd'*Q*Nd; f=f+Nd'*Q*(Md*zk+fc-dr);
    end

    % ★ SMPC: Σ_y = P Σ_z P^T + Σ_v (已知噪声)
    n_cc=2*ell*Nv; A_ch=zeros(n_cc,Nv*m); b_ch=zeros(n_cc,1); ri=0;
    for jj=1:min(Nv,N_pred)
        N_fj=Nk{jj}(:,1:Nv*m); mu0j=Mk{jj}*zk;
        Sig_j=Pk*Sz(:,:,jj)*Pk'+Sigma_v;  % 已知静态噪声
        for ch=1:ell
            tv=norminv(0.84)*sqrt(max(Sig_j(ch,ch),1e-6));
            ri=ri+1; A_ch(ri,:)=N_fj(ch,:); b_ch(ri)=y_max-mu0j(ch)-tv;
            ri=ri+1; A_ch(ri,:)=-N_fj(ch,:); b_ch(ri)=mu0j(ch)-tv+3.0;
            if jj==1, tight_rec(ch,k)=tv; end
        end
    end
    A_ch=A_ch(1:ri,:); b_ch=b_ch(1:ri);
    Hf=H+kron(eye(Nv),Rw); Hf=(Hf+Hf')/2+1e-8*eye(size(Hf));
    lb=u_min*ones(Nv*m,1); ub=u_max*ones(Nv*m,1);
    try [co0,~,~]=quadprog(Hf,f,[],[],[],[],lb,ub,cp(1:Nv*m),oq); catch, co0=cp(1:Nv*m); end
    try [co,~,ef]=quadprog(Hf,f,A_ch,b_ch,[],[],lb,ub,co0,oq);
        bnd_rec(k)=double(ef==1 && any(abs(A_ch*co-b_ch)<1e-4));
        if ef<0, co=co0; end
    catch, co=co0; end
    uk=co(1:m); uS(:,k)=uk; costS(k)=0.5*co'*Hf*co+f'*co;
    xt=A*xt+B*uk+sw0*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0; cp(1:Nv*m)=max(min(cp(1:Nv*m),ub),lb);
    xk_prev=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && k>0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,beta_u,2,A,B,C,n,m,p);
        PbrS=null(Pk'); Ak=Akh;
        for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb_eps=zeros(ell,Nw); ec_eps=0; tl=k;
    end
end
eS=mean(abs(yS(2,500:end)-Rf(2,500:end))); vS=sum(yS(2,:)>y_max);

fprintf('C=%.3f(%d) SM=%.3f(%d) imp=%+.0f%% bnd=%d\n',eC,vC,eS,vS,(eC-eS)/eC*100,sum(bnd_rec));
fprintf('σ: C_eps=%.3f C_xi=%.3f S_eps=%.3f S_xi=%.3f\n',mean(sC_eps(500:end)),mean(sC_xi(500:end)),mean(sS_eps(500:end)),mean(sS_xi(500:end)));

% ==== 5行图 ====
figure('Position',[50,50,2000,1600],'Name','C vs SM');
subplot(5,1,1); hold on; grid on;
plot(1:T_cl,Rf(2,:),'k:',1:T_cl,yC(2,:),'b-',1:T_cl,yS(2,:),'r-','LineWidth',1);
yline(y_max,'g--'); for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
legend('ref',sprintf('C(%.3f,%d)',eC,vC),sprintf('SM(%.3f,%d)',eS,vS),'y_{max}','Location','best');
ylabel('y_2'); title(sprintf('y tracking: C=%.3f(%d) SM=%.3f(%d) bnd=%d',eC,vC,eS,vS,sum(bnd_rec)));

subplot(5,1,2); hold on; grid on;
plot(500:T_cl,sC_eps(500:end),'b-',500:T_cl,sS_eps(500:end),'r-','LineWidth',1);
yline(sw0,'k-'); for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
legend('C','SM','true sw','Location','best'); ylabel('sigma_eps');

subplot(5,1,3); hold on; grid on;
plot(500:T_cl,sC_xi(500:end),'b-',500:T_cl,sS_xi(500:end),'r-','LineWidth',1);
yline(se0,'k-'); for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
legend('C','SM','true se','Location','best'); ylabel('sigma_xi');

subplot(5,1,4); hold on; grid on;
plot(1:T_cl,uC(1,:),'b-',1:T_cl,uC(2,:),'c-',1:T_cl,uC(3,:),'g-','LineWidth',1);
for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end; ylabel('u'); title('C u'); ylim([-6 6]);

subplot(5,1,5); hold on; grid on;
semilogy(500:T_cl,max(costC(500:end),1e-10),'b-',500:T_cl,max(costS(500:end),1e-10),'r-','LineWidth',1);
for ss=2000:2000:T_cl, xline(ss,'Color',[.85 .85 .85]); end
legend('C','SM','Location','best'); ylabel('J_k (log)'); xlabel('k');

sgtitle(sprintf('SMPC: C vs SM  known Σ_ξ=%.2f  Rw=I',se0));
saveas(gcf,'final2_smpc.png');
fprintf('图: final2_smpc.png\n');

% 保存
sv.yC=yC; sv.yS=yS; sv.uC=uC; sv.uS=uS;
sv.sC_eps=sC_eps; sv.sC_xi=sC_xi; sv.sS_eps=sS_eps; sv.sS_xi=sS_xi;
sv.costC=costC; sv.costS=costS; sv.tight=tight_rec; sv.bnd=bnd_rec;
sv.Rf=Rf; sv.y_max=y_max; sv.eC=eC; sv.eS=eS; sv.vC=vC; sv.vS=vS;
save('final2_data.mat','-struct','sv');
fprintf('done: final2_data.mat\n');