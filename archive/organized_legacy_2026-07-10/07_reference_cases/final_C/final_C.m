%% 副本 FINAL — 线性系统 + PredVAR对齐 + 已知噪声 + 平稳
clear; clc; addpath(fileparts(mfilename('fullpath')));

n=6; m=3; p=15; ell=2; N_pred=30; T_cl=10000; T_off=200;
beta_u=0.5; sw0=0.1; se0=0.1; u_min=-5; u_max=5;
Nr=50; Nw=50; Ns=5; Nm=12; Nx=30; gb=0.02;
y_max=3.2;  Sigma_xi_known=se0^2;  % ★ 已知静态噪声

rng(42,'twister');
A=diag([.85,.70,.55,.40,.30,.20]); A(1,2)=.1; A(3,4)=.01;
B=randn(n,m); C=randn(p,n); [C,~]=qr(C,0);

Rf=zeros(p,T_cl); rv=[1,3,0.2,2.5,0.5];
for s=1:5, ks=(s-1)*2000+500; ke=min(s*2000+499,T_cl);
    if ks<=T_cl, Rf(1:ell,ks:ke)=rv(s); end; end

oq=optimset('Display','off','LargeScale','off');
Q=zeros(p); Q(1:ell,1:ell)=eye(ell); Rw=2*eye(m); nz=ell; Ex=eye(ell);

xo=zeros(n,T_off+1); yo=zeros(p,T_off); uo=10*randn(m,T_off);
for k=1:T_off, yo(:,k)=C*xo(:,k)+se0*randn(p,1); xo(:,k+1)=A*xo(:,k)+B*uo(:,k)+sw0*randn(n,1); end
[A0,B0,P0,~,R0,~,G0,Se0,~,F0,H0]=predvarx_identify(yo,uo,ell,beta_u,2,A,B,C,n,m,p);
Ac=F0; Pbar0=null(P0');
Mc=cell(1,N_pred); Nc=cell(1,N_pred);
for j=1:N_pred, Mc{j}=P0*Ex'*Ac^j; Nc{j}=zeros(p,N_pred*m);
    for i=0:j-1, Nc{j}(:,i*m+1:(i+1)*m)=P0*Ex'*Ac^(j-1-i)*H0; end; end

% ===== 副本 C (纯MPC, 平稳噪声, PredVAR对齐) =====
fprintf('C ...\n'); yC=zeros(p,T_cl); uC=zeros(m,T_cl);
sC_eps=zeros(1,T_cl); sC_xi=zeros(1,T_cl); costC=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0; PbrS=Pbar0;
cp=zeros(Nx*m,1); eb_eps=zeros(ell,Nw); ec_eps=0; eb_xi=zeros(p-ell,Nw); ec_xi=0;
Sxi=zeros(p-ell,p-ell); Sz=zeros(nz,nz,N_pred);
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc; xk_prev=zeros(ell,1);

for k=1:T_cl
    yk=C*xt+se0*randn(p,1); yC(:,k)=yk; xk=Rk'*yk;
    v_hat=Akh*xk_prev+Bk*uC(:,max(1,k-1));
    e_full=yk-Pk*v_hat;
    eps_k=Rk'*e_full; xi_k=PbrS'*e_full;
    bi=mod(k-1,Nw)+1; eb_eps(:,bi)=eps_k; ec_eps=min(ec_eps+1,Nw);
    if ec_eps>=Nw, Sk=(eb_eps*eb_eps')/Nw; Sk=(Sk+Sk')/2; end
    bi2=mod(k-1,Nw)+1; eb_xi(:,bi2)=xi_k; ec_xi=min(ec_xi+1,Nw);
    if ec_xi>=Nw, Sxi=(eb_xi*eb_xi')/Nw; Sxi=(Sxi+Sxi')/2; end
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
        eb_eps=zeros(ell,Nw); ec_eps=0; eb_xi=zeros(p-ell,Nw); ec_xi=0; tl=k;
    end
end
eC=mean(abs(yC(2,500:end)-Rf(2,500:end))); vC=sum(yC(2,:)>y_max);

fprintf('C=%.3f(%d)  σ_eps=%.4f σ_xi=%.4f (true=%.1f)\n',eC,vC,mean(sC_eps(500:end)),mean(sC_xi(500:end)),se0);

sv.yC=yC; sv.uC=uC; sv.sC_eps=sC_eps; sv.sC_xi=sC_xi; sv.costC=costC;
sv.Rf=Rf; sv.eC=eC; sv.vC=vC; sv.y_max=y_max;
save('final_C.mat','-struct','sv');
fprintf('done: final_C.mat\n');