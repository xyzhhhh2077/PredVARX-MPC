%% 副本 Q 数据采集 — 只跑核心循环 + 保存数据
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
[A0,B0,P0,~,R0,~,G0,Se0,~,F0,H0]=predvarx_identify(yo,uo,ell,beta_u,2,A,B,C,n,m,p);
Ac=F0; Pbr=null(P0');
Mc=cell(1,N_pred); Nc=cell(1,N_pred);
for j=1:N_pred, Mc{j}=P0*Ex'*Ac^j; Nc{j}=zeros(p,N_pred*m);
    for i=0:j-1, Nc{j}(:,i*m+1:(i+1)*m)=P0*Ex'*Ac^(j-1-i)*H0; end; end

fprintf('副本Q ...\n');
yQ=zeros(p,T_cl); uQ=zeros(m,T_cl); sQe=zeros(1,T_cl); sQeb=zeros(1,T_cl); costQ=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0; PbrS=Pbr;
cp=zeros(Nx*m,1);
eb=zeros(ell,Nw); ec=0; ebb=zeros(p-ell,Nw); ecb=0; SebE=zeros(p-ell,p-ell);
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc;
for k=1:T_cl
    sg=min(floor((k-1)/ni)+1,length(nf)); sw=sw0*nf(sg); se=se0*nf(sg);
    yk=C*xt+se*randn(p,1); yQ(:,k)=yk; xk=Rk'*yk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
    if k>=2, ek=xk-Akh*xk_prev-Bk*uQ(:,k-1);
        bi=mod(k-1,Nw)+1; eb(:,bi)=ek; ec=min(ec+1,Nw);
        if ec>=Nw, Sk=(eb*eb')/Nw; Sk=(Sk+Sk')/2; end; end
    sQe(k)=sqrt(trace(Sk)/ell);
    ebar=PbrS'*(yk-Pk*xk); bi2=mod(k-1,Nw)+1; ebb(:,bi2)=ebar; ecb=min(ecb+1,Nw);
    if ecb>=Nw, SebE=(ebb*ebb')/Nw; SebE=(SebE+SebE')/2; end
    sQeb(k)=sqrt(trace(SebE)/(p-ell+1e-30));
    zk=xk; if k<500, Nv=Ns; elseif k<=8000, Nv=Nm; else Nv=Nx; end
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
    uk=co(1:m); uQ(:,k)=uk; costQ(k)=0.5*co'*H*co+f'*co;
    xt=A*xt+B*uk+sw*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0; cp(1:Nv*m)=max(min(cp(1:Nv*m),ub),lb);
    xk_prev=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && k>0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,beta_u,2,A,B,C,n,m,p);
        PbrS=null(Pk'); Ak=Akh; for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb=zeros(ell,Nw); ec=0; ebb=zeros(p-ell,Nw); ecb=0; tl=k;
    end
end
eQ=mean(abs(yQ(2,500:end)-Rf(2,500:end)));

% 越界检查
y_max_val=3.2; viol=sum(yQ(2,:)>y_max_val);
fprintf('eQ=%.3f max_y=%.4f viol=%d/%d\n', eQ, max(yQ(2,:)), viol, T_cl);

% 保存
sv.yQ=yQ; sv.uQ=uQ; sv.sQe=sQe; sv.sQeb=sQeb; sv.costQ=costQ;
sv.Rf=Rf; sv.eQ=eQ; sv.T_cl=T_cl; sv.y_max=y_max_val;
sv.nf=nf; sv.ni=ni; sv.sw0=sw0; sv.rv=rv;
save('copyQ_data.mat','-struct','sv');
fprintf('done: copyQ_data.mat\n');