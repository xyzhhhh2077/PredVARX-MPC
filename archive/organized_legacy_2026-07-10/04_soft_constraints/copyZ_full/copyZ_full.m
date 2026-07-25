%% 副本 Z — C + Z1(条件) + Z2(非对称) 3系统, 6噪声
clear; clc; addpath(fileparts(mfilename('fullpath')));

n=6; m=3; p=15; ell=2; N_pred=30; T_cl=10000; T_off=200;
beta_u=0.5; sw0=0.1; se0=0.1; u_min=-5; u_max=5;
Nr=50; Nw=50; Ns=5; Nm=12; Nx=30; gb=0.02;
ni=1000; nf=[1,1.5,0.8,2,0.5,1.2,3,0.6,1,1.8];
rv=[1,3,0.2,2.5,0.5]; y_max_cc=3.10;

rng(42,'twister');
A=diag([.85,.70,.55,.40,.30,.20]); A(1,2)=.1; A(3,4)=.01;
B=randn(n,m); C=randn(p,n); [C,~]=qr(C,0);

Rf=zeros(p,T_cl);
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

% ===== 副本 C (纯MPC) =====
fprintf('C ...\n'); yC=zeros(p,T_cl); uC=zeros(m,T_cl);
sC_eps=zeros(1,T_cl); sC_ebar=zeros(1,T_cl); costC=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0;
cp=zeros(Nx*m,1); eb=zeros(ell,Nw); ec=0; ebb=zeros(p-ell,Nw); ecb=0; SebE=zeros(p-ell,p-ell);
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc;
for k=1:T_cl
    sg=min(floor((k-1)/ni)+1,length(nf)); sw=sw0*nf(sg); se=se0*nf(sg);
    yk=C*xt+se*randn(p,1); yC(:,k)=yk; xk=Rk'*yk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
    if k>=2, ek=xk-Akh*xk_prev-Bk*uC(:,k-1);
        bi=mod(k-1,Nw)+1; eb(:,bi)=ek; ec=min(ec+1,Nw);
        if ec>=Nw, Sk=(eb*eb')/Nw; Sk=(Sk+Sk')/2; end; end
    sC_eps(k)=sqrt(trace(Sk)/ell);
    ebar=Pbar0'*(yk-Pk*xk); bi2=mod(k-1,Nw)+1; ebb(:,bi2)=ebar; ecb=min(ecb+1,Nw);
    if ecb>=Nw, SebE=(ebb*ebb')/Nw; SebE=(SebE+SebE')/2; end
    sC_ebar(k)=sqrt(trace(SebE)/(p-ell+1e-30));
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
    uk=co(1:m); uC(:,k)=uk; costC(k)=0.5*co'*H*co+f'*co;
    xt=A*xt+B*uk+sw*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0; cp(1:Nv*m)=max(min(cp(1:Nv*m),ub),lb);
    xk_prev=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && k>0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,beta_u,2,A,B,C,n,m,p);
        Ak=Akh; for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb=zeros(ell,Nw); ec=0; ebb=zeros(p-ell,Nw); ecb=0; tl=k;
    end
end
eC=mean(abs(yC(2,500:end)-Rf(2,500:end)));

% ===== Z1: 条件soft =====
fprintf('Z1 ...\n'); yZ1=zeros(p,T_cl); uZ1=zeros(m,T_cl);
sZ1_eps=zeros(1,T_cl); sZ1_ebar=zeros(1,T_cl); costZ1=zeros(1,T_cl); softZ1=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0;
cp=zeros(Nx*m,1); eb=zeros(ell,Nw); ec=0; ebb=zeros(p-ell,Nw); ecb=0; SebE=zeros(p-ell,p-ell);
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc;
for k=1:T_cl
    sg=min(floor((k-1)/ni)+1,length(nf)); sw=sw0*nf(sg); se=se0*nf(sg);
    yk=C*xt+se*randn(p,1); yZ1(:,k)=yk; xk=Rk'*yk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
    if k>=2, ek=xk-Akh*xk_prev-Bk*uZ1(:,k-1);
        bi=mod(k-1,Nw)+1; eb(:,bi)=ek; ec=min(ec+1,Nw);
        if ec>=Nw, Sk=(eb*eb')/Nw; Sk=(Sk+Sk')/2; end; end
    sig=sqrt(trace(Sk)/ell); sZ1_eps(k)=sig;
    ebar=Pbar0'*(yk-Pk*xk); bi2=mod(k-1,Nw)+1; ebb(:,bi2)=ebar; ecb=min(ecb+1,Nw);
    if ecb>=Nw, SebE=(ebb*ebb')/Nw; SebE=(SebE+SebE')/2; end
    sZ1_ebar(k)=sqrt(trace(SebE)/(p-ell+1e-30));
    % 条件soft
    margin=y_max_cc-yk(2);
    if margin<3*sig, soft=sig; else, soft=0; end
    softZ1(k)=soft; zk=xk;
    if k<500, Nv=Ns; elseif k<=8000, Nv=Nm; else Nv=Nx; end
    H=zeros(Nv*m); f=zeros(Nv*m,1); ya=Pk*zk;
    for j=1:Nv
        rj=Rf(:,min(k+j,T_cl))-bk; rp=Rf(:,min(k+j-1,T_cl))-bk;
        sft=zeros(p,1); sft(1:ell)=soft;
        if j==1, dr=rj-ya-sft; Md=Mk{j}; Nd=Nk{j}(:,1:Nv*m); fc=-ya;
        else, dr=rj-rp; Md=Mk{j}-Mk{j-1}; Nd=Nk{j}(:,1:Nv*m)-Nk{j-1}(:,1:Nv*m); fc=zeros(p,1); end
        H=H+Nd'*Q*Nd; f=f+Nd'*Q*(Md*zk+fc-dr);
    end
    H=H+kron(eye(Nv),Rw); H=(H+H')/2+1e-8*eye(size(H));
    lb=u_min*ones(Nv*m,1); ub=u_max*ones(Nv*m,1);
    try [co,~,~]=quadprog(H,f,[],[],[],[],lb,ub,cp(1:Nv*m),oq); catch, co=max(min(-H\f,ub),lb); end
    uk=co(1:m); uZ1(:,k)=uk; costZ1(k)=0.5*co'*H*co+f'*co;
    xt=A*xt+B*uk+sw*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0; cp(1:Nv*m)=max(min(cp(1:Nv*m),ub),lb);
    xk_prev=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && k>0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,beta_u,2,A,B,C,n,m,p);
        Ak=Akh; for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb=zeros(ell,Nw); ec=0; ebb=zeros(p-ell,Nw); ecb=0; tl=k;
    end
end
eZ1=mean(abs(yZ1(2,500:end)-Rf(2,500:end)));

% ===== Z2: 非对称soft =====
fprintf('Z2 ...\n'); yZ2=zeros(p,T_cl); uZ2=zeros(m,T_cl);
sZ2_eps=zeros(1,T_cl); sZ2_ebar=zeros(1,T_cl); costZ2=zeros(1,T_cl); softZ2=zeros(1,T_cl);
xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0;
cp=zeros(Nx*m,1); eb=zeros(ell,Nw); ec=0; ebb=zeros(p-ell,Nw); ecb=0; SebE=zeros(p-ell,p-ell);
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc;
for k=1:T_cl
    sg=min(floor((k-1)/ni)+1,length(nf)); sw=sw0*nf(sg); se=se0*nf(sg);
    yk=C*xt+se*randn(p,1); yZ2(:,k)=yk; xk=Rk'*yk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
    if k>=2, ek=xk-Akh*xk_prev-Bk*uZ2(:,k-1);
        bi=mod(k-1,Nw)+1; eb(:,bi)=ek; ec=min(ec+1,Nw);
        if ec>=Nw, Sk=(eb*eb')/Nw; Sk=(Sk+Sk')/2; end; end
    sig=sqrt(trace(Sk)/ell); sZ2_eps(k)=sig;
    ebar=Pbar0'*(yk-Pk*xk); bi2=mod(k-1,Nw)+1; ebb(:,bi2)=ebar; ecb=min(ecb+1,Nw);
    if ecb>=Nw, SebE=(ebb*ebb')/Nw; SebE=(SebE+SebE')/2; end
    sZ2_ebar(k)=sqrt(trace(SebE)/(p-ell+1e-30));
    % 非对称soft
    osp=yk(2)+sig-y_max_cc; soft=max(0,min(osp,sig));
    softZ2(k)=soft; zk=xk;
    if k<500, Nv=Ns; elseif k<=8000, Nv=Nm; else Nv=Nx; end
    H=zeros(Nv*m); f=zeros(Nv*m,1); ya=Pk*zk;
    for j=1:Nv
        rj=Rf(:,min(k+j,T_cl))-bk; rp=Rf(:,min(k+j-1,T_cl))-bk;
        sft=zeros(p,1); sft(1:ell)=soft;
        if j==1, dr=rj-ya-sft; Md=Mk{j}; Nd=Nk{j}(:,1:Nv*m); fc=-ya;
        else, dr=rj-rp; Md=Mk{j}-Mk{j-1}; Nd=Nk{j}(:,1:Nv*m)-Nk{j-1}(:,1:Nv*m); fc=zeros(p,1); end
        H=H+Nd'*Q*Nd; f=f+Nd'*Q*(Md*zk+fc-dr);
    end
    H=H+kron(eye(Nv),Rw); H=(H+H')/2+1e-8*eye(size(H));
    lb=u_min*ones(Nv*m,1); ub=u_max*ones(Nv*m,1);
    try [co,~,~]=quadprog(H,f,[],[],[],[],lb,ub,cp(1:Nv*m),oq); catch, co=max(min(-H\f,ub),lb); end
    uk=co(1:m); uZ2(:,k)=uk; costZ2(k)=0.5*co'*H*co+f'*co;
    xt=A*xt+B*uk+sw*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0; cp(1:Nv*m)=max(min(cp(1:Nv*m),ub),lb);
    xk_prev=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && k>0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,beta_u,2,A,B,C,n,m,p);
        Ak=Akh; for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb=zeros(ell,Nw); ec=0; ebb=zeros(p-ell,Nw); ecb=0; tl=k;
    end
end
eZ2=mean(abs(yZ2(2,500:end)-Rf(2,500:end)));

% 越界
vC=sum(yC(2,:)>y_max_cc); vZ1=sum(yZ1(2,:)>y_max_cc); vZ2=sum(yZ2(2,:)>y_max_cc);
fprintf('C=%.3f(%d) Z1=%.3f(%d) Z2=%.3f(%d)\n',eC,vC,eZ1,vZ1,eZ2,vZ2);

% ===== 保存 (6噪声+全部数据) =====
sv.yC=yC; sv.yZ1=yZ1; sv.yZ2=yZ2;
sv.uC=uC; sv.uZ1=uZ1; sv.uZ2=uZ2;
sv.sC_eps=sC_eps; sv.sC_ebar=sC_ebar;
sv.sZ1_eps=sZ1_eps; sv.sZ1_ebar=sZ1_ebar;
sv.sZ2_eps=sZ2_eps; sv.sZ2_ebar=sZ2_ebar;
sv.costC=costC; sv.costZ1=costZ1; sv.costZ2=costZ2;
sv.softZ1=softZ1; sv.softZ2=softZ2;
sv.Rf=Rf; sv.y_max=y_max_cc; sv.T_cl=T_cl;
sv.eC=eC; sv.eZ1=eZ1; sv.eZ2=eZ2;
sv.violC=vC; sv.violZ1=vZ1; sv.violZ2=vZ2;
sv.nf=nf; sv.ni=ni; sv.sw0=sw0; sv.se0=se0; sv.rv=rv;
sv.Rw_val=2;
save('copyZ_data.mat','-struct','sv');
fprintf('done: copyZ_data.mat (6噪声)\n');