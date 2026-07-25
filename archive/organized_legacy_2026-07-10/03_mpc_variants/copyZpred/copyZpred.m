%% 副本 Zpred — PredVAR原文对齐: e_k = P*eps_k + P⊥*xi_k
%  修复: 显式构造全创新e_k, σ_ξ用已知测量噪声
clear; clc; addpath(fileparts(mfilename('fullpath')));

n=6; m=3; p=15; ell=2; N_pred=30; T_cl=10000; T_off=200;
beta_u=0.5; sw0=0.1; se0=0.1; u_min=-5; u_max=5;
Nr=50; Nw=50; Ns=5; Nm=12; Nx=30; gb=0.02;
ni=1000; nf=[1,1.5,0.8,2,0.5,1.2,3,0.6,1,1.8];
rv=[1,3,0.2,2.5,0.5]; y_max=3.10;
Sigma_v=se0^2*eye(p);  % ★ 已知测量噪声 (PredVAR Σ_ξ的理论值)

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

% ===== Zpred: PredVAR对齐 + Z2 soft =====
fprintf('Zpred ...\n'); yZ=zeros(p,T_cl); uZ=zeros(m,T_cl);
sZ_eps=zeros(1,T_cl); sZ_xi=zeros(1,T_cl); sZ_e=zeros(1,T_cl);
costZ=zeros(1,T_cl); softZ=zeros(1,T_cl);

xt=zeros(n,1); Sk=Se0;  % Σ_ε (ℓ×ℓ)
Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0; PbrS=Pbar0;
cp=zeros(Nx*m,1);
eb_eps=zeros(ell,Nw); ec_eps=0;         % Dinkla for ε_k
eb_xi=zeros(p-ell,Nw); ec_xi=0;          % Dinkla for ξ_k
Sxi=zeros(p-ell,p-ell);                   % Σ_ξ
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc;
Sz=zeros(nz,nz,N_pred);
xk_prev=zeros(ell,1);  % 初始化

for k=1:T_cl
    sg=min(floor((k-1)/ni)+1,length(nf)); sw=sw0*nf(sg); se=se0*nf(sg);
    yk=C*xt+se*randn(p,1); yZ(:,k)=yk;
    xk=Rk'*yk;  % v_k = R y_k  (Eq.7)
    
    % === PredVAR 对齐: 构造全创新 e_k ===
    v_hat = Akh*xk_prev + Bk*uZ(:,max(1,k-1));  % v̂_k = Â v_{k-1} + B̂ u_{k-1}
    e_full = yk - Pk*v_hat;  % ★ e_k = y_k - P v̂_k  (Eq.4-5, p维)
    
    % 分解: e_k = P ε_k + P_⊥ ξ_k  (Eq.5)
    eps_k = Rk' * e_full;       % ε_k = R e_k  (Eq.7, ℓ维)
    xi_k = PbrS' * e_full;     % ξ_k = R_⊥ e_k  ((p-ℓ)维)
    % 注: PbrS'*Pk=0, PbrS'*yk = PbrS'*e_full, 所以 xi_k = PbrS'*yk 等价
    
    % Dinkla 窗口
    bi_eps=mod(k-1,Nw)+1; eb_eps(:,bi_eps)=eps_k; ec_eps=min(ec_eps+1,Nw);
    if ec_eps>=Nw, Sk=(eb_eps*eb_eps')/Nw; Sk=(Sk+Sk')/2; end  % Σ_ε
    
    bi_xi=mod(k-1,Nw)+1; eb_xi(:,bi_xi)=xi_k; ec_xi=min(ec_xi+1,Nw);
    if ec_xi>=Nw, Sxi=(eb_xi*eb_xi')/Nw; Sxi=(Sxi+Sxi')/2; end  % Σ_ξ (数据估)
    
    sZ_eps(k)=sqrt(trace(Sk)/ell);
    sZ_xi(k)=sqrt(trace(Sxi)/(p-ell+1e-30));
    sZ_e(k)=sqrt(trace(Pk*Sk*Pk'+PbrS*Sxi*PbrS')/p);  % ★ 全创新σ
    
    % Z2非对称soft
    sig=sqrt(trace(Sk)/ell);
    osp=yk(2)+sig-y_max; soft=max(0,min(osp,sig));
    softZ(k)=soft;
    
    % MPC
    zk=xk; Qa=G0*Sk*G0';
    if k>1, Sp=Ak*Sz(:,:,N_pred)*Ak'+Qa; Sj=Sp; else Sp=zeros(nz); Sj=zeros(nz); end
    Sz(:,:,1)=Sj; for j=2:N_pred, Sj=Ak*Sp*Ak'+Qa; Sp=Sj; Sz(:,:,j)=Sj; end
    if k<500, Nv=Ns; elseif k<=8000, Nv=Nm; else Nv=Nx; end
    H=zeros(Nv*m); f=zeros(Nv*m,1); ya=Pk*zk;
    
    % 代价 (Δ+绝对混合, 含bk)
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
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
    uk=co(1:m); uZ(:,k)=uk; costZ(k)=0.5*co'*H*co+f'*co;
    xt=A*xt+B*uk+sw*randn(n,1);
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
eZ=mean(abs(yZ(2,500:end)-Rf(2,500:end))); vZ=sum(yZ(2,:)>y_max);

fprintf('Zpred=%.3f(%d)  sigma_e=%.4f sigma_eps=%.4f sigma_xi=%.4f\n',...
    eZ,vZ,mean(sZ_e(500:end)),mean(sZ_eps(500:end)),mean(sZ_xi(500:end)));

% 保存
sv.yZ=yZ; sv.uZ=uZ; sv.sZ_eps=sZ_eps; sv.sZ_xi=sZ_xi; sv.sZ_e=sZ_e;
sv.costZ=costZ; sv.softZ=softZ; sv.Rf=Rf; sv.y_max=y_max; sv.T_cl=T_cl;
sv.eZ=eZ; sv.vZ=vZ; sv.nf=nf; sv.sw0=sw0; sv.se0=se0; sv.rv=rv;
save('copyZpred_data.mat','-struct','sv');
fprintf('done: copyZpred_data.mat\n');