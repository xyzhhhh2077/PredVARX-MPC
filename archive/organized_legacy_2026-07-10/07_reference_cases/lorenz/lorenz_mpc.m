%% 副本 LORENZ — PredVAR原文仿真复现: Lorenz系统 + 正交P + 平稳噪声
%  Lorenz 3维DLV → P(15×3)正交 → y_k(15维) + 噪声 → PredVAR辨识 → MPC
clear; clc; addpath(fileparts(mfilename('fullpath')));

% ===== Lorenz 参数 =====
dt=0.01; sigma=10; rho=28; beta_l=8/3;
n_lorenz=3;  % 天然3维

% ===== 系统参数 =====
ell=3;  % ℓ = n = 3, 完美匹配原文
p=15; m=3;  % 观测+控制
N_total=10000; N_train=7000; N_test=3000;
sw0=0.1; se0=0.1;  % 平稳噪声
u_min=-5; u_max=5;

rng(42,'twister');

% ===== 1. 生成 Lorenz DLV + 控制 =====
v=zeros(n_lorenz,N_total);
v(:,1)=[1;1;1]*0.1;  % 初始点
u_seq=zeros(m,N_total);
B_ctrl=0.5*eye(n_lorenz,m);  % 控制输入矩阵 (3×3)
w_seq=sw0*randn(n_lorenz,N_total);

for k=1:N_total-1
    x=v(1,k); y=v(2,k); z=v(3,k);
    dx=sigma*(y-x);
    dy=x*(rho-z)-y;
    dz=x*y-beta_l*z;
    v(:,k+1)=v(:,k)+dt*[dx;dy;dz]+B_ctrl*u_seq(:,k)+w_seq(:,k);
end

% ===== 2. 正交P, P⊥ (15×3, 15×12) =====
[P_true,~]=qr(randn(p,ell),0);  % 正交列
P_perp=null(P_true');            % 正交补
e_seq=se0*randn(p,N_total);
y_all=P_true*v+P_perp*(se0*randn(p-ell,N_total))+e_seq;  % 含噪声观测
% 简化: y = P*v + total_noise
y_all=P_true*v+se0*randn(p,N_total);

% ===== 3. 离线辨识 (前7000步) =====
yo=y_all(:,1:N_train);
uo=u_seq(:,1:N_train);
[A0,B0,P0,~,R0,~,G0,Se0,~,F0,H0]=predvarx_identify(yo,uo,ell,0.5,2,...
    eye(n_lorenz),B_ctrl,eye(p),n_lorenz,m,p);
Ac=F0; Pbar0=null(P0');
fprintf('=== Lorenz辨识 ===\n');
fprintf('P̂相似度 vs P_true: 待验证\n');

% ===== 4. MPC 测试 (后3000步) =====
T_cl=N_test; N_pred=30;
Ns=5; Nm=12; Nx=30; Nr=50; Nw=50; gb=0.02; nz=ell;
Ex=eye(ell); Q=zeros(p); Q(1:ell,1:ell)=eye(ell); Rw=5*eye(m);  % 增大正则化, 平滑控制
y_max=12;  % 物理约束 (y范围~±15) (Lorenz输出范围~±20, 选一个合适值)

% 参考信号: 台阶 [0.2, 0.8, 0.2]
Rf=zeros(p,T_cl);
for s=1:3, ks=(s-1)*1000+300; ke=min(s*1000+299,T_cl);
    refs=[3,8,3];  % 参考信号 (在Lorenz范围内)
    if ks<=T_cl, Rf(1:ell,ks:ke)=refs(s); end
end

oq=optimset('Display','off','LargeScale','off','Algorithm','interior-point-convex');
Mc=cell(1,N_pred); Nc=cell(1,N_pred);
for j=1:N_pred, Mc{j}=P0*Ex'*Ac^j; Nc{j}=zeros(p,N_pred*m);
    for i=0:j-1, Nc{j}(:,i*m+1:(i+1)*m)=P0*Ex'*Ac^(j-1-i)*H0; end; end

% 初始化
yZ=zeros(p,T_cl); uZ=zeros(m,T_cl); costZ=zeros(1,T_cl);
sZ_eps=zeros(1,T_cl); sZ_xi=zeros(1,T_cl); softZ=zeros(1,T_cl);

v_test=v(:,N_train+1:end);  % 测试集真实DLV
u_test=u_seq(:,N_train+1:end);
w_test=w_seq(:,N_train+1:end);

xt=zeros(n_lorenz,1);  % 真实Lorenz状态
Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0; PbrS=Pbar0;
cp=zeros(Nx*m,1); eb_eps=zeros(ell,Nw); ec_eps=0; eb_xi=zeros(p-ell,Nw); ec_xi=0;
Sxi=zeros(p-ell,p-ell); Sz=zeros(nz,nz,N_pred);
bk=zeros(p,1); yc=yo(:,end-200:end); uc=uo(:,end-200:end); tl=-Nr; Mk=Mc; Nk=Nc;
xk_prev=zeros(ell,1);

for k=1:T_cl
    % 真实系统演化
    x=v_test(1,k); y=v_test(2,k); z_=v_test(3,k);
    dx=sigma*(y-x); dy=x*(rho-z_)-y; dz=x*y-beta_l*z_;
    xt=v_test(:,k);
    yk=y_all(:,N_train+k);  % 从预生成的y_all取

    yZ(:,k)=yk; xk=Rk'*yk;

    % PredVAR对齐
    v_hat=Akh*xk_prev+Bk*uZ(:,max(1,k-1));
    e_full=yk-Pk*v_hat;
    eps_k=Rk'*e_full;
    xi_k=PbrS'*e_full;

    bi_eps=mod(k-1,Nw)+1; eb_eps(:,bi_eps)=eps_k; ec_eps=min(ec_eps+1,Nw);
    if ec_eps>=Nw, Sk=(eb_eps*eb_eps')/Nw; Sk=(Sk+Sk')/2; end
    bi_xi=mod(k-1,Nw)+1; eb_xi(:,bi_xi)=xi_k; ec_xi=min(ec_xi+1,Nw);
    if ec_xi>=Nw, Sxi=(eb_xi*eb_xi')/Nw; Sxi=(Sxi+Sxi')/2; end
    sZ_eps(k)=sqrt(trace(Sk)/ell); sZ_xi(k)=sqrt(trace(Sxi)/max(p-ell,1));

    % Z2 soft
    sig=sqrt(trace(Sk)/ell);
    osp=yk(2)+sig-y_max; soft=max(0,min(osp,sig));
    softZ(k)=soft;

    % MPC
    zk=xk; Qa=G0*Sk*G0';
    if k>1, Sp=Ak*Sz(:,:,N_pred)*Ak'+Qa; Sj=Sp; else Sp=zeros(nz); Sj=zeros(nz); end
    Sz(:,:,1)=Sj; for j=2:N_pred, Sj=Ak*Sp*Ak'+Qa; Sp=Sj; Sz(:,:,j)=Sj; end
    if k<300, Nv=Ns; elseif k<=2400, Nv=Nm; else Nv=Nx; end
    H=zeros(Nv*m); f=zeros(Nv*m,1); ya=Pk*zk;

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
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0; cp(1:Nv*m)=max(min(cp(1:Nv*m),ub),lb);
    xk_prev=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && k>0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,0.5,2,...
            eye(n_lorenz),B_ctrl,eye(p),n_lorenz,m,p);
        PbrS=null(Pk'); Ak=Akh;
        for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb_eps=zeros(ell,Nw); ec_eps=0; eb_xi=zeros(p-ell,Nw); ec_xi=0; tl=k;
    end
end
eZ=mean(abs(yZ(2,300:end)-Rf(2,300:end))); vZ=sum(yZ(2,:)>y_max);

fprintf('LorenzMPC=%.3f(%d)  σ_eps=%.4f σ_xi=%.4f\n',eZ,vZ,mean(sZ_eps(300:end)),mean(sZ_xi(300:end)));

% 保存
sv.yZ=yZ; sv.uZ=uZ; sv.costZ=costZ; sv.softZ=softZ;
sv.sZ_eps=sZ_eps; sv.sZ_xi=sZ_xi; sv.Rf=Rf; sv.y_max=y_max;
sv.eZ=eZ; sv.vZ=vZ; sv.T_cl=T_cl;
save('lorenz_mpc.mat','-struct','sv');
fprintf('done: lorenz_mpc.mat\n');