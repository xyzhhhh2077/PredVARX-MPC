%% 副本 Zse: 信号/噪声分离后的 sigma_ebar
%  思路: ebar_low = EMA(ebar) → 信号
%        ebar_high = ebar - ebar_low → 噪声
%        sigma_ebar_clean = Dinkla(ebar_high)
clear; clc; addpath(fileparts(mfilename('fullpath')));

load('copyZ_data.mat');
beta_u=0.5;  % 未保存, 手动补

% 重跑Z2干净版, 加入ebar分离
n=6; m=3; p=15; ell=2; N_pred=30;
T_off=200; u_min=-5; u_max=5; Nr=50; Nw=50; Ns=5; Nm=12; Nx=30; gb=0.02;
ni=1000; nf_v=nf; sw0_v=sw0; se0_v=se0;
A=diag([.85,.70,.55,.40,.30,.20]); A(1,2)=.1; A(3,4)=.01;
Rw=2*eye(m); Q=zeros(p); Q(1:ell,1:ell)=eye(ell);
oq=optimset('Display','off','LargeScale','off');

rng(42,'twister'); B=randn(n,m); C=randn(p,n); [C,~]=qr(C,0);
Rf=zeros(p,T_cl);
for s=1:5, ks=(s-1)*2000+500; ke=min(s*2000+499,T_cl);
    if ks<=T_cl, Rf(1:ell,ks:ke)=rv(s); end; end

xo=zeros(n,T_off+1); yo=zeros(p,T_off); uo=10*randn(m,T_off);
for k=1:T_off, yo(:,k)=C*xo(:,k)+se0*randn(p,1); xo(:,k+1)=A*xo(:,k)+B*uo(:,k)+sw0*randn(n,1); end
[A0,B0,P0,~,R0,~,G0,Se0,~,F0,H0]=predvarx_identify(yo,uo,ell,beta_u,2,A,B,C,n,m,p);
Ac=F0; Pbar0=null(P0');
Ex=eye(ell); Mc=cell(1,N_pred); Nc=cell(1,N_pred);
for j=1:N_pred, Mc{j}=P0*Ex'*Ac^j; Nc{j}=zeros(p,N_pred*m);
    for i=0:j-1, Nc{j}(:,i*m+1:(i+1)*m)=P0*Ex'*Ac^(j-1-i)*H0; end; end

% ===== Z2_clean: ebar信号/噪声分离 =====
fprintf('Z2_clean ...\n');
alpha_sep=0.05;  % 分离的EMA系数(与bk一致)
yZc=zeros(p,T_cl); uZc=zeros(m,T_cl);
sZc_eps=zeros(1,T_cl); sZc_ebar=zeros(1,T_cl); sZc_ebar_raw=zeros(1,T_cl);
costZc=zeros(1,T_cl); softZc=zeros(1,T_cl);
ebar_low=zeros(p-ell,1);  % EMA跟踪信号

xt=zeros(n,1); Sk=Se0; Ak=Ac; Pk=P0; Rk=R0; Akh=A0; Bk=B0;
cp=zeros(Nx*m,1); eb=zeros(ell,Nw); ec=0; ebb_raw=zeros(p-ell,Nw); ebb_clean=zeros(p-ell,Nw); ecb=0; SebE_clean=zeros(p-ell,p-ell);
bk=zeros(p,1); yc=yo; uc=uo; tl=-Nr; Mk=Mc; Nk=Nc;
for k=1:T_cl
    sg=min(floor((k-1)/ni)+1,length(nf_v)); sw=sw0_v*nf_v(sg); se=se0_v*nf_v(sg);
    yk=C*xt+se*randn(p,1); yZc(:,k)=yk; xk=Rk'*yk;
    rc=Rf(:,min(k,T_cl)); bk=(1-gb)*bk+gb*(yk-rc);
    if k>=2, ek=xk-Akh*xk_prev-Bk*uZc(:,k-1);
        bi=mod(k-1,Nw)+1; eb(:,bi)=ek; ec=min(ec+1,Nw);
        if ec>=Nw, Sk=(eb*eb')/Nw; Sk=(Sk+Sk')/2; end; end
    sig=sqrt(trace(Sk)/ell); sZc_eps(k)=sig;
    
    % ★ 信号/噪声分离
    ebar_raw=Pbar0'*(yk-Pk*xk);  % p-ell维原始残差
    ebar_low=(1-alpha_sep)*ebar_low+alpha_sep*ebar_raw;  % EMA=信号
    ebar_high=ebar_raw-ebar_low;  % 残差=噪声
    
    bi2=mod(k-1,Nw)+1; ebb_clean(:,bi2)=ebar_high; ecb=min(ecb+1,Nw);
    if ecb>=Nw, SebE_clean=(ebb_clean*ebb_clean')/Nw; SebE_clean=(SebE_clean+SebE_clean')/2; end
    sZc_ebar(k)=sqrt(trace(SebE_clean)/(p-ell+1e-30));  % 干净的sigma_ebar
    
    % 非对称soft(和Z2一致)
    osp=yk(2)+sig-y_max; soft=max(0,min(osp,sig));
    softZc(k)=soft; zk=xk;
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
    uk=co(1:m); uZc(:,k)=uk; costZc(k)=0.5*co'*H*co+f'*co;
    xt=A*xt+B*uk+sw*randn(n,1);
    tmp=[co(m+1:end);co(end-m+1:end)]; cp(1:length(tmp))=tmp; cp(length(tmp)+1:end)=0; cp(1:Nv*m)=max(min(cp(1:Nv*m),ub),lb);
    xk_prev=xk; yc=[yc,yk]; uc=[uc,uk];
    if mod(k,Nr)==0 && k>0 && (k-tl)>=Nr
        [Akh,Bk,Pk,~,Rk,~,G0,~,~,~,Ha]=predvarx_identify(yc,uc,ell,beta_u,2,A,B,C,n,m,p);
        Ak=Akh; for j=1:N_pred, Mk{j}=Pk*Ex'*Ak^j;
            for i=0:j-1, Nk{j}(:,i*m+1:(i+1)*m)=Pk*Ex'*Ak^(j-1-i)*Ha; end; end
        eb=zeros(ell,Nw); ec=0; ebb_clean=zeros(p-ell,Nw); ecb=0; ebar_low=zeros(p-ell,1); tl=k;
    end
end
eZc=mean(abs(yZc(2,500:end)-Rf(2,500:end))); vZc=sum(yZc(2,:)>y_max);

% 对比: Z2原始 vs Z2_clean
fprintf('\\n=== ebar分离效果 ===\\n');
fprintf('Z2原始 s_ebar mean=%.4f max=%.4f\\n', mean(sZ2_ebar(500:end)), max(sZ2_ebar(500:end)));
fprintf('Z2clean s_ebar mean=%.4f max=%.4f\\n', mean(sZc_ebar(500:end)), max(sZc_ebar(500:end)));
fprintf('\\nZ2=%.3f(%d) Z2c=%.3f(%d)\\n', eZ2, violZ2, eZc_sav, vZc_sav);

% 追加保存到copyZ_data
load('copyZ_data.mat');
sv.sZc_eps=sZc_eps; sv.sZc_ebar=sZc_ebar; sv.yZc=yZc; sv.uZc=uZc; sv.eZc=eZc_sav; sv.vZc=vZc_sav;
save('copyZ_data.mat','-struct','sv');
fprintf('copyZ_data.mat 已更新 (新增 Z2_clean)\\n');

% 噪声对比图
figure('Position',[50,50,2000,600]);
subplot(1,3,1); hold on; grid on;
plot(500:T_cl,sZ2_ebar(500:end),'Color',[0.7,0,0.7],'LineWidth',1);
yline(se0,'k-'); ylabel('sigma_ebar'); title(sprintf('Z2原始 mean=%.2f',mean(sZ2_ebar(500:end))));
subplot(1,3,2); hold on; grid on;
plot(500:T_cl,sZc_ebar(500:end),'Color',[0,0.7,0.4],'LineWidth',1);
yline(se0,'k-'); ylabel('sigma_ebar'); title(sprintf('Z2clean mean=%.2f',mean(sZc_ebar(500:end))));
subplot(1,3,3); hold on; grid on;
plot(500:T_cl,sZ2_eps(500:end),'b-','LineWidth',1);
yline(sw0,'k-'); ylabel('sigma_eps'); title(sprintf('Z2 eps mean=%.2f',mean(sZ2_eps(500:end))));
saveas(gcf,'ebar_compare.png');
fprintf('图: ebar_compare.png\\n');