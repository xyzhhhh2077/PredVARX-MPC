function D=copyOPQSTR_generate_common_data(cfg)
% Generate one plant/offline dataset/reference/noise realization for all copies.
rng(cfg.seed,'twister'); C=zeros(cfg.p,cfg.n); C(1,1)=1; C(1,3)=0.16; C(2,2)=1; C(2,4)=-0.12;
for i=3:cfg.p, C(i,:)=0.45*randn(1,cfg.n); C(i,:)=C(i,:)/max(norm(C(i,:)),1e-12); end
u_off=1.20*randn(cfg.m,cfg.T_off); x_off=zeros(cfg.n,cfg.T_off+1); y_off=zeros(cfg.p,cfg.T_off);
for k=1:cfg.T_off, y_off(:,k)=C*x_off(:,k)+cfg.se*randn(cfg.p,1); x_off(:,k+1)=cfg.A*x_off(:,k)+cfg.B*u_off(:,k)+cfg.sw*randn(cfg.n,1); end
phase=2*pi*(0:cfg.T_cl-1)/cfg.noise_cycle;
swp=cfg.sw_min+(cfg.sw_max-cfg.sw_min)*0.5*(1-cos(phase));
sep=cfg.se_min+(cfg.se_max-cfg.se_min)*0.5*(1-cos(phase+cfg.noise_phase_e));
Wstd=randn(cfg.n,cfg.T_cl); Vstd=randn(cfg.p,cfg.T_cl);
Rf=zeros(cfg.p,cfg.T_cl); levels=[0.25 1.50 0.65 1.85 0.45;0.35 1.25 1.75 0.80 1.55]; seg=floor(cfg.T_cl/5);
for s=1:5, ix=(s-1)*seg+1:min(s*seg,cfg.T_cl); Rf(1,ix)=levels(1,s); Rf(2,ix)=levels(2,s); end
D=struct('C',C,'u_off',u_off,'x_off',x_off,'y_off',y_off,'sigma_w_profile',swp,'sigma_e_profile',sep,'Wstd',Wstd,'Vstd',Vstd,'Rf',Rf);
end
