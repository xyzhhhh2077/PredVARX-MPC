%% 快速验证 predvarx_identify 接口兼容性
clear; clc; addpath(fileparts(mfilename('fullpath')));

n=6; m=3; p=15; ell=2; beta_u=0.5; L_u=2;

rng(42,'twister');
A=diag([.85,.70,.55,.40,.30,.20]); A(1,2)=.1; A(3,4)=.01;
B=randn(n,m); C=randn(p,n); [C,~]=qr(C,0);
sw0=0.1; se0=0.1;

% 离线数据
T_off=200;
xo=zeros(n,T_off+1); yo=zeros(p,T_off); uo=10*randn(m,T_off);
for k=1:T_off, yo(:,k)=C*xo(:,k)+se0*randn(p,1); xo(:,k+1)=A*xo(:,k)+B*uo(:,k)+sw0*randn(n,1); end

% 调用
tic;
[A0,B0,P0,~,R0,~,G0,Se0,Sebar0,F0,H0,~,~]=predvarx_identify(yo,uo,ell,beta_u,2,A,B,C,n,m,p);
t_ivr=toc;

% 维度检查
fprintf('\n=== 维度检查 ===\n');
fprintf('A0: %dx%d (期望 %dx%d)\n', size(A0), ell, ell);
fprintf('B0: %dx%d (期望 %dx%d)\n', size(B0), ell, m);
fprintf('P0: %dx%d (期望 %dx%d)\n', size(P0), p, ell);
fprintf('R0: %dx%d (期望 %dx%d)\n', size(R0), p, ell);  % 注意: R0 = P0 (QR正交化后)
fprintf('F0: %dx%d (期望 %dx%d)\n', size(F0), ell, ell);
fprintf('H0: %dx%d (期望 %dx%d)\n', size(H0), ell, m);
fprintf('Se0: %dx%d (期望 %dx%d)\n', size(Se0), ell, ell);
fprintf('Sebar0: %dx%d (期望 %dx%d)\n', size(Sebar0), p-ell, p-ell);

% 正交性检查
fprintf('\n=== 正交性 ===\n');
fprintf('P0''*P0 - I: %.2e\n', norm(P0'*P0 - eye(ell), 'fro'));
fprintf('R0''*P0 - I: %.2e\n', norm(R0'*P0 - eye(ell), 'fro'));

% IVR 子空间质量
C_true = C;
cos_theta = svd(P0' * C_true);
fprintf('cos(theta): [%.4f, %.4f]\n', cos_theta(1), cos_theta(2));

% 特征值
fprintf('A0 eig: [%.4f, %.4f]\n', sort(real(eig(A0)), 'descend'));

% Sigma 检查
fprintf('Sigma_eps diag: [%.4f, %.4f]\n', diag(Se0));
fprintf('Sigma_ebar diag norm: %.4f\n', norm(diag(Sebar0)));

fprintf('\nIVR time: %.3f sec\n', t_ivr);
fprintf('=== 接口验证完成 ===\n');
