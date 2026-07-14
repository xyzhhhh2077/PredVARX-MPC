%% Select regularized oblique strength on a held-out offline segment.
clear; clc; here=fileparts(mfilename('fullpath')); addpath(here);
rng(20260710,'twister');
n=6; m=3; p=30; ell=5; tracked=[1 2]; T=1500; sw=.045; se=.055;
A=diag([.94 .88 .78 .64 .50 .35]); A(1,2)=.10; A(2,3)=-.06; A(3,4)=.05; A(4,5)=.04;
B=[.34 -.10 .05;.12 .28 -.06;.05 .12 .24;-.05 .06 .18;.02 -.10 .14;.08 .02 -.08];
C=zeros(p,n); C(1,1)=1; C(1,3)=.16; C(2,2)=1; C(2,4)=-.12;
for i=3:p, C(i,:)=.45*randn(1,n); end
for i=3:p, C(i,:)=C(i,:)/max(norm(C(i,:)),1e-12); end
u=1.2*randn(m,T); x=zeros(n,T+1); y=zeros(p,T);
for k=1:T, y(:,k)=C*x(:,k)+se*randn(p,1); x(:,k+1)=A*x(:,k)+B*u(:,k)+sw*randn(n,1); end
Tfit=1100; alphas=[0 .02 .05 .10 .20 .40 .70 1.0];
fprintf('alpha dual asym recon latentRMSE outputRMSE trackedRMSE spectralRadius\n');
rows=zeros(numel(alphas),8);
for i=1:numel(alphas)
    a=alphas(i);
    [Ah,Bh,P,R,~,s]=control_aware_oblique_ivr_varx(y(:,1:Tfit),u(:,1:Tfit),ell,tracked,a);
    ym=mean(y(:,1:Tfit),2); um=mean(u(:,1:Tfit),2);
    z0=R'*(y(:,Tfit:T-1)-ym); z1=R'*(y(:,Tfit+1:T)-ym); ur=u(:,Tfit:T-1)-um;
    zp=Ah*z0+Bh*ur; yp=ym+P*zp;
    latent=sqrt(mean((z1-zp).^2,'all'));
    output=sqrt(mean((y(:,Tfit+1:T)-yp).^2,'all'));
    tracked_rmse=sqrt(mean((y(tracked,Tfit+1:T)-yp(tracked,:)).^2,'all'));
    rho=max(abs(eig(Ah)));
    rows(i,:)=[a s.dual_error s.pr_asymmetry s.reconstruction_residual latent output tracked_rmse rho];
    fprintf('%.2f %.2e %.3f %.3f %.4f %.4f %.4f %.4f\n',rows(i,:));
end
save(fullfile(here,'results','copyX_alpha_screen.mat'),'rows','alphas');
