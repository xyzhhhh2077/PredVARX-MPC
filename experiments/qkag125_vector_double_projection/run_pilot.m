function run_pilot(nseeds)
% Vector double-projection analogy, not Chang et al.'s tensor algorithm.
if nargin < 1, nseeds = 5; end
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(here));
addpath(fullfile(repo, 'experiments', 'copyR_moqin_oblique'));
results = fullfile(here, 'results');
if ~exist(results,'dir'), mkdir(results); end
fid = fopen(fullfile(results, 'paired_metrics.csv'), 'w');
assert(fid > 0);
fprintf(fid, 'seed,method,subspace_error,prediction_rmse,dual_error,projector_error,train_trace\n');
for seed = 1:nseeds
  rand('state',10000+seed); randn('state',10000+seed);
  p=30; ell=2; Tid=1200; Ttest=500;
  A=[.92 .12;-.04 .78]; B=[.30;.08];
  C=zeros(p,ell);
  for j=1:p
    i=1+mod(j-1,ell);
    C(j,i)=.7+.3*rand();
    if rand()<.20, C(j,3-i)=.1*randn(); end
  end
  C=bsxfun(@rdivide,C,sqrt(sum(C.^2,1)));
  uid=randn(1,Tid); [yid,~]=sim(A,B,C,uid,.08);
  ute=.8*randn(1,Ttest); [yte,~]=sim(A,B,C,ute,.08);
  [~,~,P,Pbar,R,Rbar,~,~,~,~,~,~,info] = predvarx_identify_moqin(yid,uid,ell);
  yc=bsxfun(@minus,yid,mean(yid,2));
  Sigma=(yc*yc')/Tid; [U,D]=eig((Sigma+Sigma')/2);
  [d,ix]=sort(real(diag(D)),'descend'); U=real(U(:,ix)); d=max(d,1e-10);
  Ystar=bsxfun(@times,U'*yc,1./sqrt(d));
  % Map the original oblique readouts back into whitening coordinates.
  Q0=diag(sqrt(d))*U'*R;
  [Q0,~]=qr(Q0,0);
  for method=1:3
    if method==1
      Pm=P; Rm=R; Pbm=Pbar; Rbm=Rbar; label='MoQin_IVR';
    else
      Q=refine(Ystar,Q0,method==3,.08);
      Qbar=null(Q');
      Pm=U*diag(sqrt(d))*Q; Rm=U*diag(1./sqrt(d))*Q;
      Pbm=U*diag(sqrt(d))*Qbar; Rbm=U*diag(1./sqrt(d))*Qbar;
      if method==2, label='lag_projection_only'; else, label='double_projection'; end
    end
    assert(all(isfinite(Pm(:))) && all(isfinite(Rm(:))));
    dual=norm([Rm Rbm]'*[Pm Pbm]-eye(p),'fro');
    proj=norm((Pm*Rm')^2-Pm*Rm','fro');
    assert(dual<1e-5 && proj<1e-5);
    z=Rm'*yc; zm=mean(uid,2); uc=uid-zm;
    phi=[z(:,1:end-1);uc(:,1:end-1)];
    theta=(phi*phi'+1e-9*eye(ell+1))\(phi*z(:,2:end)');
    Ah=theta(1:ell,:)'; Bh=theta(ell+1,:)';
    ztest=Rm'*bsxfun(@minus,yte,mean(yid,2));
    zpred=Ah*ztest(:,1:end-1)+Bh*(ute(:,1:end-1)-zm);
    ypred=bsxfun(@plus,Pm*zpred,mean(yid,2));
    rmse=sqrt(mean((ypred(:)-yte(:,2:end)(:)).^2));
    [Cq,~]=qr(C,0); [Pq,~]=qr(Pm,0);
    sv=svd(Pq'*Cq); suberr=sqrt(max(0,ell-sum(min(sv,1).^2)));
    fprintf(fid,'%d,%s,%.12g,%.12g,%.12g,%.12g,%.12g\n',10000+seed,label,suberr,rmse,dual,proj,info.predicted_trace);
  end
  if mod(seed,5)==0, fprintf('completed %d/%d seeds\n',seed,nseeds); end
end
fclose(fid);
end

function Q=refine(Y,Q,second,sigma)
[p,T]=size(Y); N=T-1;
for iter=1:8
  z=Q'*Y; W=zeros(p,size(Q,2));
  for i=1:size(Q,2)
    target=z(i,1:N);
    if second
      others=setdiff(1:size(Q,2),i);
      X=z(others,2:T);
      target=target-(target*X')/(X*X'+1e-8*eye(numel(others)))*X;
    end
    w=Y(:,2:T)*target'/N;
    tau=1.5*sigma*sqrt(log(p)/N);
    w(abs(w)<tau)=0;
    if norm(w)<1e-10, w=Q(:,i); end
    W(:,i)=w/norm(w);
  end
  [Qnew,~]=qr(W,0);
  if rank(W)<size(Q,2), break; end
  if norm(Qnew*Qnew'-Q*Q','fro')<1e-7, Q=Qnew; break; end
  Q=Qnew;
end
end

function [y,x]=sim(A,B,C,u,sigma)
T=size(u,2); x=zeros(size(A,1),T); y=zeros(size(C,1),T);
for k=1:T
  if k>1, x(:,k)=A*x(:,k-1)+B*u(k-1)+.05*randn(size(A,1),1); end
  y(:,k)=C*x(:,k)+sigma*randn(size(C,1),1);
end
end
