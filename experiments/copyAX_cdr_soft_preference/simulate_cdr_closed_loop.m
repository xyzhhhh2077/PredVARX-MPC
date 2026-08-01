function out = simulate_cdr_closed_loop(cfg)
%SIMULATE_CDR_CLOSED_LOOP Apply the AU SMPC QP to the true 200-state CDR plant.
rng(cfg.seed,'twister');
T=cfg.T; p=cfg.p; m=cfg.m; q=cfg.q;
x=cfg.x0; y=zeros(p,T); s=zeros(q,T); u=zeros(m,T);
reference=build_reference(cfg,T); exitflag=zeros(1,T);
maxcc=nan(1,T); cost=nan(1,T); fallback=0;
u_min_std=(-1-cfg.u_offset)./cfg.u_scale;
u_max_std=( 1-cfg.u_offset)./cfg.u_scale;
opt=struct('N',cfg.N,'Q',cfg.Q,'Ru',cfg.Ru,'u_min',u_min_std, ...
    'u_max',u_max_std,'H',cfg.H,'h',cfg.task_limit*ones(q,1), ...
    'alpha_joint',cfg.alpha_joint);
model=cfg.model_control;
Gram=cfg.Etask'*cfg.Etask;
for k=1:T
    yk=cfg.C_plant*x+cfg.sensor_noise_std*randn(p,1);
    ystd=(yk-cfg.y_offset)./cfg.y_scale;
    y(:,k)=yk; s(:,k)=cfg.Etask'*ystd;
    desired=reference(:,min(k+1,T));
    r_full=cfg.Etask*(Gram\desired);
    try
        [~,~,U,info]=centered_smpc_step(ystd,r_full,model,opt);
        ustd=U(1:m); exitflag(k)=info.exitflag;
        maxcc(k)=max(info.A_ch*U-info.b_ch); cost(k)=info.cost;
    catch ME
        fallback=fallback+1; exitflag(k)=-1;
        ustd=min(max(model.u_mean,u_min_std),u_max_std);
        if fallback<=3, fprintf('CDR fallback k=%d: %s\n',k,ME.message); end
    end
    uk=cfg.u_offset+cfg.u_scale.*ustd;
    u(:,k)=uk;
    x=cfg.A_plant*x+cfg.B_plant*uk;
end
start=min(max(cfg.warmup+1,1),T);
warm=start:T; err=s(:,warm)-reference(:,warm);
out=struct('y',y,'s',s,'u',u,'reference',reference, ...
    'exitflag',exitflag,'maxcc',maxcc,'cost',cost,'fallback',fallback, ...
    'MAE',mean(abs(err),2),'RMSE',sqrt(mean(err.^2,2)), ...
    'Bias',mean(err,2),'qp_success',mean(exitflag>0), ...
    'max_qp',max(maxcc,[],'omitnan'), ...
    'input_saturation_rate',mean(abs(u)>=(cfg.u_max-1e-8),'all'), ...
    'task_violation_rate',mean(abs(s)>cfg.task_limit,2));
end

function reference=build_reference(cfg,T)
q=cfg.q; reference=zeros(q,T); nseg=4;
levels=[0.00 0.55 -0.35 0.30; 0.00 -0.35 0.55 0.25];
levels=cfg.reference_amplitude*levels;
for j=1:nseg
    ix=floor((j-1)*T/nseg)+1:floor(j*T/nseg);
    reference(:,ix)=repmat(cfg.reference_basis*levels(:,j),1,numel(ix));
end
end
