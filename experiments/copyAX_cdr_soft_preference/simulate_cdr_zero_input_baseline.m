function out = simulate_cdr_zero_input_baseline(cfg,reference)
%SIMULATE_CDR_ZERO_INPUT_BASELINE Same CDR/noise/reference with physical u=0.
rng(cfg.seed,'twister');
T=cfg.T; p=cfg.p; q=cfg.q; x=cfg.x0;
y=zeros(p,T); s=zeros(q,T); u=zeros(cfg.m,T);
for k=1:T
    yk=cfg.C_plant*x+cfg.sensor_noise_std*randn(p,1);
    ystd=(yk-cfg.y_offset)./cfg.y_scale;
    y(:,k)=yk; s(:,k)=cfg.Etask'*ystd;
    x=cfg.A_plant*x;
end
start=min(max(cfg.warmup+1,1),T); warm=start:T;
err=s(:,warm)-reference(:,warm);
out=struct('y',y,'s',s,'u',u,'reference',reference, ...
    'MAE',mean(abs(err),2),'RMSE',sqrt(mean(err.^2,2)), ...
    'Bias',mean(err,2),'task_violation_rate',mean(abs(s)>cfg.task_limit,2));
end
