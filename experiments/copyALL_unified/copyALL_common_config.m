function cfg=copyALL_common_config()
% One authoritative plant/controller fixture for H..Z unified copies.
cfg.seed=20260710; cfg.n=6; cfg.m=3; cfg.p=30; cfg.ell=5; cfg.tracked=[1 2];
cfg.T_off=1500; cfg.T_cl=1200; cfg.N=18; cfg.warm_start=151;
cfg.sw=0.045; cfg.se=0.055; cfg.noise_cycle=400;
cfg.sw_min=0.020; cfg.sw_max=0.090; cfg.se_min=0.025; cfg.se_max=0.100; cfg.noise_phase_e=pi/3;
cfg.u_min=-3; cfg.u_max=3; cfg.y_max=2.00; cfg.alpha_joint=0.10;
cfg.noise_window=40; cfg.reidentify_period=30; cfg.reidentify_recent=100;
% Historical H/K restored on the same fixture: 500/1200 is too sparse for a
% controlled comparison, so use one unified 300-step cadence/window.
cfg.hk_reidentify_period=300; cfg.hk_reidentify_window=300;
cfg.hk_dinkla_window=50; cfg.hk_bias_gamma=0.02;
cfg.oblique_alpha=0.02;
cfg.Q=zeros(cfg.p); cfg.Q(1,1)=80; cfg.Q(2,2)=80; cfg.Ru=0.18*eye(cfg.m);
cfg.A=diag([0.94 0.88 0.78 0.64 0.50 0.35]);
cfg.A(1,2)=0.10; cfg.A(2,3)=-0.06; cfg.A(3,4)=0.05; cfg.A(4,5)=0.04;
cfg.B=[0.34 -0.10 0.05;0.12 0.28 -0.06;0.05 0.12 0.24;-0.05 0.06 0.18;0.02 -0.10 0.14;0.08 0.02 -0.08];
end
