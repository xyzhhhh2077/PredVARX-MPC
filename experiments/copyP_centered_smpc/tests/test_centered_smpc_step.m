function test_centered_smpc_step
% Required behavior: prediction uses identified y/u means and a positive
% Boole-allocated Gaussian tightening.
A = 0.8; B = 1; P = 1; R = 1;
model.A=A; model.B=B; model.P=P; model.R=R;
model.y_mean=2; model.u_mean=3; model.Sigma_eps=0.04; model.Sigma_obs=0.01;
opt.N=3; opt.Q=10; opt.Ru=0.1; opt.u_min=-10; opt.u_max=10;
opt.H=[1;-1]; opt.h=[5;5]; opt.alpha_joint=0.10;

[z_next, y_pred, U, diag] = centered_smpc_step(2, 2, model, opt);
assert(abs(z_next) < 1e-12, 'Centered state must subtract y_mean.');
assert(abs(y_pred(1) - (2 + U(1)-model.u_mean)) < 1e-10, ...
    'Output prediction must restore y_mean and use centered input.');
assert(all(U >= opt.u_min-1e-10 & U <= opt.u_max+1e-10));
assert(diag.z_quantile > 0, 'Chance tightening requires positive quantile.');
assert(abs(diag.risk_each - opt.alpha_joint/(2*size(opt.H,1)*opt.N)) < 1e-12);
assert(max(diag.A_ch*U-diag.b_ch) < 1e-7, 'Returned plan must satisfy chance constraints.');
fprintf('PASS: centered prediction and Boole SMPC tightening hold.\n');
end
