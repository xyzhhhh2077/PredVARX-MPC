function tests = test_opinion9_terminal
% TEST_OPINION9_TERMINAL Default OFF + optional terminal cost/set smoke tests.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root,'lib'));
testCase.TestData.root = root;
end

function [model, opt] = local_toy(seed)
rng(seed,'twister');
p=6; ell=3; m=2; N=4;
[U,~]=qr(randn(p,ell),0);
model.A=diag([0.9 0.7 0.5]);
model.B=0.25*randn(ell,m);
model.P=U; model.R=U;
model.y_mean=zeros(p,1); model.u_mean=zeros(m,1);
model.Sigma_eps=1e-3*eye(ell); model.Sigma_obs=1e-3*eye(p);
opt.N=N; opt.Q=eye(p); opt.Ru=0.1*eye(m);
opt.H=zeros(2,p); opt.H(1,1)=1; opt.H(2,2)=1; opt.h=[3;3];
opt.u_min=-3; opt.u_max=3; opt.alpha_joint=0.2;
end

function test_default_flags_off(testCase)
[model,opt]=local_toy(1);
y=0.1*randn(size(model.P,1),1); r=zeros(size(model.P,1),1);
[~,~,U,out]=centered_smpc_step(y,r,model,opt);
verifyTrue(testCase, out.exitflag>0);
verifyFalse(testCase, out.use_terminal_cost);
verifyFalse(testCase, logical(out.terminal_cost_applied));
verifyFalse(testCase, out.use_terminal_set);
verifyFalse(testCase, logical(out.terminal_set_applied));
verifyEqual(testCase, numel(U), opt.N*size(model.B,2));
end

function test_terminal_cost_on(testCase)
[model,opt]=local_toy(2);
opt.use_terminal_cost=true;
y=0.05*randn(size(model.P,1),1); r=zeros(size(model.P,1),1);
[~,~,~,out]=centered_smpc_step(y,r,model,opt);
verifyTrue(testCase, out.terminal_cost_applied);
verifyTrue(testCase, ~isempty(out.Pterm));
verifyTrue(testCase, min(eig(out.Pterm))>0);
end

function test_terminal_set_on_requires_alpha(testCase)
[model,opt]=local_toy(3);
opt.use_terminal_set=true;  % missing alpha_term
y=0.05*randn(size(model.P,1),1); r=zeros(size(model.P,1),1);
threw=false;
try
    centered_smpc_step(y,r,model,opt);
catch
    threw=true;
end
verifyTrue(testCase, threw);

opt.alpha_term=50; opt.use_terminal_cost=true;
[~,~,~,out]=centered_smpc_step(y,r,model,opt);
verifyTrue(testCase, out.terminal_set_applied);
verifyEqual(testCase, out.alpha_term, 50);
verifyTrue(testCase, isfinite(out.V_term));
end

function test_alpha_calibrate_and_rf_flags(testCase)
A=diag([0.9 0.8 0.6]);
Qf=eye(3); Pt=dlyap(A',Qf);
cal=calibrate_alpha_term(A,Pt,'Nmc',200,'seed',1);
verifyTrue(testCase, cal.alpha_recommend>0);
verifyFalse(testCase, cal.is_Xf_invariance_proof);
verifyTrue(testCase, cal.free_level_set_invariant_proxy);
end
