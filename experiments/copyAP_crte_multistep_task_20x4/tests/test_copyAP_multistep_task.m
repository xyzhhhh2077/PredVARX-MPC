function tests=test_copyAP_multistep_task
%TEST_COPYAP_MULTISTEP_TASK Focused task-stack, unknown-noise, and smoke tests.
tests=functiontests(localfunctions);
end

function testHorizonOneMatchesCopyAO(testCase)
[root,y,u,ell,tracked,base]=fixtureData(testCase);
ao=fullfile(fileparts(root),'copyAO_crte_teacher_profiled_unknown_noise'); addpath(ao);
[~,~,~,~,~,stAO]=crte_profiled_teacher_unknown_noise(y,u,ell,tracked,base);
opt=base; opt.task_horizon=1; opt.task_omega=1;
[~,~,~,~,~,st1]=crte_profiled_teacher_unknown_noise_multitask(y,u,ell,tracked,opt);
testCase.verifyLessThan(norm(st1.A_T-stAO.A_T,'fro')/max(norm(stAO.A_T,'fro'),1),1e-11);
testCase.verifyLessThan(norm(st1.B_T-stAO.B_T,'fro')/max(norm(stAO.B_T,'fro'),1),1e-11);
end

function testFutureTaskStackUsesAllHorizons(testCase)
[~,y,u,ell,tracked,base]=fixtureData(testCase);
opt=base; opt.task_horizon=3; opt.task_omega=ones(1,3)/3;
[~,~,~,~,~,st3]=crte_profiled_teacher_unknown_noise_multitask(y,u,ell,tracked,opt);
T=size(y,2); ntr=T-max(50,round(0.25*T)); idx=1:ntr-3;
ybar=mean(y(:,1:ntr),2); yc=y-ybar;
expected=[sqrt(1/3)*yc(tracked,idx+1);sqrt(1/3)*yc(tracked,idx+2);sqrt(1/3)*yc(tracked,idx+3)];
testCase.verifySize(st3.task_future_stack,[numel(tracked)*3,numel(idx)]);
testCase.verifyEqual(st3.task_future_stack,expected,'AbsTol',1e-11);
testCase.verifyEqual(st3.task_horizon,3);
testCase.verifyEqual(st3.task_future_rows,6);
end

function testTaskHorizonChangesOperator(testCase)
[~,y,u,ell,tracked,base]=fixtureData(testCase);
opt1=base; opt1.task_horizon=1; opt1.task_omega=1;
[~,~,~,~,~,st1]=crte_profiled_teacher_unknown_noise_multitask(y,u,ell,tracked,opt1);
opt3=base; opt3.task_horizon=3; opt3.task_omega=ones(1,3)/3;
[~,~,~,~,~,st3]=crte_profiled_teacher_unknown_noise_multitask(y,u,ell,tracked,opt3);
testCase.verifyGreaterThan(norm(st3.A_T-st1.A_T,'fro')/max(norm(st1.A_T,'fro'),1),1e-6);
end

function testUnknownNoiseContract(testCase)
[root,y,u,ell,tracked,base]=fixtureData(testCase);
src=fileread(fullfile(root,'crte_profiled_teacher_unknown_noise_multitask.m'));
sig=regexp(src,'function[^\n]+','match','once');
testCase.verifyFalse(contains(sig,'Sigma_n'));
opt=base; opt.task_horizon=3; opt.task_omega=ones(1,3)/3;
[~,~,~,~,~,st]=crte_profiled_teacher_unknown_noise_multitask(y,u,ell,tracked,opt);
testCase.verifyFalse(st.uses_true_Sigma_n);
testCase.verifyTrue(contains(st.noise_object,'cross-fitted'));
testCase.verifyEqual(st.noise_crossfit_scheme,'two-fold forward-chaining blocked split');
testCase.verifyGreaterThan(st.noise_crossfit_num_residuals,0);
end

function testSmokeLocksMuOne(testCase)
[root,~,~,~,~,~]=fixtureData(testCase);
smokeDir=fullfile(root,'results','smoke_test');
run=copyAP_run_seed(3,1,'Smoke',true,'Overwrite',true,'ResultsDir',smokeDir, ...
    'TOff',420,'TCl',60,'NumRandomSubspaces',0,'MuGrid',1,'PredictionHorizon',3);
testCase.verifyTrue(run.completed);
testCase.verifyEqual(run.config.task_horizon,3);
testCase.verifyEqual(run.seed.seed_id,1);
testCase.verifyEqual(run.teacher.selected_mu,1,'AbsTol',0);
testCase.verifyFalse(run.algorithm_contract.uses_true_Sigma_n);
testCase.verifyTrue(all(isfield(run.metrics,{'MAE','RMSE','Bias','qp_success_rate', ...
    'upper_violation_rate','abs_violation_rate','cost_mean','cost_sum'})));
testCase.verifyTrue(all(isfield(run.teacher,{'objective','prediction_term','task_term','noise_term'})));
testCase.verifyEqual(exist(run.result_file,'file'),2);
fid=H5F.open(run.result_file,'H5F_ACC_RDONLY','H5P_DEFAULT'); H5F.close(fid);
end

function [root,y,u,ell,tracked,base]=fixtureData(testCase)
here=fileparts(mfilename('fullpath')); root=fileparts(here);
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
originalRng=rng; testCase.addTeardown(@() rng(originalRng)); rng(2208,'twister');
p=10; m=2; ell=5; tracked=[1 2]; T=450;
F=diag([0.92 0.79 0.61]); Gp=[0.40 -0.12;0.15 0.30;-0.05 0.22];
C=randn(p,3); C(1,:)=[1 0 0]; C(2,:)=[0 1 0];
u=randn(m,T); x=zeros(3,T+1); y=zeros(p,T); D=diag(linspace(.03,.15,p));
for k=1:T
    y(:,k)=C*x(:,k)+D*randn(p,1);
    x(:,k+1)=F*x(:,k)+Gp*u(:,k)+.03*randn(3,1);
end
base=struct('mu_grid',1,'alpha',.5,'beta',.5,'prediction_horizon',3, ...
    'Ru',diag([.2 .8]),'num_random_subspaces',1,'seed',77,'reach_tau',1e-12);
end
