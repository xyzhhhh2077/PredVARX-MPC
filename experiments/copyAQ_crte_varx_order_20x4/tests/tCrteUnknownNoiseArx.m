classdef tCrteUnknownNoiseArx < matlab.unittest.TestCase
    %TCRTEUNKNOWNNOISEARX Integration contract for CRTE with true ARX(s).

    methods (TestClassSetup)
        function addExperimentPath(testCase)
            sourceFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(sourceFolder));
        end
    end

    methods (Test)
        function testOrderTwoProfilesCandidateSpecificArx(testCase)
            rng(2207,'twister');
            p = 10;
            nu = 2;
            ell = 5;
            tracked = [1 2];
            sampleCount = 420;
            order = 2;
            F = diag([0.90 0.72 0.48]);
            G = [0.4 -0.1; 0.1 0.3; -0.05 0.2];
            C = randn(p,3);
            C(1,:) = [1 0 0];
            C(2,:) = [0 1 0];
            u = randn(nu,sampleCount);
            x = zeros(3,sampleCount+1);
            y = zeros(p,sampleCount);
            for k = 1:sampleCount
                y(:,k) = C*x(:,k)+0.04*randn(p,1);
                x(:,k+1) = F*x(:,k)+G*u(:,k)+0.03*randn(3,1);
            end
            opt = struct('mu_grid',[0 1], 'alpha',0.5, 'beta',0.5, ...
                'prediction_horizon',3, 'Ru',diag([0.2 0.8]), ...
                'num_random_subspaces',1, 'seed',77, 'reach_tau',1e-12);

            [ABlocks,BBlocks,P,R,SigmaEps,stats] = ...
                crteProfiledTeacherUnknownNoiseArx(y,u,ell,tracked,order,opt);

            testCase.verifySize(ABlocks,[ell ell order]);
            testCase.verifySize(BBlocks,[ell nu order]);
            testCase.verifySize(SigmaEps,[ell ell]);
            testCase.verifyLessThan(norm(R'*P-eye(ell),'fro'),1e-8);
            testCase.verifyEqual(stats.arx_order,order);
            testCase.verifyFalse(stats.uses_true_Sigma_n);
            testCase.verifySubstring(stats.noise_object,'cross-fitted');
            testCase.verifyEqual(stats.companion_state_dimension, ...
                order*ell+(order-1)*nu);
            testCase.verifyTrue(all([stats.rows.arx_order] == order));
            testCase.verifyGreaterThan(stats.num_feasible,0);
        end

        function testSignatureForbidsTrueSigmaN(testCase)
            sourceFolder = fileparts(fileparts(mfilename('fullpath')));
            source = fileread(fullfile(sourceFolder, ...
                'crteProfiledTeacherUnknownNoiseArx.m'));
            signature = regexp(source,'function[^\n]+','match','once');

            testCase.verifyFalse(contains(signature,'Sigma_n'));
            testCase.verifySubstring(source,'uses_true_Sigma_n=false');
        end
    end
end
