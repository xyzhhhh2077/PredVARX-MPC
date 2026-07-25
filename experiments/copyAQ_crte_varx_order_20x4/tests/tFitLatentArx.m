classdef tFitLatentArx < matlab.unittest.TestCase
    %TFITLATENTARX Tests true multi-lag latent ARX identification.

    methods (TestClassSetup)
        function addExperimentPath(testCase)
            sourceFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(sourceFolder));
        end
    end

    methods (Test)
        function testRecoversOrderTwoBlocks(testCase)
            rng(41, 'twister');
            ell = 2;
            nu = 1;
            sampleCount = 1200;
            order = 2;
            ATrue = cat(3, [0.55 0.08; -0.03 0.45], ...
                [-0.12 0.02; 0.04 -0.08]);
            BTrue = cat(3, [0.35; -0.18], [0.09; 0.06]);
            u = randn(nu,sampleCount);
            z = zeros(ell,sampleCount);
            for k = order:sampleCount-1
                z(:,k+1) = ATrue(:,:,1)*z(:,k) + ...
                    ATrue(:,:,2)*z(:,k-1) + BTrue(:,:,1)*u(:,k) + ...
                    BTrue(:,:,2)*u(:,k-1);
            end

            [Ahat, Bhat, SigmaEps, stats] = fitLatentArx( ...
                z, u, order, 1e-12, order:sampleCount-1);

            testCase.verifyEqual(Ahat, ATrue, AbsTol=1e-8);
            testCase.verifyEqual(Bhat, BTrue, AbsTol=1e-8);
            testCase.verifyLessThan(norm(SigmaEps,'fro'), 1e-16);
            testCase.verifyEqual(stats.order, order);
            testCase.verifyEqual(stats.num_transitions, sampleCount-order);
            testCase.verifyEqual(stats.regressor_dimension, order*(ell+nu));
            testCase.verifyGreaterThanOrEqual(stats.design_rank, ...
                stats.regressor_dimension);
        end

        function testOrderChangesRegressorAndCompanionDimension(testCase)
            rng(9, 'twister');
            z = randn(3,80);
            u = randn(2,80);

            [Ahat, Bhat, ~, stats] = fitLatentArx(z,u,4,1e-8,4:79);
            [Ac, ~, ~, ~, meta] = buildArxCompanion( ...
                Ahat,Bhat,eye(3),eye(3));

            testCase.verifySize(Ahat,[3 3 4]);
            testCase.verifySize(Bhat,[3 2 4]);
            testCase.verifyEqual(stats.regressor_dimension,20);
            testCase.verifyEqual(meta.state_dimension,18);
            testCase.verifySize(Ac,[18 18]);
        end
    end
end
