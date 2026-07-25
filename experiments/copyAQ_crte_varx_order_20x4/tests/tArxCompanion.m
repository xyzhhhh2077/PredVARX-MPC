classdef tArxCompanion < matlab.unittest.TestCase
    %TARXCOMPANION Contract tests for an ARX(s) companion realization.

    methods (TestClassSetup)
        function addExperimentPath(testCase)
            sourceFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(sourceFolder));
        end
    end

    methods (Test)
        function testOrderThreeOneStepPrediction(testCase)
            ABlocks = cat(3, [0.6 0.1; 0 0.5], ...
                [-0.2 0; 0.1 -0.1], [0.05 0; 0 0.04]);
            BBlocks = cat(3, [0.4; -0.1], [0.2; 0.05], [-0.08; 0.03]);
            P = [1 0; 0 1; 0.5 -0.2];
            SigmaEps = diag([0.04 0.09]);

            [Ac, Bc, Cc, Qc, meta] = buildArxCompanion( ...
                ABlocks, BBlocks, P, SigmaEps);
            z0 = [0.7; -0.3];
            z1 = [-0.2; 0.4];
            z2 = [0.1; 0.2];
            vPrev1 = 0.25;
            vPrev2 = -0.15;
            vNow = 0.6;
            xi = [z0; z1; z2; vPrev1; vPrev2];
            expectedZ = ABlocks(:,:,1)*z0 + ABlocks(:,:,2)*z1 + ...
                ABlocks(:,:,3)*z2 + BBlocks(:,:,1)*vNow + ...
                BBlocks(:,:,2)*vPrev1 + BBlocks(:,:,3)*vPrev2;

            xiNext = Ac*xi + Bc*vNow;
            testCase.verifyEqual(xiNext(1:2), expectedZ, AbsTol=1e-12);
            testCase.verifyEqual(Cc*xiNext, P*expectedZ, AbsTol=1e-12);
            testCase.verifySize(Ac, [8 8]);
            testCase.verifySize(Bc, [8 1]);
            testCase.verifyEqual(Qc(1:2,1:2), SigmaEps, AbsTol=1e-12);
            testCase.verifyEqual(Qc(3:end,:), zeros(6,8), AbsTol=1e-12);
            testCase.verifyEqual(meta.order, 3);
            testCase.verifyEqual(meta.state_dimension, 8);
        end

        function testOrderOneReducesToFirstOrder(testCase)
            A = [0.8 0.1; 0 0.6];
            B = [0.3; -0.2];
            P = eye(2);
            SigmaEps = 0.05*eye(2);

            [Ac, Bc, Cc, Qc, meta] = buildArxCompanion(A, B, P, SigmaEps);

            testCase.verifyEqual(Ac, A, AbsTol=1e-12);
            testCase.verifyEqual(Bc, B, AbsTol=1e-12);
            testCase.verifyEqual(Cc, P, AbsTol=1e-12);
            testCase.verifyEqual(Qc, SigmaEps, AbsTol=1e-12);
            testCase.verifyEqual(meta.order, 1);
        end
    end
end
