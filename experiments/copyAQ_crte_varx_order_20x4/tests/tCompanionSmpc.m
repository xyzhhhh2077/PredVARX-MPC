classdef tCompanionSmpc < matlab.unittest.TestCase
    %TCOMPANIONSMPC Verifies the QP predicts with the companion state.

    methods (TestClassSetup)
        function addExperimentPath(testCase)
            sourceFolder = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(sourceFolder));
        end
    end

    methods (Test)
        function testPredictionUsesFullCompanionState(testCase)
            ABlocks = cat(3,[0.65 0.05;0 0.55],[-0.18 0;0.06 -0.10]);
            BBlocks = cat(3,[0.30;-0.12],[0.08;0.04]);
            P = eye(2);
            [Ac,Bc,Cc,Qc,meta] = buildArxCompanion( ...
                ABlocks,BBlocks,P,0.01*eye(2));
            model = struct('Ac',Ac,'Bc',Bc,'Cc',Cc,'Qc',Qc, ...
                'y_mean',[0.1;-0.2],'u_mean',0.15, ...
                'Sigma_obs',0.02*eye(2),'companion',meta);
            opt = struct('N',3,'Q',diag([4 3]),'Ru',0.2, ...
                'u_min',-1,'u_max',1,'H',eye(2),'h',10*ones(2,1), ...
                'alpha_joint',0.1);
            xi = [0.4;-0.3;-0.1;0.2;0.25];
            reference = [0.5;-0.4];

            [zNow,yPred,U,out] = companionSmpcStep(xi,reference,model,opt);
            vFirst = U(1)-model.u_mean;
            expectedXiNext = Ac*xi+Bc*vFirst;

            testCase.verifyEqual(zNow,xi(1:2),AbsTol=1e-12);
            testCase.verifyEqual(yPred,model.y_mean+Cc*expectedXiNext,AbsTol=1e-10);
            testCase.verifyEqual(out.prediction_state_dimension,5);
            testCase.verifyTrue(out.used_companion_prediction);
            testCase.verifyEqual(numel(out.lb),3);
            testCase.verifyGreaterThan(out.exitflag,0);
        end
    end
end
