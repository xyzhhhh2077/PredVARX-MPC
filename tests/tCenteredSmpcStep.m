classdef tCenteredSmpcStep < matlab.unittest.TestCase
    %tCenteredSmpcStep Centered prediction and Boole SMPC tightening.

    properties
        Paths
    end

    methods (TestClassSetup)
        function addCanonicalPaths(testCase)
            testCase.Paths = localpaths();
            addpath(testCase.Paths.copyP);
            testCase.addTeardown(@() rmpath(testCase.Paths.copyP));
        end
    end

    methods (Test)
        function centeredPredictionAndChanceTightening(testCase)
            % Arrange
            model.A = 0.8;
            model.B = 1;
            model.P = 1;
            model.R = 1;
            model.y_mean = 2;
            model.u_mean = 3;
            model.Sigma_eps = 0.04;
            model.Sigma_obs = 0.01;
            opt.N = 3;
            opt.Q = 10;
            opt.Ru = 0.1;
            opt.u_min = -10;
            opt.u_max = 10;
            opt.H = [1; -1];
            opt.h = [5; 5];
            opt.alpha_joint = 0.10;

            % Act
            [z_next, y_pred, U, diagOut] = centered_smpc_step(2, 2, model, opt);

            % Assert
            testCase.verifyLessThan(abs(z_next), 1e-12, ...
                'Centered state must subtract y_mean.');
            testCase.verifyLessThan( ...
                abs(y_pred(1) - (2 + U(1) - model.u_mean)), 1e-10, ...
                'Output prediction must restore y_mean and use centered input.');
            testCase.verifyGreaterThanOrEqual(min(U), opt.u_min - 1e-10);
            testCase.verifyLessThanOrEqual(max(U), opt.u_max + 1e-10);
            testCase.verifyGreaterThan(diagOut.z_quantile, 0, ...
                'Chance tightening requires positive quantile.');
            expectedRisk = opt.alpha_joint / (2 * size(opt.H, 1) * opt.N);
            testCase.verifyEqual(diagOut.risk_each, expectedRisk, 'AbsTol', 1e-12);
            testCase.verifyLessThan(max(diagOut.A_ch * U - diagOut.b_ch), 1e-7, ...
                'Returned plan must satisfy chance constraints.');
        end

        function alphaGivesPositiveQuantile(testCase)
            % Arrange: per-row risk must stay < 0.5 so norminv(1-risk) > 0.
            model.A = 0.5; model.B = 1; model.P = 1; model.R = 1;
            model.y_mean = 0; model.u_mean = 0;
            model.Sigma_eps = 0.01; model.Sigma_obs = 0.01;
            opt.N = 2; opt.Q = 1; opt.Ru = 0.1;
            opt.u_min = -5; opt.u_max = 5;
            opt.H = 1; opt.h = 10; opt.alpha_joint = 0.20;

            % Act
            [~, ~, ~, diagOut] = centered_smpc_step(0, 0, model, opt);

            % Assert
            testCase.verifyGreaterThan(diagOut.z_quantile, 0);
            testCase.verifyLessThan(diagOut.risk_each, 0.5);
        end
    end
end
