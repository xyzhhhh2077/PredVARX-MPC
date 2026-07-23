classdef tPredvarxIdentifyOblique < matlab.unittest.TestCase
    %tPredvarxIdentifyOblique Dual-basis invariants for oblique PredVARX.

    properties
        Paths
    end

    methods (TestClassSetup)
        function addCanonicalPaths(testCase)
            testCase.Paths = localpaths();
            addpath(testCase.Paths.copyO);
            testCase.addTeardown(@() rmpath(testCase.Paths.copyO));
        end
    end

    methods (Test)
        function dualBasisIdentitiesHold(testCase)
            % Arrange
            rng(7, 'twister');
            p = 6; ell = 2; m = 2; T = 120;
            A = diag([0.82, 0.61, 0.37]);
            B = randn(3, m);
            C = randn(p, 3); [C, ~] = qr(C, 0);
            u = randn(m, T);
            x = zeros(3, T + 1);
            y = zeros(p, T);
            for k = 1:T
                y(:, k) = C * x(:, k) + 0.05 * randn(p, 1);
                x(:, k + 1) = A * x(:, k) + B * u(:, k) + 0.03 * randn(3, 1);
            end

            % Act
            [~, ~, P, Pbar, R, Rbar, ~, Sigma_eps, Sigma_ebar] = ...
                predvarx_identify_oblique(y, u, ell, 0.5, 2, A, B, C, 3, m, p);

            % Assert
            testCase.verifySize(P, [p, ell]);
            testCase.verifySize(R, [p, ell]);
            testCase.verifyLessThan( ...
                norm(R' * P - eye(ell), 'fro'), 1e-8, ...
                'Oblique dual bases must satisfy R''P = I.');
            testCase.verifyLessThan( ...
                norm(R' * Pbar, 'fro'), 1e-8, ...
                'Static complement must be annihilated by R''.');
            testCase.verifyLessThan( ...
                norm(Rbar' * P, 'fro'), 1e-8, ...
                'Dynamic subspace must be annihilated by Rbar''.');
            testCase.verifyGreaterThanOrEqual( ...
                min(eig((Sigma_eps + Sigma_eps') / 2)), -1e-10);
            testCase.verifyGreaterThanOrEqual( ...
                min(eig((Sigma_ebar + Sigma_ebar') / 2)), -1e-10);
        end
    end
end
