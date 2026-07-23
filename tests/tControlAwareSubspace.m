classdef tControlAwareSubspace < matlab.unittest.TestCase
    %tControlAwareSubspace Tracked-output coverage for reduced orthogonal model.

    properties
        Paths
    end

    methods (TestClassSetup)
        function addCanonicalPaths(testCase)
            testCase.Paths = localpaths();
            addpath(testCase.Paths.copyQ);
            testCase.addTeardown(@() rmpath(testCase.Paths.copyQ));
        end
    end

    methods (Test)
        function trackedAxesLieInSubspace(testCase)
            % Arrange
            rng(11, 'twister');
            y = randn(8, 300);
            u = randn(2, 300);
            tracked = [1 2];
            ell = 4;

            % Act
            [Ahat, Bhat, P, R, Sigma_eps, stats] = ...
                control_aware_subspace_varx(y, u, ell, tracked);

            % Assert
            E = zeros(8, 2); E(tracked, :) = eye(2);
            testCase.verifySize(P, [8, ell]);
            testCase.verifySize(Ahat, [ell, ell]);
            testCase.verifySize(Bhat, [ell, 2]);
            testCase.verifyLessThan(norm(P' * P - eye(ell), 'fro'), 1e-10);
            testCase.verifyLessThan(norm(P * P' * E - E, 'fro'), 1e-10, ...
                'Tracked output axes must lie in the retained subspace.');
            testCase.verifyLessThan(norm(R - P, 'fro'), 1e-12);
            testCase.verifyGreaterThanOrEqual(min(eig((Sigma_eps + Sigma_eps') / 2)), -1e-10);
            testCase.verifyLessThan(stats.tracked_projection_error, 1e-10);
        end

        function ellBelowTrackedCountErrors(testCase)
            y = randn(5, 50);
            u = randn(1, 50);
            testCase.verifyError( ...
                @() control_aware_subspace_varx(y, u, 1, [1 2]), ...
                'control_aware_subspace_varx:ellTooSmall');
        end
    end
end
