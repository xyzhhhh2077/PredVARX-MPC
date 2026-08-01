classdef boptestSegmentContractTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addSource(testCase)
            src = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(src));
        end
    end

    methods (Test)
        function excludesTransitionsAcrossSegments(testCase)
            rng(20260801,'twister');
            p = 30; m = 6; q = 2; ell = 5; samplesPerSegment = 24;
            segmentId = repelem(1:3,samplesPerSegment);
            u = randn(m,numel(segmentId));
            y = randn(p,numel(segmentId));
            [Eanchor,~] = qr(randn(p,q),0);

            [A,B,P,R,S,stats] = fit_segmented_anchored_varx( ...
                y,u,segmentId,Eanchor,ell,struct('ridge',1e-8,'mu',0.10, ...
                'ntr_epsilon',10e-6));

            testCase.verifyEqual(stats.transition_count,3*(samplesPerSegment-1));
            testCase.verifyEqual(stats.segment_count,3);
            testCase.verifySize(A,[ell ell]);
            testCase.verifySize(B,[ell m]);
            testCase.verifySize(S,[ell ell]);
            testCase.verifyLessThan(norm(R'*P-eye(ell),'fro'),1e-7);
        end

        function alignsAdvanceControlWithReturnedMeasurement(testCase)
            % Synthetic BOPTEST contract: record k stores u_k and returned y_k.
            rng(8,'twister');
            n = 5000; u = randn(1,n); y = zeros(1,n);
            for k = 2:n
                y(k) = 0.4*y(k-1) + 1.7*u(k);
            end
            segmentId = ones(1,n);
            [A,B,~,~,~,~] = fit_segmented_anchored_varx( ...
                y,u,segmentId,1,1,struct('ridge',1e-12,'mu',0.10, ...
                'ntr_epsilon',10e-6));
            testCase.verifyEqual(A,0.4,'AbsTol',2e-3);
            testCase.verifyEqual(B,1.7,'AbsTol',2e-3);
        end
    end
end