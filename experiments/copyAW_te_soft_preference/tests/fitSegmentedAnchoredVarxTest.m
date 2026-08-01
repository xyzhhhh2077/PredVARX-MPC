classdef fitSegmentedAnchoredVarxTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addSource(testCase)
            src = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(src));
        end
    end

    methods (Test)
        function excludesTransitionsAcrossRunBoundaries(testCase)
            rng(20260801,'twister');
            p = 41; m = 11; q = 2; ell = 5; samplesPerRun = 60;
            runId = repelem(1:4,samplesPerRun);
            u = randn(m,numel(runId));
            y = randn(p,numel(runId));
            [Eanchor,~] = qr(randn(p,q),0);

            [A,B,P,R,S,stats] = fit_segmented_anchored_varx( ...
                y,u,runId,Eanchor,ell,struct('ridge',1e-8,'mu',0.10, ...
                'ntr_epsilon',10e-6));

            testCase.verifyEqual(stats.transition_count,4*(samplesPerRun-1));
            testCase.verifyEqual(stats.segment_count,4);
            testCase.verifySize(A,[ell ell]);
            testCase.verifySize(B,[ell m]);
            testCase.verifySize(S,[ell ell]);
            testCase.verifyLessThan(norm(R'*P-eye(ell),'fro'),1e-7);
        end
    end
end