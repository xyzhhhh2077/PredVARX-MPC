classdef validateCrteSourceTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addSource(testCase)
            src = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(src));
        end
    end

    methods (Test)
        function acceptsCorrectedCrteArtifact(testCase)
            D.model = struct('A',eye(3),'B',ones(3,1));
            D.stats = struct('selected_mu',0.10, ...
                'ntr_mode','paper_trace_normalize', ...
                'S_yu',diag([2 4 6]), ...
                'A_T',diag([3 6 9]), ...
                'C_n',diag([1 2 3]));
            [model,provenance] = validate_crte_source(D,0.10);

            testCase.verifyEqual(model,D.model);
            testCase.verifyEqual(provenance.selected_mu,0.10,AbsTol=eps);
            testCase.verifyEqual(provenance.ntr_epsilon,10e-6,AbsTol=eps);
            testCase.verifyEqual(provenance.ntr_mode, ...
                'paper_trace_normalize_by_mean_trace');
            testCase.verifyEqual(provenance.ntr_formula, ...
                'A/max(abs(trace(A))/d,10e-6)');
            testCase.verifyEqual(provenance.Ntr_A_T, ...
                D.stats.A_T/(abs(trace(D.stats.A_T))/3),AbsTol=1e-12);
        end

        function rejectsStaleMu(testCase)
            D.model = struct('A',eye(3),'B',ones(3,1));
            D.stats = struct('selected_mu',0.25, ...
                'ntr_mode','paper_trace_normalize', ...
                'S_yu',eye(3),'A_T',eye(3),'C_n',eye(3));
            testCase.verifyError(@()validate_crte_source(D,0.10), ...
                'copyAV:SourceMuMismatch');
        end

        function rejectsMissingRawNtrMatrices(testCase)
            D.model = struct('A',eye(3),'B',ones(3,1));
            D.stats = struct('selected_mu',0.10, ...
                'ntr_mode','paper_trace_normalize');
            testCase.verifyError(@()validate_crte_source(D,0.10), ...
                'copyAV:SourceNtrMismatch');
        end
    end
end