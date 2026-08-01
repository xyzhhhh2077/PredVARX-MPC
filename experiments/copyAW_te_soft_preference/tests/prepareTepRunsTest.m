classdef prepareTepRunsTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addSource(testCase)
            src = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(src));
        end
    end

    methods (Test)
        function separatesMeasurementsInputsAndRuns(testCase)
            raw = local_fixture();

            data = prepare_tep_runs(raw, 1:2, 3);

            testCase.verifySize(data.y_train, [41 8]);
            testCase.verifySize(data.u_train, [11 8]);
            testCase.verifySize(data.y_validation, [41 4]);
            testCase.verifyEqual(data.train_run_id, [ones(1,4), 2*ones(1,4)]);
            testCase.verifyEqual(data.validation_run_id, 3*ones(1,4));
            testCase.verifyEqual(data.y_train(1,:), [101:104, 201:204]);
            testCase.verifyEqual(data.u_train(1,:), [1101:1104, 1201:1204]);
            testCase.verifyFalse(any(data.y_train(1,:) == raw.faultNumber(1)));
        end

        function rejectsOverlappingRunSplit(testCase)
            raw = local_fixture();

            testCase.verifyError(@() prepare_tep_runs(raw, 1:2, 2:3), ...
                'prepare_tep_runs:OverlappingRuns');
        end
    end
end

function raw = local_fixture()
raw = table();
raw.faultNumber = zeros(12,1);
raw.simulationRun = repelem((1:3)',4);
raw.sample = repmat((1:4)',3,1);
for j = 1:41
    raw.(sprintf('xmeas_%d',j)) = 100*raw.simulationRun + raw.sample + j - 1;
end
for j = 1:11
    raw.(sprintf('xmv_%d',j)) = 1000 + 100*raw.simulationRun + raw.sample + j - 1;
end
end