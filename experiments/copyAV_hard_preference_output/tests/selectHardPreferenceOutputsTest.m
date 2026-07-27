classdef selectHardPreferenceOutputsTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addSource(testCase)
            src=fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(src));
        end
    end
    methods (Test)
        function locksHighestWeightOriginalAxes(testCase)
            w=[1 .8 .5 .3 .2]';
            [E,st]=select_hard_preference_outputs(w,2);
            testCase.verifyEqual(st.selected_indices,[1 2]);
            testCase.verifyEqual(E,[1 0;0 1;0 0;0 0;0 0]);
            testCase.verifyTrue(st.hard_locked);
            testCase.verifyFalse(st.uses_new_training_data);
        end
    end
end
