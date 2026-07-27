classdef learnPreferredOutputDirectionsTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addSource(testCase)
            src=fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(src));
        end
    end
    methods (Test)
        function softPreferenceFavorsFrontWithoutHardLock(testCase)
            rng(20260730,'twister');
            p=10; m=2; T=800; u=randn(m,T); y=zeros(p,T); x=zeros(3,1);
            A=diag([0.86 0.72 0.55]); B=[.4 0;.1 .25;0 .12]; C=randn(p,3);
            for k=1:T
                y(:,k)=C*x+.03*randn(p,1); x=A*x+B*u(:,k)+.02*randn(3,1);
            end
            w=exp(-0.25*(0:p-1))';
            [E,st]=learn_preferred_output_directions(y,u,2,struct( ...
                'weights',w,'preference_strength',0.75,'reach_horizon',12,'Ru',eye(m)));
            testCase.verifySize(E,[p 2]);
            testCase.verifyEqual(E'*E,eye(2),AbsTol=1e-10);
            testCase.verifyGreaterThan(sum(st.contribution(1:3)),sum(st.contribution(end-2:end)));
            testCase.verifyGreaterThan(norm(E(3:end,:),'fro'),1e-3);
            testCase.verifyFalse(st.uses_new_training_data);
        end
    end
end
