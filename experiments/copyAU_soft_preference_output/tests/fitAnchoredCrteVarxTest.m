classdef fitAnchoredCrteVarxTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addSource(testCase)
            src = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(src));
        end
    end

    methods (Test)
        function usesFixedInteriorMuAndDimensionScaledNtr(testCase)
            rng(20260731,'twister');
            p = 9; m = 2; q = 2; ell = 5; T = 700;
            F = diag([0.88 0.72 0.55]);
            G = [0.35 -0.08; 0.10 0.27; -0.05 0.19];
            C = randn(p,3);
            u = randn(m,T); x = zeros(3,1); y = zeros(p,T);
            for k = 1:T
                y(:,k) = C*x + 0.04*randn(p,1);
                x = F*x + G*u(:,k) + 0.03*randn(3,1);
            end
            [Eanchor,~] = qr(randn(p,q),0);

            [A,B,P,R,S,st] = fit_anchored_varx(y,u,Eanchor,ell,struct( ...
                'ridge',1e-8,'mu',0.10,'ntr_epsilon',10e-6));

            testCase.verifyEqual(st.selected_mu,0.10,AbsTol=eps);
            testCase.verifyEqual(st.ntr_epsilon,10e-6,AbsTol=eps);
            testCase.verifyEqual(st.ntr_mode,'paper_trace_normalize_by_mean_trace');
            testCase.verifyEqual(st.Ntr_S_yu,local_ntr(st.S_yu,10e-6),AbsTol=1e-11);
            testCase.verifyEqual(st.Ntr_A_T,local_ntr(st.A_T,10e-6),AbsTol=1e-11);
            testCase.verifyEqual(st.Ntr_C_n,local_ntr(st.C_n,10e-6),AbsTol=1e-11);
            testCase.verifyLessThan(norm(R'*P-eye(ell),'fro'),1e-7);
            testCase.verifyLessThan(st.anchor_preservation_error,1e-7);
            testCase.verifySize(A,[ell ell]);
            testCase.verifySize(B,[ell m]);
            testCase.verifyGreaterThanOrEqual(min(eig((S+S')/2)),-1e-8);
        end
    end
end

function A = local_ntr(A,epsilon_ntr)
A = (A+A')/2;
d = size(A,1);
A = A/max(abs(trace(A))/d,epsilon_ntr);
end