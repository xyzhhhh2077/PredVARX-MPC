classdef cdrClosedLoopContractTest < matlab.unittest.TestCase
    methods (Test)
        function testConfigHasAUControlSemantics(testCase)
            cfg = build_cdr_closed_loop_config();
            testCase.verifyEqual(cfg.p,30);
            testCase.verifyEqual(cfg.m,8);
            testCase.verifyEqual(cfg.q,2);
            testCase.verifyEqual(cfg.ell,5);
            testCase.verifyEqual(cfg.N,18);
            testCase.verifyLessThan(cfg.u_min,cfg.u_max);
            testCase.verifyGreaterThan(cfg.alpha_joint,0);
            testCase.verifyLessThan(cfg.alpha_joint,1);
            testCase.verifySize(cfg.Etask,[30 2]);
            testCase.verifySize(cfg.model_control.B,[5 8]);
            testCase.verifyGreaterThan(cfg.reference_amplitude,0);
            testCase.verifyGreaterThan(cfg.task_limit,cfg.reference_amplitude);
            testCase.verifyGreaterThanOrEqual(cfg.max_chance_tightening,0);
        end

        function testReferenceDesignUsesIdentifiedModel(testCase)
            cfg=build_cdr_closed_loop_config();
            Gss=cfg.H*cfg.model_control.P*((eye(cfg.ell)-cfg.model_control.A)\cfg.model_control.B);
            [basis,gains,~]=svd(Gss,'econ');
            expected_amplitude=min(0.20,0.35*max(diag(gains)));
            testCase.verifyEqual(cfg.reference_amplitude,expected_amplitude,'AbsTol',1e-14);
            testCase.verifyEqual(abs(cfg.reference_basis),abs(basis(:,1:cfg.q)),'AbsTol',1e-14);
        end

        function testClosedLoopRecordsRealControlActions(testCase)
            cfg = build_cdr_closed_loop_config();
            cfg.T = 12;
            cfg.warmup = 0;
            out = simulate_cdr_closed_loop(cfg);
            testCase.verifySize(out.u,[cfg.m,cfg.T]);
            testCase.verifySize(out.y,[cfg.p,cfg.T]);
            testCase.verifyEqual(numel(out.exitflag),cfg.T);
            testCase.verifyEqual(out.fallback,0);
            testCase.verifyGreaterThanOrEqual(min(out.u,[],'all'),cfg.u_min-1e-9);
            testCase.verifyLessThanOrEqual(max(out.u,[],'all'),cfg.u_max+1e-9);
            testCase.verifyLessThanOrEqual(max(out.maxcc,[],'omitnan'),1e-7);
        end
    end
end
