function results = run_all_tests()
%RUN_ALL_TESTS Run PredVARX-MPC canonical matlab.unittest suite.
here = fileparts(mfilename('fullpath'));
addpath(here);
import matlab.unittest.TestSuite
import matlab.unittest.TestRunner
import matlab.unittest.plugins.DiagnosticsOutputPlugin

suite = TestSuite.fromFolder(here, 'IncludingSubfolders', false);
runner = TestRunner.withTextOutput('Verbosity', 3);
runner.addPlugin(DiagnosticsOutputPlugin( ...
    'LoggingLevel', 3, 'IncludingPassingDiagnostics', false));
results = runner.run(suite);
disp(table(results));
assert(all([results.Passed]), 'One or more PredVARX-MPC unit tests failed.');
end
