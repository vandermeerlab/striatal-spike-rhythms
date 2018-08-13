import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.TAPPlugin;
import matlab.unittest.plugins.ToFile;
import('matlab.unittest.plugins.CodeCoveragePlugin');
import('matlab.unittest.plugins.codecoverage.CoberturaFormat');

try
    % set up paths
    ws = getenv('WORKSPACE');
    src = fullfile(ws, 'shared');
    addpath(genpath(src));
    
    tests = fullfile(ws, 'tests');
    suite = testsuite(tests);
    
    runner = TestRunner.withTextOutput();
    
    % add TAP
    tapFile = fullfile(getenv('WORKSPACE'), 'testResults.tap');
    runner.addPlugin(TAPPlugin.producingOriginalFormat(ToFile(tapFile)));
    
    % Add Cobertura
    coverageFile = fullfile(getenv('WORKSPACE'), 'coverage.xml');
    runner.addPlugin(CodeCoveragePlugin.forFolder(src,'Producing', CoberturaFormat(coverageFile)));
    
    % Run the tests
    results = runner.run(suite);
    display(results);
catch e
    disp(getReport(e,'extended'));
    exit(1);
end
exit;
