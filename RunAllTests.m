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
    p = genpath(src);
    addpath(p); % this returns a single string with ; as separator between folders
    
    tests = fullfile(ws, 'tests');
    suite = testsuite(tests);
    
    runner = TestRunner.withTextOutput();
    
    % add TAP
    tapFile = fullfile(getenv('WORKSPACE'), 'testResults.tap');
    runner.addPlugin(TAPPlugin.producingOriginalFormat(ToFile(tapFile)));
    
    % Add Cobertura
    coverageFile = fullfile(getenv('WORKSPACE'), 'coverage.xml');
    
    % need to add each tracked folder separately, apparently
    sep = strfind(p,';'); idx = 1; % idx tracks position in path string
    for iF = 1:length(sep)
        this_folder = p(idx:sep(iF)-1);
        runner.addPlugin(CodeCoveragePlugin.forFolder(this_folder,'Producing', CoberturaFormat(coverageFile)));
        idx = sep(iF)+1;
    end

    % Run the tests
    results = runner.run(suite);
    display(results);
catch e
    disp(getReport(e,'extended'));
    exit(1);
end
exit;
