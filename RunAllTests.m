import matlab.unittest.TestSuite;
import matlab.unittest.TestRunner;
import matlab.unittest.plugins.TAPPlugin;
import matlab.unittest.plugins.ToFile;
import('matlab.unittest.plugins.CodeCoveragePlugin');
import('matlab.unittest.plugins.codecoverage.CoberturaFormat');

if ispc
    pathsep = ';';
elseif isunix
    pathsep = ':';
else
   error ('Undefined path separator.');
end

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

    
    % need to add each tracked folder separately, apparently
    fprintf('\nRunAllTests.m: Adding folders to cover:\n')
    sep = strfind(p,pathsep); idx = 1; % idx "cursor" tracks position in path string
    for iF = 1:length(sep)
        this_folder = p(idx:sep(iF) - 1);
        disp(this_folder);
        
        coverageFile = fullfile(getenv('WORKSPACE'), sprintf('coverage%d.xml',iF));
        runner.addPlugin(CodeCoveragePlugin.forFolder(this_folder,'Producing', CoberturaFormat(coverageFile)));
        
        idx = sep(iF)+1; % update cursor to start of next path
    end

    % Run the tests
    results = runner.run(suite);
    display(results);
catch e
    fprintf('\n*********************\nRunAllTests.m failed!\n*********************\n');
    disp(getReport(e,'extended'));
    exit(1);
end
exit;
