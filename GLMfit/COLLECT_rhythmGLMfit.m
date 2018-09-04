% COLLECT_rhythmGLMfit.m
%
% Collector script for output of ALL_rhythmGLMfit

cfg.inputDir = 'C:\temp'; % where the files to load are
cfg.input_prefix = 'R0_'; 
cfg.models = {'allphi', 'all'}; % models to be compared to baseline
cfg.nMaxCells = 100;
cfg.nTimeBins = 100;
cfg.nSpaceBins = 100;

%%
pushdir(cfg.inputDir);

fd = FindFiles(cat(2,cfg.input_prefix,'R*.mat'));
nSessions = length(fd);

popdir;

%% initialize variables
% for all non-baseline models:
nModels = length(cfg.models);

ALL_avgErr = nan(nModels, cfg.nMaxCells); % average error for each cell

ALL_ttrErr = nan(nModels, cfg.nTimeBins, cfg.nMaxCells); % nTimeBins x nCells matrix (time to reward-binned model improvement)
ALL_spaceErr = nan(nModels, cfg.nSpaceBins, cfg.nMaxCells); % nLocBins x nCells matrix (linearized position)
ALL_tstat = nan(nModels, 16, cfg.nMaxCells); % nVariables x nCells matrix (t-stats)

%for iM = 1:nModels
%end % over models

% counters
cellCount = 1; % MATLAB indexing starts at 1

%% loop across sessions to build variables

for iS = 1:nSessions

    load(fd{iS});
    baseline_err = sd.m.baseline.err;
    
    % indexing for populating ALL matrices
    nCells = size(baseline_err, 1);
    start_idx = cellCount; end_idx = cellCount + nCells - 1;
    
    for iM = 1:nModels
    
        this_diff = baseline_err - sd.m.(cfg.models{iM}).err; % error difference
        
        ALL_avgErr(iM,start_idx:end_idx) = nanmean(this_diff, 2);
        
        % ttr err
        this_f = ones(cfg_master.smooth, 1) ./ cfg_master.smooth;
        mdiffmean = conv(nanmean(this_err), this_f, 'same');
        [~, lambda_mttr, ~] = MakeTC_1D(cfg_ttr, sd.TVECc, t_to_reward, sd.TVECc, mdiffmean);
        
        % space err 
        
        % t-stats
           
        
    end % of models
    
    
    cellCount = cellCount + nCells;
end

%% plot