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

ALL_ttrErr = nan(nModels, cfg.nMaxCells, cfg.nTimeBins); % time to reward-binned model improvement
ALL_spaceErr = nan(nModels, cfg.nMaxCells, cfg.nSpaceBins); % linearized position-binned model improvement
ALL_tstat = nan(nModels, cfg.nMaxCells, 16); % t-stats

% counters
cellCount = 1; % MATLAB indexing starts at 1

%% loop across sessions to build variables

for iS = 1:nSessions

    load(fd{iS});
    
    % indexing for populating ALL matrices
    nCells = length(sd.S.t);
    start_idx = cellCount; end_idx = cellCount + nCells - 1;
    
    for iM = 1:nModels
    
        % average error
        this_diff = sd.m.baseline.err - sd.m.(cfg.models{iM}).err;
        ALL_avgErr(iM, start_idx:end_idx) = this_diff;
        
        % ttr err
        ALL_ttrErr(iM, start_idx:end_idx, :) = sd.m.(cfg.models{iM}).ttr_err;
        
        % space err 
        ALL_spaceErr(iM, start_idx:end_idx, :) = sd.m.(cfg.models{iM}).linpos_err;
        
        % t-stats
        nT = length(sd.m.(cfg.models{iM}).varnames);
        this_t = sd.m.(cfg.models{iM}).tstat(:,2:end);
        
        ALL_tstat(iM, start_idx:end_idx, 1:size(this_t,2)) = this_t;
        ALL_tstat_varnames{iM} = sd.m.(cfg.models{iM}).varnames;
        
    end % of models
    
    
    cellCount = cellCount + nCells;
end

%% plot
for iM = 1:nModels
    
    figure;
    
    % average improvement across cells
    subplot(221);
    this_err = ALL_avgErr(iM,:);
    
    plot(this_err,'LineWidth',1); hold on;
    keep = find(this_err > 0); plot(keep,this_err(keep),'.g','MarkerSize',20);
    keep = find(this_err < 0); plot(keep,this_err(keep),'.r','MarkerSize',20);
    set(gca,'LineWidth',1,'FontSize',18); box off;
    ylabel('Prediction improvement'); xlabel('cell #');
    title(cat(2,'Model ', cfg.models{iM}));
    
    % t-stat across cells
    subplot(222);
    
    this_tstat = sq(ALL_tstat(iM,:,:));
    this_vars = ALL_tstat_varnames{iM}; nVars = length(ALL_tstat_varnames{iM});
    this_tstat = this_tstat(1:cellCount-1,1:nVars);
    
    imagesc(this_tstat); colorbar; caxis([0 20]);
    
    trunc = @(x) x(1:2); lbl = cellfun(trunc, this_vars, 'UniformOutput', false);
    set(gca,'LineWidth',1,'FontSize',12,'XTick',1:nVars,'XTickLabel',lbl); box off;
    ylabel('cell #');
    
    % t-stat correlations
    subplot(223);
    
    tcorr = corrcoef(this_tstat);
    imagesc(tcorr);
    set(gca,'LineWidth',1,'FontSize',12,'XTick',1:nVars,'XTickLabel',lbl,'YTick',1:nVars,'YTickLabel',lbl); box off;
    
    figure;
    
    this_ttr = sq(ALL_ttrErr(iM,:,:));
    
    subplot(221);
    plot(nanmean(this_ttr));
    
    title(cat(2,'Model ', cfg.models{iM}));
   
    subplot(223);
    imagesc(this_ttr(1:cellCount - 1,:));
    
    %
    this_space = sq(ALL_spaceErr(iM,:,:));
    
    subplot(222);
    plot(nanmean(this_space));
    
    subplot(224);
    imagesc(this_space(1:cellCount - 1,:));
    
end