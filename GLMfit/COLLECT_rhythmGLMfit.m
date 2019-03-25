% COLLECT_rhythmGLMfit.m
%
% Collector script for output of ALL_rhythmGLMfit

clear;

cfg.inputDir = 'C:\temp\GLMfit'; % where the files to load are
cfg.input_prefix = 'R0_'; 
cfg.models = {'dphi','tphi','lgphi','hgphi','allphi'}; % models to be compared to baseline
cfg.nMaxCells = 1000;
cfg.nTimeBins = 100;
cfg.nSpaceBins = 100;

%%
pushdir(cfg.inputDir);

fd = FindFiles(cat(2,cfg.input_prefix,'*.mat'), 'CheckSubdirs', 0);
nSessions = length(fd);

popdir;

%% initialize variables
% for all non-baseline models:
nModels = length(cfg.models);

ALL_avgErr = nan(nModels, cfg.nMaxCells); % average error for each cell

ALL_ttrErr = nan(nModels, cfg.nMaxCells, cfg.nTimeBins); % time to reward-binned model improvement
ALL_spaceErr = nan(nModels, cfg.nMaxCells, cfg.nSpaceBins); % linearized position-binned model improvement
ALL_tstat = nan(nModels, cfg.nMaxCells, 16); % t-stats

ALL_celltype = nan(nModels, cfg.nMaxCells); % redundant across models but will be easier to deal with later

% envelope variables (but this may not make sense for all models...)
load(fd{1}); % load one file to get access to fieldnames
envnames = fieldnames(sd.env_ttr);
for iEnv = 1:length(envnames)
    
   ALL_env.(envnames{iEnv}) = nan(cfg.nMaxCells, cfg.nTimeBins);
    
end

% also need to add speed! include ttr and linpos versions in sd

% counters
cellCount = 1; % MATLAB indexing starts at 1

%% loop across sessions to build variables

for iS = 1:nSessions

    load(fd{iS});
    
    % indexing for populating ALL matrices
    keep = ~all(isnan(sd.m.baseline.tstat),2); % find cells that were actually fit
    nCells = sum(keep);
    start_idx = cellCount; end_idx = cellCount + nCells - 1;
    
    for iM = 1:nModels
    
        % average error
        this_diff = sd.m.baseline.err(keep, :) - sd.m.(cfg.models{iM}).err(keep, :);
        ALL_avgErr(iM, start_idx:end_idx) = this_diff;
        
        % ttr err
        ALL_ttrErr(iM, start_idx:end_idx, :) = sd.m.(cfg.models{iM}).ttr_err(keep, :);
        
        % space err 
        ALL_spaceErr(iM, start_idx:end_idx, :) = sd.m.(cfg.models{iM}).linpos_err(keep, :);
        
        % t-stats
        nT = length(sd.m.(cfg.models{iM}).varnames);
        this_t = sd.m.(cfg.models{iM}).tstat(keep, 2:end);
        
        if any(isnan(this_t))
           this_t = this_t(1:nT);
        end
        
        ALL_tstat(iM, start_idx:end_idx, 1:size(this_t, 2)) = this_t;
        ALL_tstat_varnames{iM} = sd.m.(cfg.models{iM}).varnames;
        
        % cell type
        ALL_cellType(iM, start_idx:end_idx) = sd.cellType;
        
    end % of models

    % collect session-wide variables (envelopes, speed...)
    % so need to repmat across cells
    
    for iE = 1:length(envnames)
       
        this_envdata = sd.env_ttr.(envnames{iE});
        this_envdata = repmat(this_envdata, [nCells 1]);
        ALL_env.(envnames{iE})(start_idx:end_idx, :) = this_envdata;
        
    end
    
    cellCount = cellCount + nCells;

end

%% plot
cellTypes = {'MSNs', 'FSIs', 'all'};

for iC = 1:length(cellTypes)

    switch cellTypes{iC}
       
        case 'MSNs'
            cell_keep = find(ALL_cellType(1,:) == 1);
        case 'FSIs'
            cell_keep = find(ALL_cellType(1,:) == 2);
        otherwise
            cell_keep = find(~isnan(ALL_cellType(1,:)));
    end
    
    cellCount = length(cell_keep);
    
for iM = 1:nModels
    
    figure;
    
    % average improvement across cells
    subplot(221);
    this_err = ALL_avgErr(iM, cell_keep);
    
    plot(this_err, 'LineWidth', 1); hold on;
    keep = find(this_err > 0); plot(keep, this_err(keep), '.g', 'MarkerSize', 20);
    keep = find(this_err < 0); plot(keep, this_err(keep), '.r', 'MarkerSize', 20);
    set(gca, 'LineWidth', 1, 'FontSize', 18); box off;
    ylabel('Prediction improvement'); xlabel('cell #');
    title(sprintf('Model %s (%s)', cfg.models{iM}, cellTypes{iC}));
    
    % t-stat across cells
    subplot(222);
    
    this_tstat = sq(ALL_tstat(iM, cell_keep, :));
    this_vars = ALL_tstat_varnames{iM}; nVars = length(ALL_tstat_varnames{iM});
    this_tstat = this_tstat(1:cellCount-1,1:nVars);
    
    imagesc(this_tstat); colorbar; caxis([0 20]);
    
    trunc = @(x) x(1:2); lbl = cellfun(trunc, this_vars, 'UniformOutput', false);
    set(gca, 'LineWidth', 1, 'FontSize', 12, 'XTick', 1:nVars, 'XTickLabel', lbl); box off;
    ylabel('cell #');
    
    % t-stat correlations
    subplot(223);
    
    tcorr = corrcoef(this_tstat);
    imagesc(tcorr);
    set(gca,'LineWidth',1,'FontSize',12,'XTick',1:nVars,'XTickLabel',lbl,'YTick',1:nVars,'YTickLabel',lbl); box off;
    
    figure;
    
    this_ttr = sq(ALL_ttrErr(iM, cell_keep, :));
    xvec = sd.cfg.ttr_bins(1:end-1) + diff(sd.cfg.ttr_bins)/2;
    
    env_idx = strmatch(cfg.models{iM}(1), envnames); % find envelope that shares first letter with model name...
    if length(env_idx) == 1
        this_env = ALL_env.(envnames{env_idx});
        envname = envnames{env_idx};
    else
        this_env = nan(size(xvec,2));
        envname = 'none';
    end
    
    subplot(221);
    [ax h1 h2] = plotyy(xvec, nanmean(this_ttr), xvec, nanmean(this_env));
    hold on;
    set(h1,'LineWidth',2);
    %set(gca,'XTick',-5:5,'LineWidth',1,'FontSize',18,'YLim',[0 1.5e-6],'YTick',0:0.5e-6:1.5e-6); box off;
    set(gca,'XTick',-5:5,'LineWidth',1,'FontSize',18); box off;
    ylabel('Prediction improvement'); xlabel('time from reward (s)');
    
    title(sprintf('Model %s (%s), env %s', cfg.models{iM}, cellTypes{iC}, envname));
    
    subplot(223);
    imagesc(this_ttr(1:cellCount - 1,:));
    
    % normalize within each cell first, then average
    this_ttr = normalizeM(this_ttr);
    this_env = normalizeM(this_env);
    
    subplot(222);
    [ax h1 h2] = plotyy(xvec,nanmean(this_ttr),xvec,nanmean(this_env));
    hold on;
    set(h1,'LineWidth',2);
    set(gca,'XTick',-5:5,'LineWidth',1,'FontSize',18); box off;
    ylabel('Prediction improvement'); xlabel('time from reward (s)');
    
    title(sprintf('Model %s (%s), env %s', cfg.models{iM}, cellTypes{iC}, envname));
    
    subplot(224);
    imagesc(this_ttr(1:cellCount - 1,:));
    
end

end % of cell types