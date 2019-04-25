% COLLECT_rhythmGLMfit.m
%
% Collector script for output of ALL_rhythmGLMfit

cfg.inputDir = 'C:\temp\GLMfit\HC'; % where the files to load are
cfg.input_prefix = 'R0_'; 
%cfg.models = {'dphi','tphi','lgphi','hgphi','allphi'}; % models to be compared to baseline
cfg.models = {'hctheta'};
cfg.nMaxCells = 1000;
cfg.nTimeBins = 100;
cfg.nSpaceBins = 100;

%%
pushdir(cfg.inputDir);

fd = FindFiles(cat(2,cfg.input_prefix,'*.mat'));
nSessions = length(fd);

popdir;

%% initialize variables
% for all non-baseline models:
nModels = length(cfg.models);

ALL_avgErr = nan(nModels, cfg.nMaxCells); % average error for each cell

ALL_ttrErr = nan(nModels, cfg.nMaxCells, cfg.nTimeBins); % time to reward-binned model improvement
ALL_spaceErr = nan(nModels, cfg.nMaxCells, cfg.nSpaceBins); % linearized position-binned model improvement
ALL_tstat = nan(nModels, cfg.nMaxCells, 16); % t-stats

% envelope variables (but this may not make sense for all models...)
load(fd{1}); % load one file to get access to fieldnames
envnames = fieldnames(sd.env_ttr);
for iEnv = 1:length(envnames)
    
   ALL_env.(envnames{iEnv}) = nan(cfg.nMaxCells, cfg.nTimeBins);
    
end

% also need to add speed! include in sd

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
skipCells = 144;
ALL_avgErr(:, skipCells) = []; ALL_tstat(:, skipCells, :) = []; ALL_ttrErr(:, skipCells, :) = [];
cellCount = cellCount - length(skipCells);
cell_keep = 1:cellCount;

%%
for iM = 1:nModels
    
    figure;
    
    % improvement across cells
    subplot(221);
    this_err = ALL_avgErr(iM, cell_keep);
    
    errM = nanmean(this_err); errSD = nanstd(this_err);
    errZ = errM ./ errSD;
    p = signrank(this_err);
    fprintf('Model %s dErr %.2e +/- %.2e (z = %.2f), p = %.3e (nCells = %d)\n', cfg.models{iM}, errM, errSD, errZ, p, cellCount);
    
    high_idx = find(this_err > 0);
    low_idx = find(this_err <= 0);
    fprintf('%d (%.2f %%) cells improved\n', length(high_idx), (length(high_idx) ./ length(cell_keep)) .* 100);

    nHistBins = 100; histLims = [-4e-6 1e-5]; 
    histBinEdges = linspace(histLims(1), histLims(2), nHistBins);
    histBinCenters = histBinEdges + median(diff(histBinEdges))/2; histBinCenters = histBinCenters(1:end-1);
    low_idx = histBinCenters <= 0; high_idx = histBinCenters > 0;
    
    hData = histc(this_err, histBinEdges); hData = hData(1:end-1);
   
    h1 = bar(histBinCenters(low_idx), hData(low_idx)); set(h1, 'FaceColor', 'r', 'EdgeColor', 'none');
    hold on;
    h2 = bar(histBinCenters(high_idx), hData(high_idx)); set(h2, 'FaceColor', 'g', 'EdgeColor', 'none');
    
    set(gca,'LineWidth', 1, 'FontSize', 18, 'TickDir', 'out'); box off;
    ylabel('count'); xlabel('prediction improvement');
    title(sprintf('Model %s', cfg.models{iM}));
    
    % t-stat across cells
    subplot(222);
    
    this_tstat = sq(ALL_tstat(iM, cell_keep, :));
    this_vars = ALL_tstat_varnames{iM}; nVars = length(ALL_tstat_varnames{iM});
    this_tstat = this_tstat(1:cellCount-1,1:nVars);
    
    imagesc(this_tstat); colorbar; caxis([0 20]);
    
    trunc = @(x) x(1:2); lbl = cellfun(trunc, this_vars, 'UniformOutput', false);
    set(gca, 'LineWidth', 1, 'FontSize', 12, 'XTick', 1:nVars, 'XTickLabel', lbl, 'TickDir', 'out'); box off;
    ylabel('cell #');
    
    % average t-stats
    subplot(223);
    errorbar(1:nVars, nanmean(this_tstat), nanstd(this_tstat), 'k');
    hold on;
    plot(nanmean(this_tstat), '.k', 'MarkerSize', 20);
    set(gca, 'LineWidth', 1, 'FontSize', 18, 'XTick', 1:nVars, 'XTickLabel', lbl, 'TickDir', 'out'); box off;
    ylabel('mean t-stat');
    
    % t-stat correlations
    subplot(224);
    
    tcorr = corrcoef(this_tstat);
    imagesc(tcorr); colorbar;
    set(gca,'LineWidth',1,'FontSize',12,'XTick',1:nVars,'XTickLabel',lbl,'YTick',1:nVars,'YTickLabel',lbl, 'TickDir', 'out'); box off;
    
    %%%
    figure;
    
    this_ttr = sq(ALL_ttrErr(iM, cell_keep, :));
    xvec = sd.cfg.ttr_bins(1:end-1) + diff(sd.cfg.ttr_bins)/2;
    
    env_idx = strmatch(cfg.models{iM}(1:2), envnames); % find envelope that shares first letter with model name...

    if length(env_idx) == 1
        this_env = ALL_env.(envnames{env_idx});
        this_env = this_env(cell_keep, :);
        envname = envnames{env_idx};
    else
        this_env = nan(size(xvec,2));
        envname = 'none';
    end
    
    subplot(221);
    [ax h1 h2] = plotyy(xvec, nanmean(this_ttr), xvec, nanmean(this_env));
    hold on;
    set(h1,'LineWidth',2);
    set(gca,'XTick',-5:5,'LineWidth',1,'FontSize',18,'YLim',[0 3e-6],'YTick',0:1e-6:3e-6, 'TickDir', 'out'); box off;
    set(ax(2),'FontSize',18,'LineWidth',1,'YTick',[]);

    ylabel('Prediction improvement'); xlabel('time from reward (s)');
    
    title(sprintf('Model %s, env %s', cfg.models{iM}, envname));
    
    subplot(223);
    imagesc(xvec, 1:cellCount - 1, this_ttr(1:cellCount - 1,:));
    set(gca,'LineWidth', 1, 'XTick', -5:5, 'FontSize', 18, 'TickDir', 'out'); 
    xlabel('time from reward (s)'); ylabel('cell#'); box off;
    
    % normalize within each cell first, then average
    this_ttr = normalizeM(this_ttr);
    this_env = normalizeM(this_env);
    
    subplot(222);
    [ax h1 h2] = plotyy(xvec,nanmean(this_ttr),xvec,nanmean(this_env));
    hold on;
    set(h1,'LineWidth',2);
    set(gca,'XTick', -5:5, 'LineWidth', 1, 'FontSize',18, 'TickDir', 'out'); box off;
    set(ax(2), 'FontSize', 18,'LineWidth', 1, 'TickDir', 'out');
    ylabel('Prediction improvement'); xlabel('time from reward (s)');
    
    title(sprintf('Model %s (%s), env %s', cfg.models{iM}, cellTypes{iC}, envname));
    
    subplot(224);
    imagesc(xvec, 1:cellCount - 1, this_ttr(1:cellCount - 1,:));
    set(gca,'LineWidth',1,'XTick',-5:5,'FontSize',18); xlabel('time from reward (s)'); ylabel('cell#'); box off;
    
end