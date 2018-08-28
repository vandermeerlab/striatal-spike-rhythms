%% setup
clear;
% paths: striatal-spike-rhythms

% load data
%cd('D:\data\adrlab\R117\R117-2007-06-20');
cd('D:\data\adrlab\R132\R132-2007-10-20');
LoadExpKeys;

% spikes -- should remove cells with < 100 spikes on the track here
% (currently in main cell loop)
S = LoadSpikes([]);
S = KeepCells(S);
%S.t = S.t(setxor(1:27,[12 16]));

% LFP
cfg = []; cfg.fc = ExpKeys.goodGamma_vStr;
csc = LoadCSC(cfg);

% position
pos = LoadPos([]);

% reward deliveries
evt = LoadEvents([]);
reward_t = evt.t{1}; % should generalize this with known labels etc.. and remove double labels

%% params
cfg_master = []; % overall params
cfg_master.dt = 0.001;
cfg_master.f.delta = [2 5];
cfg_master.f.theta = [6.5 9.5];
cfg_master.f.beta = [14 25];
cfg_master.f.lowGamma = [40 65];
cfg_master.f.highGamma = [70 100];
cfg_master.nPleats = 5;
cfg_master.kFold = 2;

cfg_tcx = []; % tuning curve for x position
cfg_tcx.bins = 50:10:600;

cfg_phi = []; % LFP features
cfg_phi.dt = median(diff(csc.tvec));
cfg_phi.ord = 100;
cfg_phi.bins = -pi:pi/18:pi;
cfg_phi.interp = 'nearest';

cfg_env = [];
cfg_env.bins = 100;
cfg_env.interp = 'nearest';

%% initialize variables
clear mdiff; % session-wide variable for model differences

%% establish common timebase and other variables common across cells in this session
TVECe = ExpKeys.TimeOnTrack:cfg_master.dt:ExpKeys.TimeOffTrack; % edges
TVECc = TVECe(1:end-1)+cfg_master.dt/2; % centers

% time-to-reward variable is used to restrict samples
t_to_reward = ComputeTTR(reward_t,TVECc);
cfg_ttr = []; cfg_ttr.interp = 'nearest'; cfg_ttr.bins = [-5:0.1:5];
spk_binned = logical(zeros(size(TVECc))); % dummy spike train just to get binned TTR
[ttr, lambda_ttr, ttr_binned] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, spk_binned);

MASTER_keep = ~isnan(ttr); % NEED THIS TO RESTRICT BINNED SPIKE TRAIN LATER!

t_to_reward = t_to_reward(MASTER_keep);
TVECc = TVECc(MASTER_keep);

for iPleat = 1:cfg_master.nPleats
    C{iPleat} = cvpartition(length(TVECc),'KFold',cfg_master.kFold);
end

% LFP features
disp('Computing session-wide LFP features...');
fb_names = fieldnames(cfg_master.f);
cfg_phi.debug = 0;
for iF = 1:length(fb_names)
   
    cfg_phi.fpass = cfg_master.f.(fb_names{iF});
    fprintf('%s: [%d %d]\n', fb_names{iF}, cfg_phi.fpass(1), cfg_phi.fpass(2));
    
    [phase, env] = ComputePhase(cfg_phi, csc);
    
    FF.(fb_names{iF}).phase = phase; 
    FF.(fb_names{iF}).env = env;

end

%
spd = getLinSpd([],pos);

%%% TODO: time / nTrials predictor

%% loop over all cells 
for iC = length(S.t):-1:1

% dependent variable: binned spike train
spk_binned = histc(S.t{iC}, TVECe); spk_binned = spk_binned(1:end - 1);

if sum(spk_binned) < 100
    fprintf(' *** Cell %d skipped because of low number of spikes.\n',iC);
    continue;
end

spk_binned = logical(spk_binned > 0); % binarize
%k = gausskernel(50,2); % smoothing works too
%spk_binned = conv2(spk_binned,k,'same');

spk_binned = spk_binned(MASTER_keep);

%%% PREDICTOR: x- and y-location %%%
cfg_tc = cfg_tcx;
cfg_tc.interp = 'linear';
[predicted_x, lambda_x, x_binned] = MakeTC_1D(cfg_tc, pos.tvec, pos.data(1,:), TVECc, spk_binned);
predicted_x = predicted_x-nanmean(predicted_x);

cfg_tc.bins = 0:10:480;
[predicted_y, lambda_y, y_binned] = MakeTC_1D(cfg_tc, pos.tvec, pos.data(2,:), TVECc, spk_binned);
predicted_y = predicted_y-nanmean(predicted_y);

%%% TODO: linpos predictor %%%

% how to deal with left vs. right laps? fit all data separately with L and R coords, and then
% combine (UnionTSD -- too far samples will be NaN, and duplicates can be removed)

%%% PREDICTOR: speed %%%
cfg_spd = cfg_tc;
cfg_spd.bins = 0:10:200;
[pred_spd, ~, spd_binned] = MakeTC_1D(cfg_spd, pos.tvec, spd.data, TVECc, spk_binned);
pred_spd = pred_spd-nanmean(pred_spd);

%%% PREDICTOR: LFP phases and envelopes (e.g. gamma) %%%
for iF = 1:length(fb_names) % loop over frequency bands
    lfp_pred.(fb_names{iF}) = MakeTC_1D(cfg_phi, csc.tvec, FF.(fb_names{iF}).phase, TVECc, spk_binned);
    lfp_pred_env.(fb_names{iF}) = MakeTC_1D(cfg_env, csc.tvec, FF.(fb_names{iF}).env, TVECc, spk_binned);
end

%%% POSSIBLE PREDICTOR: LFP envelope

%%% PREDICTOR: basic conditional intensity function %%%
cif = conv(spk_binned,[0 0 0 1 1],'same'); % 2 ms refractory period
cif = (cif > 0);

%%% PREDICTOR: better CIF based on acorr
% why does this work so well? 
% what happens for a symmetric acorr?
% what happens when x-val is done on "blocked" rather than random folds?
% how do results depend on this thing being present or not?
cfg_acf = []; cfg_acf.maxlag = 100; cfg_acf.binsize = cfg_master.dt;
cfg_acf.sided = 'onezero';
[acf, acf_tvec] = ComputeACF(cfg_acf, spk_binned);

cif_full = conv(spk_binned, acf ./ sum(acf), 'same');
cif_full = cif_full-nanmean(cif_full);
cif_full(cif) = 0; % absolute refractory period

%%% PREDICTOR: time to reward
[ttr, lambda_ttr, ttr_binned] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, spk_binned);
ttr = ttr-nanmean(ttr);

%% x-val loop
mdiff.full(iC,:) = nan(size(spk_binned)); % model fit difference for full model
for iF = 1:length(fb_names)
   mdiff.(fb_names{iF})(iC,:) = nan(size(spk_binned));
   mdiff_env.(fb_names{iF})(iC,:) = nan(size(spk_binned));
   mdiff_both.(fb_names{iF})(iC,:) = nan(size(spk_binned));
end
    
base_predictors = cat(1,predicted_x,predicted_y,ttr,cif_full',pred_spd);

LFP_predictors = []; LFP_predictors_env = [];
for iF = 1:length(fb_names)
    LFP_predictors = cat(1,LFP_predictors,lfp_pred.(fb_names{iF}));
    LFP_predictors_env = cat(1,LFP_predictors_env,lfp_pred_env.(fb_names{iF}));
end
   
% temporary variables to store error across pleats and folds
this_mdiff.full = zeros(size(spk_binned)); % running total of model diffs (across pleats and folds)
for iF = 1:length(fb_names)
   this_mdiff.(fb_names{iF}) = zeros(size(spk_binned));
   this_mdiff_env.(fb_names{iF}) = zeros(size(spk_binned));
   this_mdiff_both.(fb_names{iF}) = zeros(size(spk_binned));
end
this_mdiff.count = zeros(size(spk_binned)); % count of test idxs (across pleats and folds, to compute mean later)    

for iPleat = 1:cfg_master.nPleats
    
for iFold = 1:C{iPleat}.NumTestSets
    
    fprintf('Pleat %d fold %d...\n', iPleat, iFold);
    
    % first, train models
    tr_idx = C{iPleat}.training(iFold);
    
    % baseline model without LFP features
    m2 = fitglm(base_predictors(:,tr_idx)', spk_binned(tr_idx)', 'Distribution', 'binomial')
    
    % big model with all the LFP features so we can correlate them
    m3 = fitglm(cat(1, base_predictors(:,tr_idx), LFP_predictors(:,tr_idx), LFP_predictors_env(:,tr_idx))', spk_binned(tr_idx)', 'Distribution', 'binomial')
    OUT_tstat(iC,:) = m3.Coefficients.tStat;
    
    % baseline + specific LFP feature models
    for iF = 1:length(fb_names)
        mLFP{iF} = fitglm(cat(1, base_predictors(:,tr_idx), LFP_predictors(iF,tr_idx))', spk_binned(tr_idx)', 'Distribution', 'binomial');
        mLFP_env{iF} = fitglm(cat(1, base_predictors(:,tr_idx), LFP_predictors_env(iF,tr_idx))', spk_binned(tr_idx)', 'Distribution', 'binomial');
        mLFP_both{iF} = fitglm(cat(1, base_predictors(:,tr_idx), LFP_predictors(iF,tr_idx), LFP_predictors_env(iF,tr_idx))', spk_binned(tr_idx)', 'Distribution', 'binomial');
    end
    
    %%% test models %%%
    te_idx = C{iPleat}.test(iFold);
    
    m2p = m2.predict(base_predictors(:,te_idx)');
    sse2 = (m2p-spk_binned(te_idx)).^2;
    
    m3p = m3.predict(cat(1, base_predictors(:,te_idx), LFP_predictors(:,te_idx), LFP_predictors_env(:,te_idx))');
    sse3 = (m3p-spk_binned(te_idx)).^2;
    
    this_mdiff.full(te_idx) = this_mdiff.full(te_idx) + (sse2 - sse3);
    this_mdiff.count(te_idx) = this_mdiff.count(te_idx) + 1;
    %mdiff.full(iC,te_idx) = sse2-sse3; % overall model improv
    
    for iF = 1:length(fb_names)
    
        % phase
        this_p = mLFP{iF}.predict(cat(1, base_predictors(:,te_idx), LFP_predictors(iF,te_idx))');
        this_sse = (this_p-spk_binned(te_idx)).^2;
        
        this_mdiff.(fb_names{iF})(te_idx) = this_mdiff.(fb_names{iF})(te_idx) + (sse2 - this_sse);
        %mdiff.(fb_names{iF})(iC,te_idx) = sse2-this_sse;
        
        % env
        this_p = mLFP_env{iF}.predict(cat(1, base_predictors(:,te_idx), LFP_predictors_env(iF, te_idx))');
        this_sse = (this_p-spk_binned(te_idx)).^2;
        
        this_mdiff_env.(fb_names{iF})(te_idx) = this_mdiff_env.(fb_names{iF})(te_idx) + (sse2 - this_sse);
        
        % both
        this_p = mLFP_both{iF}.predict(cat(1, base_predictors(:,te_idx), LFP_predictors(iF, te_idx), LFP_predictors_env(iF, te_idx))');
        this_sse = (this_p-spk_binned(te_idx)).^2;
        
        this_mdiff_both.(fb_names{iF})(te_idx) = this_mdiff_both.(fb_names{iF})(te_idx) + (sse2 - this_sse);
    end

end % over folds

plot(this_mdiff.full ./ this_mdiff.count);
title(sprintf('cell %d pleat %d', iC, iPleat));
drawnow;

end % over pleats

% compute average over pleats and folds
mdiff.full(iC,:) = this_mdiff.full ./ this_mdiff.count;
for iF = 1:length(fb_names)
    mdiff.(fb_names{iF})(iC,:) = this_mdiff.(fb_names{iF}) ./ this_mdiff.count;
    mdiff_env.(fb_names{iF})(iC,:) = this_mdiff_env.(fb_names{iF}) ./ this_mdiff.count;
    mdiff_both.(fb_names{iF})(iC,:) = this_mdiff_both.(fb_names{iF}) ./ this_mdiff.count;
end

end % over cells

%% make correlation matrix of GLM predictors
all_pred = cat(1,base_predictors,LFP_predictors,LFP_predictors_env);
cm = corrcoef(all_pred');

figure;
imagesc(cm);

set(gca,'XTick',1:16,'XTickLabel',{'x','y','tr','ac','s','d','t','b','lg','hg','de','te','be','lge','hge'});
set(gca,'YTick',1:16,'YTickLabel',{'x','y','tr','ac','s','d','t','b','lg','hg','de','te','be','lge','hge'});

%% plot performance increase of full model across all bands
figure;
celldiffmean = nanmean(mdiff.full,2);

cellskip = [];
%cellskip = [6 7 19]; %celldiffmean < 0; %[]; % TODO: figure out why some model diffs end up being NaN -- model fit fail?
%mdiff.full(cellskip,:) = NaN;

celldiffmean = nanmean(mdiff.full,2);
subplot(221);
plot(celldiffmean,'LineWidth',1); hold on; 
keep = find(celldiffmean > 0);
plot(keep,celldiffmean(keep),'.g','MarkerSize',20);
keep = find(celldiffmean < 0);
plot(keep,celldiffmean(keep),'.r','MarkerSize',20);
set(gca,'LineWidth',1,'FontSize',18); box off;
ylabel('Prediction improvement with LFP features added');
xlabel('cell #');

%%% TODO: better ways of filtering this. Shouldn't go across trial breaks.
%mdiffmean = (medfilt1(nanmean(mdiff.full),1001));
nF = 501; this_f = ones(nF,1)./nF;
mdiffmean = conv(nanmean(mdiff.full),this_f,'same');

[~, lambda_mttr, ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, mdiffmean);
[~, lambda_spd, ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, spd_binned);

subplot(222)
[ax h1 h2] = plotyy(cfg_ttr.bins,lambda_mttr,cfg_ttr.bins,lambda_spd);
hold on;
set(h1,'LineWidth',2);
set(gca,'XTick',-5:5,'LineWidth',1,'FontSize',18);
ylabel('Prediction improvement with LFP features added');
xlabel('time from reward (s)');
box off;

subplot(223);
OUT_tstat(cellskip,:) = NaN;
imagesc(OUT_tstat); caxis([0 10]);
set(gca,'XTick',1:11,'XTickLabel',{'i','x','y','tr','ac','s','d','t','b','lg','hg'});

subplot(224);
imagesc(corrcoef(OUT_tstat(setxor(cellskip,1:length(celldiffmean)),:)));
set(gca,'XTick',1:16,'XTickLabel',{'i','x','y','tr','ac','s','d','t','b','lg','hg','de','te','be','lge','hge'});
set(gca,'YTick',1:16,'YTickLabel',{'i','x','y','tr','ac','s','d','t','b','lg','hg','de','te','be','lge','hge'});

%% normalized by firing rate
if 0
    for iC = length(S.t):-1:1
        this_spk_binned = histc(S.t{iC}, TVECe); spk_binned = spk_binned(1:end - 1);
        [~, lambda_ttr(iC,:), ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, this_spk_binned);
        
        nF = 501; this_f = ones(nF,1)./nF;
        mdiffmean = conv(mdiff.full(iC,:),this_f,'same');
        [~, lambda_ttr_mdiff(iC,:), ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, mdiffmean);
    end
    mdiff_fullN = lambda_ttr_mdiff ./ lambda_ttr;
    mdiff_fullN(isinf(mdiff_fullN)) = NaN;
    
    figure; subplot(222);
    [ax h1 h2] = plotyy(cfg_ttr.bins,nanmean(mdiff_fullN),cfg_ttr.bins,lambda_spd);
    hold on;
    set(h1,'LineWidth',2);
    set(gca,'XTick',-5:5,'LineWidth',1,'FontSize',18);
    ylabel('Prediction improvement with LFP features added');
    xlabel('time from reward (s)');
    box off;
    
    subplot(223);
    imagesc(lambda_ttr);
    
    subplot(224);
    plot(cfg_ttr.bins,nanmean(lambda_ttr));
    set(gca,'XTick',-5:5,'LineWidth',1,'FontSize',18); box off;
    ylabel('frate');
    xlabel('time from reward (s)');
end
%% different freq bands - contribution of phase
for iF = 1:length(fb_names)

    mdiff.(fb_names{iF})(cellskip,:) = NaN;
    
    nF = 501; this_f = ones(nF,1)./nF;
    mdiffmean = conv(nanmean(mdiff.(fb_names{iF})),this_f,'same');

    %mdiffmean = (medfilt1(nanmean(mdiff.(fb_names{iF})),1001));
    [~, this_mttr, ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, mdiffmean);
    [~, this_env, ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, csc.tvec, FF.(fb_names{iF}).env);
        
    figure;
    
    subplot(221);
    plot(nanmean(mdiff.(fb_names{iF}),2));
    hold on;
    plot(nanmean(mdiff.(fb_names{iF}),2),'.b','MarkerSize',20);
    
    subplot(222); % average improvement
    [ax h1 h2] = plotyy(cfg_ttr.bins,this_mttr,cfg_ttr.bins,this_env);
    set(h1,'LineWidth',2);
    hold on;
    set(gca,'XTick',-5:5,'LineWidth',1,'FontSize',18);
    ylabel('Prediction improvement with LFP features added');
    xlabel('time from reward (s)');
    title([fb_names{iF} ' phase']);
    box off;

    nCells = size(mdiff.full, 1);
    subplot(223); % imagesc across cells
    for iC = nCells:-1:1
        fprintf('band %d, filtering cell %d...\n', iF, iC);
        
        this_flt = ones(nF,1)./nF;
        this_f = conv(mdiff.(fb_names{iF})(iC,:),this_flt,'same');

        %this_f = medfilt1(mdiff.(fb_names{iF})(iC,:), 1001);
        
        [~, mdiffm_f(iC,:), ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, this_f);
    end
    imagesc(cfg_ttr.bins, 1:nCells, mdiffm_f);
    
    subplot(224);
    % normalized version of above matrix
    imagesc(cfg_ttr.bins, 1:nCells, zmat(mdiffm_f, 'z'));

end

%% different freq bands - contribution of env
for iF = 1:length(fb_names)

    mdiff_env.(fb_names{iF})(cellskip,:) = NaN;
    
    nF = 501; this_f = ones(nF,1)./nF;
    mdiffmean = conv(nanmean(mdiff_env.(fb_names{iF})),this_f,'same');

    %mdiffmean = (medfilt1(nanmean(mdiff.(fb_names{iF})),1001));
    [~, this_mttr, ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, mdiffmean);
    [~, this_env, ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, csc.tvec, FF.(fb_names{iF}).env);
        
    figure;
    
    subplot(221);
    plot(nanmean(mdiff_env.(fb_names{iF}),2));
    hold on;
    plot(nanmean(mdiff_env.(fb_names{iF}),2),'.b','MarkerSize',20);
    
    subplot(222); % average improvement
    [ax h1 h2] = plotyy(cfg_ttr.bins,this_mttr,cfg_ttr.bins,this_env);
    set(h1,'LineWidth',2);
    hold on;
    set(gca,'XTick',-5:5,'LineWidth',1,'FontSize',18);
    ylabel('Prediction improvement with LFP features added');
    xlabel('time from reward (s)');
    title([fb_names{iF} ' envelope']);
    box off;

    nCells = size(mdiff.full, 1);
    subplot(223); % imagesc across cells
    for iC = nCells:-1:1
        fprintf('band %d, filtering cell %d...\n', iF, iC);
        
        this_flt = ones(nF,1)./nF;
        this_f = conv(mdiff_env.(fb_names{iF})(iC,:),this_flt,'same');

        %this_f = medfilt1(mdiff.(fb_names{iF})(iC,:), 1001);
        
        [~, mdiffm_f(iC,:), ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, this_f);
    end
    imagesc(cfg_ttr.bins, 1:nCells, mdiffm_f);
    
    subplot(224);
    % normalized version of above matrix
    imagesc(cfg_ttr.bins, 1:nCells, zmat(mdiffm_f, 'z'));

end

%% different freq bands - contribution of both
for iF = 1:length(fb_names)

    mdiff_both.(fb_names{iF})(cellskip,:) = NaN;
    
    nF = 501; this_f = ones(nF,1)./nF;
    mdiffmean = conv(nanmean(mdiff_both.(fb_names{iF})),this_f,'same');

    %mdiffmean = (medfilt1(nanmean(mdiff.(fb_names{iF})),1001));
    [~, this_mttr, ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, mdiffmean);
    [~, this_env, ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, csc.tvec, FF.(fb_names{iF}).env);
        
    figure;
    
    subplot(221);
    plot(nanmean(mdiff_both.(fb_names{iF}),2));
    hold on;
    plot(nanmean(mdiff_both.(fb_names{iF}),2),'.b','MarkerSize',20);
    
    subplot(222); % average improvement
    [ax h1 h2] = plotyy(cfg_ttr.bins,this_mttr,cfg_ttr.bins,this_env);
    set(h1,'LineWidth',2);
    hold on;
    set(gca,'XTick',-5:5,'LineWidth',1,'FontSize',18);
    ylabel('Prediction improvement with LFP features added');
    xlabel('time from reward (s)');
    title([fb_names{iF} ' both']);
    box off;

    nCells = size(mdiff.full, 1);
    subplot(223); % imagesc across cells
    for iC = nCells:-1:1
        fprintf('band %d, filtering cell %d...\n', iF, iC);
        
        this_flt = ones(nF,1)./nF;
        this_f = conv(mdiff_both.(fb_names{iF})(iC,:),this_flt,'same');

        %this_f = medfilt1(mdiff.(fb_names{iF})(iC,:), 1001);
        
        [~, mdiffm_f(iC,:), ~] = MakeTC_1D(cfg_ttr, TVECc, t_to_reward, TVECc, this_f);
    end
    imagesc(cfg_ttr.bins, 1:nCells, mdiffm_f);
    
    subplot(224);
    % normalized version of above matrix
    imagesc(cfg_ttr.bins, 1:nCells, zmat(mdiffm_f, 'z'));

end

%% what should go into an output structure, to be saved for this session?
% mdiff, and all the predictors -- need to make tuning curves for them later
% spike train

%% correlate GLM results with things like spike spectrum, STA spectrum (across cells)

%% correlate tuning curves with GLM results
% why are the TTR GLM fits so "streaky"?
% is there some correlation between TTR, space, movement onset TCs

%% is there some way to normalize for number of spikes?

%% some functions
function [predicted_x, lambda_x, x_binned] = MakeTC_1D(cfg, tvec, data, TVECc, spk_binned)
% construct tuning curve (mean value of spk_binned for binned values of
% data
%
% INPUTS:
%
% tvec, data: describe tuning variable (e.g. position, phase)
% TVECc, spk_binned: describe spike train
%
% cfg.interp: option to use for interpolation, 'linear', 'nearest' (for
%   things like phase)
% cfg.bins: bin edges (vector) for constructing tuning curves, or number of bins (scalar) for autoscale
%
% OUTPUTS:
%
% x_binned: binned version of data input
% lambda_x: tuning curve (size cfg.bins), gives spike probability as
%   function of binned tuning variable
% predicted_x: predicted value of tuning variable for each TVECc input
%   sample (simply looked up in lambda_x)

if length(cfg.bins) == 1 % number of bins specified
    cfg.bins = linspace(min(data),max(data),cfg.bins+1);
end

x_interp = interp1(tvec, data, TVECc, cfg.interp); % get variable of interest on common timebase 
[x_occ_hist, x_binned] = histc(x_interp, cfg.bins); % bin variable of interest
x_binned(x_interp < cfg.bins(1) | x_interp >= cfg.bins(end)) = NaN;

bin_edges = 0.5:1:length(cfg.bins) + 0.5; % edges for binned version of tuning variable

if islogical(spk_binned) % binned spike train (binary), can use fast method
    x_spk = x_binned(spk_binned); % variable of interest at spike times
    x_spk_hist = histc(x_spk, bin_edges); x_spk_hist = x_spk_hist(1:end - 1);
    lambda_x = x_spk_hist ./ x_occ_hist;
else % continuous input, need to average across samples in each bin (slow)
    for iB = length(cfg.bins):-1:1
        this_samples = (x_binned == iB);
        lambda_x(iB) = nanmean(spk_binned(this_samples));
    end
end

% based on tuning curve, get lambdas for each time bin
predicted_x = nan(size(x_binned));
keep = x_binned ~= 0 & ~isnan(x_binned);
predicted_x(keep) = lambda_x(x_binned(keep));

end

function [phase, env] = ComputePhase(cfg, csc)
% compute csc.data phase and envelope in frequency band cfg.fpass

Fs = (1 ./ cfg.dt);

d = fdesign.bandpass('N,Fp1,Fp2,Ap', 10, cfg.fpass(1), cfg.fpass(2), .5, Fs);
Hd = design(d, 'cheby1');

if cfg.debug
   figure;
   fvtool(Hd);
end

data = filtfilt(Hd.sosMatrix, Hd.ScaleValues, csc.data);
h = hilbert(data);
phase = angle(h);
env = abs(h);

end

function t_to_reward = ComputeTTR(reward_t,TVECc)
% find time to closest reward at each time point in TVECc

t_to_reward = TVECc-reward_t(1); % running total
for iR = 2:length(reward_t)
    
    this_dt = TVECc-reward_t(iR);
    swap_idx = find(abs(this_dt) < abs(t_to_reward)); % samples that are closer than running total
    t_to_reward(swap_idx) = this_dt(swap_idx);
end

end

%%
function [acf,tvec] = ComputeACF(cfg,spk_binned)

if isfield(cfg,'maxlag')
    [acf,tvec] = xcorr(spk_binned,spk_binned,cfg.maxlag);
else
    [acf,tvec] = xcorr(spk_binned,spk_binned);
end

tvec = tvec.*cfg.binsize;
acf(ceil(length(acf)/2)) = 0;

if isfield(cfg,'sided')
    cut_idx = ceil(length(acf)/2);
    switch cfg.sided
        case 'one'
            acf = acf(cut_idx:end);
            tvec = tvec(cut_idx:end);
        case 'onezero'
            acf(1:cut_idx-1) = 0;
    end
end

end

%%
function out = zmat(in,cfg)
% z-score each row of input matrix
out = nan(size(in));
for iRow = size(in, 1):-1:1
    this_row = in(iRow,:);
    switch(cfg)
        case 'z'
            out(iRow,:) = (this_row - nanmean(this_row)) ./ nanstd(this_row);
        case 'max'
            out(iRow,:) = this_row ./ max(this_row);
    end
end
end
% could speed up by repmatting

function S = KeepCells(S)
thr = 500;

nCells = length(S.t);
for iC = nCells:-1:1
   l(iC) = length(S.t{iC}); 
end
keep = l >= thr;

S.label = S.label(keep);
S.t = S.t(keep);
end