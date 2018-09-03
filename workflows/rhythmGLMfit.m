function sd = rhythmGLMfit(cfg_in)
% GLM for spike prediction with LFP features
%
% This top-level function fits a number of GLMs to single session spike train data.
%
% The overall idea is to test whether addition of LFP features, such as
% theta power and/or gamma phase, improve cross-validated model prediction.
%
% The analysis proceeds as follows:
% - load data
% - define models (in sd.m.modelspec), MUST include a baseline model
% - prepare session-wide variables (linearized position, LFP features,
%   speed etc) on common timebase ('TVECc')
% - for each cell:
%   * prepare regressors for this cell
%   - for each cross-validation run ("pleat"; a set of folds) and fold:
%     + fit models on training data
%     + test models
%     + store error relative to baseline model
% - for each model:
%   * plot error across cells
%   * plot error across time to reward (tuning curve)
%
% output error is stored in the sd variable, which can be saved for a later
% collector to aggregate across sessions.
%
% maybe storing errors is overkill though: too much space, could store
% tuning curves for various things instead
%
% run this with data session to analyze as the working folder.
%
% required path: striatal-spike-rhythms github repo

% load data
%cd('C:\data\adrlab\R117-2007-06-20');
%cd('D:\data\adrlab\R132\R132-2007-10-20');
LoadExpKeys;

% spikes
sd.S = LoadSpikes([]);
sd.S = KeepCells(sd.S); % note this has a hardcoded parameter (500 spikes minimum)
nCells = length(sd.S.t);

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
cfg_master.nPleats = 2;
cfg_master.kFold = 2;
cfg_master.writeOutput = 0;
cfg_master.output_prefix = 'R0_';
cfg_master.output_dir = 'C:\temp';

cfg_phi = []; % LFP features
cfg_phi.dt = median(diff(csc.tvec));
cfg_phi.ord = 100;
cfg_phi.bins = -pi:pi/18:pi;
cfg_phi.interp = 'nearest';

%% initialize variables
% overall plan:
%
% for each cell, p table contains regressors. these are NOT stored across cells (into session data) for memory reasons
% but maybe session-wide variables should be stored (time_to_reward,
% linpos, LFP features) for later analysis
%
% across cells, sd (session data) struct tracks:
% - TVECc: centers of time bins for data to be fitted [1 x nTimeBins]
% for each model:
% - .m.(modelName).err (nCells x nTimeBins)
% - .m.(modelName).modelspec

p = table;

% common timebase for this session
TVECe = ExpKeys.TimeOnTrack:cfg_master.dt:ExpKeys.TimeOffTrack; % edges
sd.TVECc = TVECe(1:end-1)+cfg_master.dt/2; % centers

% time-to-reward variable is used to restrict samples
t_to_reward = ComputeTTR(reward_t,sd.TVECc);
cfg_ttr = []; cfg_ttr.interp = 'nearest'; cfg_ttr.bins = [-5:0.1:5];
spk_binned = logical(zeros(size(sd.TVECc))); % dummy spike train just to get binned TTR
[ttr, lambda_ttr, ttr_binned] = MakeTC_1D(cfg_ttr, sd.TVECc, t_to_reward, sd.TVECc, spk_binned);

MASTER_keep = ~isnan(ttr); % NEED THIS TO RESTRICT BINNED SPIKE TRAIN LATER!

t_to_reward = t_to_reward(MASTER_keep);
sd.TVECc = sd.TVECc(MASTER_keep);

% baseline model MUST be defined first or things will break!
sd.m.baseline.modelspec = 'spk ~ 1 + linpos + spd + ttr + cif';
%sd.m.gphi.modelspec = 'spk ~ 1 + linpos + spd + ttr + cif + lowGamma_phase + highGamma_phase';
%sd.m.tphi.modelspec = 'spk ~ 1 + linpos + spd + ttr + cif + theta_phase';
sd.m.allphi.modelspec = 'spk ~ 1 + linpos + spd + ttr + cif + delta_phase + beta_phase + theta_phase + lowGamma_phase + highGamma_phase';
sd.m.all.modelspec = 'spk ~ 1 + linpos + spd + ttr + cif + delta_phase + beta_phase + theta_phase + lowGamma_phase + highGamma_phase + delta_env + beta_env + theta_env + lowGamma_env + highGamma_env';

% init error vars
mn = fieldnames(sd.m);
for iM = 1:length(mn)
   sd.m.(mn{iM}).err = zeros(nCells, length(sd.TVECc)); % needs to be zeros because error output will be added to this
end

% define training and testing sets
for iPleat = cfg_master.nPleats:-1:1
    C{iPleat} = cvpartition(length(sd.TVECc), 'KFold', cfg_master.kFold);
end

% LFP features
disp('Computing session-wide LFP features...');
fb_names = fieldnames(cfg_master.f);
cfg_phi.debug = 0;
for iF = 1:length(fb_names)
   
    cfg_phi.fpass = cfg_master.f.(fb_names{iF});
    fprintf('%s: [%d %d]\n', fb_names{iF}, cfg_phi.fpass(1), cfg_phi.fpass(2));
    
    [FF.(fb_names{iF}).phase, FF.(fb_names{iF}).env] = ComputePhase(cfg_phi, csc);
    
end

%
spd = getLinSpd([],pos);

%
linpos = ComputeLinPos(pos);

%%% time / nTrials predictor
p.time = sd.TVECc';

%% loop over all cells 
for iC = nCells:-1:1

    fprintf('Cell %d/%d...\n',iC,nCells);
    
    % dependent variable: binned spike train
    spk_binned = histc(sd.S.t{iC}, TVECe); spk_binned = spk_binned(1:end - 1);
    
    spk_binned = logical(spk_binned > 0); % binarize
    %k = gausskernel(50,2); % smoothing works too, should make into option
    %spk_binned = conv2(spk_binned,k,'same');
    
    %%% PREDICTOR: conditional intensity function
    % note we need to do this before restricting the spike train to avoid weird breaks
    cif = conv(spk_binned,[0 0 0 1 1],'same'); % 2 ms refractory period
    cif = (cif > 0);
    
    %%% PREDICTOR: better CIF based on acorr
    % why does this work so well?
    % what happens for a symmetric acorr?
    % what happens when x-val is done on "blocked" rather than random folds?
    % how do results depend on this thing being present or not?
    cfg_acf = []; cfg_acf.maxlag = 500; cfg_acf.binsize = cfg_master.dt;
    cfg_acf.sided = 'onezero';
    [acf, acf_tvec] = ComputeACF(cfg_acf, spk_binned);
    
    cif_full = conv(spk_binned, acf ./ sum(acf), 'same');
    cif_full = cif_full-nanmean(cif_full);
    cif_full(cif) = 0; % absolute refractory period

    p.cif = cif_full(MASTER_keep);
    
    %%% now can restrict spike train
    spk_binned = spk_binned(MASTER_keep);
    
    %%% PREDICTOR: linpos %%%
    cfg_linpos = [];
    cfg_linpos.bins = linspace(min(linpos.data), max(linpos.data), 100);
    cfg_linpos.interp = 'nearest';
    [predicted_linpos, lambda_linpos, linpos_binned] = MakeTC_1D(cfg_linpos, linpos.tvec, linpos.data, sd.TVECc, spk_binned);
    p.linpos = predicted_linpos' - nanmean(predicted_linpos);
    
    %%% PREDICTOR: speed %%%
    cfg_spd = []; cfg_spd.interp = 'linear';
    cfg_spd.bins = 0:10:200;
    [pred_spd, ~, spd_binned] = MakeTC_1D(cfg_spd, pos.tvec, spd.data, sd.TVECc, spk_binned);
    p.spd = pred_spd' - nanmean(pred_spd);

    %%% PREDICTOR: LFP phase and envelope for all defined bands %%%
    what = {'phase','env'};
    for iF = 1:length(fb_names) % loop over frequency bands
        for iW = 1:length(what)
            this_name = cat(2,fb_names{iF},'_',what{iW});
            switch iW
                case 1 % use phase bins
                    cfg_temp = cfg_phi;
                case 2 % envelope bins depend on data
                    cfg_temp = cfg_phi; cfg_temp.bins = linspace(min(FF.(fb_names{iF}).(what{iW})), 0.5*max(FF.(fb_names{iF}).(what{iW})), 100);
            end
            p.(this_name) = MakeTC_1D(cfg_temp, csc.tvec, FF.(fb_names{iF}).(what{iW}), sd.TVECc, spk_binned)';
        end
    end
    
    %%% PREDICTOR: time to reward
    [ttr, lambda_ttr, ttr_binned] = MakeTC_1D(cfg_ttr, sd.TVECc, t_to_reward, sd.TVECc, spk_binned);
    p.ttr = ttr' - nanmean(ttr);

    %% x-val loop
    p.spk = spk_binned;
    clear cif cif_full;
    
    for iPleat = 1:cfg_master.nPleats
        
        for iFold = 1:C{iPleat}.NumTestSets
            
            fprintf('Pleat %d fold %d...\n', iPleat, iFold);
            
            % get idxs for training and testing set
            tr_idx = C{iPleat}.training(iFold); te_idx = C{iPleat}.test(iFold);
            
            for iModel = 1:length(mn)
                
                % train initial model
                this_m = fitglm(p(tr_idx,:), sd.m.(mn{iModel}).modelspec, 'Distribution', 'binomial')
                sd.m.(mn{iModel}).tstat(iC,:) = this_m.Coefficients.tStat; % should initialize this, but that's a pain
                sd.m.(mn{iModel}).varnames = this_m.PredictorNames;
                
                % refine model by throwing out ineffective predictors
                toss = this_m.Coefficients.Row(abs(this_m.Coefficients.tStat) < 2);
                for iToss = 1:length(toss)
                   this_m.removeTerms(toss{iToss}); 
                end

                % test it and add resulting error to running total
                this_err = this_m.predict(p(te_idx,:));
                this_err = (this_err - spk_binned(te_idx)).^2;
                sd.m.(mn{iModel}).err(iC,te_idx) = sd.m.(mn{iModel}).err(iC,te_idx) + (this_err ./ cfg_master.nPleats)';
                
            end % across models
            
        end % over folds
        
    end % over pleats
    
end % over cells

%% analyze models
cfg_master.smooth = 501; % smoothing window (in samples) for error

[~, lambda_spd, ~] = MakeTC_1D(cfg_ttr, sd.TVECc, t_to_reward, sd.TVECc, spd_binned); % speed by time to reward

for iM = 1:length(mn)
    
    if strcmp(mn{iM},'baseline')
        continue;
    end
    
    figure;
    
    this_err = sd.m.baseline.err - sd.m.(mn{iM}).err;
    
    % mean error over cells
    celldiffmean = nanmean(this_err,2);
    
    subplot(221);
    plot(celldiffmean,'LineWidth',1); hold on;
    keep = find(celldiffmean > 0); plot(keep,celldiffmean(keep),'.g','MarkerSize',20);
    keep = find(celldiffmean < 0); plot(keep,celldiffmean(keep),'.r','MarkerSize',20);
    set(gca,'LineWidth',1,'FontSize',18); box off;
    ylabel('Prediction improvement'); xlabel('cell #');
    title(mn{iM});
    
    % as a tuning curve
    clear lambda_mttr;
    this_f = ones(cfg_master.smooth, 1) ./ cfg_master.smooth;
    mdiffmean = conv(nanmean(this_err), this_f, 'same');
    [~, lambda_mttr, ~] = MakeTC_1D(cfg_ttr, sd.TVECc, t_to_reward, sd.TVECc, mdiffmean);

    subplot(222)
    [ax h1 h2] = plotyy(cfg_ttr.bins,lambda_mttr,cfg_ttr.bins,lambda_spd);
    hold on;
    set(h1,'LineWidth',2);
    set(gca,'XTick',-5:5,'LineWidth',1,'FontSize',18); box off;
    ylabel('Prediction improvement'); xlabel('time from reward (s)');
    
    % t-stat across neurons -- HOW TO GET THE LABELS FROM THE MODELSPEC?
    subplot(223);
    imagesc(sd.m.(mn{iM}).tstat); caxis([0 10]);
    
    subplot(224);
    imagesc(corrcoef(sd.m.(mn{iM}).tstat));

    % could show each cell's model improvement across time
    
end

%% prepare and write output
if cfg_master.writeOutput

    sd.linpos = linpos;
    sd.t_to_reward = t_to_reward;
    sd.cfg = cfg_master;
    
    save(cfg_master.output_prefix,'sd'); % should add option to save in specified output dir
    
end
  
end % of top-level function

%% TODO: correlate GLM results with things like spike spectrum, STA spectrum (across cells)

%% TODO: correlate tuning curves with GLM results
% why are the TTR GLM fits so "streaky"?
% is there some correlation between TTR, space, movement onset TCs

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

%%
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

%%
function linpos_tsd = ComputeLinPos(pos)

f = FindFile('*Coord*');
if ~isempty(f)
    load(f);
else
    error('No Coord file found.');
end

Coord_out = [];
Coord_out.coord = Coord;
Coord_out.nPoints = size(Coord,2);
Coord_out.pointDist = pdist(Coord(:,1:2)','euclidean');
Coord_out.standardized = 0;
Coord_out.units = 'px';

linpos_tsd = LinearizePos([],pos,Coord_out);

end

%%
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
