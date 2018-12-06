%% script exploring spike autocorrelations, spike spectra, and spike-triggered averages
% remember to set path
clear
addpath(genpath('/Users/jericcarmichael/Documents/GitHub/EC_Multisite'))
addpath(genpath('/Users/jericcarmichael/Documents/GitHub/striatal-spike-rhythms'))%chronux_2_12\spectral_analysis'));
% fd = {'D:\data\adrlab\R117\R117-2007-06-20'};
% fd = {'C:\data\adrlab\R117-2007-06-20'};
% fd = {'C:\data\R042\R042-2013-08-17'};
fd = {'/Volumes/Fenrir/Str_rhythms/R117-2007-06-20'}; % EC dta location
cd(fd{1});

LoadExpKeys;
cfg = []; cfg.fc = ExpKeys.goodGamma_vStr;
% cfg = []; cfg.fc = ExpKeys.goodTheta;
csc = LoadCSC(cfg); csc.data = csc.data-nanmean(csc.data); % may need to locdetrend

%% load data
S = LoadSpikes([]);

%% artificial "gamma cell" which is firing at 55 Hz
S.label{end+1} = 'artificial';
S.t{end+1} = [1:(csc.cfg.hdr{1}.SamplingFrequency/55)*0.001:length(csc.tvec)/4]';
CheckTS(S)
%% detect gamma events
cfg_gamma_detect = [];
cfg_gamma_detect.detect_thr = [0.7 .7 ]; % threshold for event detection: 95th percentile of (amplitude) envelope

gamma_iv = MS_DetectEvents(cfg_gamma_detect, csc, ExpKeys);

%% detect NON-gamma events

cfg_other = [];
cfg_other.epoch = 'all'; % chewing occurs during task mostly
cfg_other.minlen = 0.005;
cfg_other.filter_cfg.f = [45 65]; % default [25 35]
cfg_other.threshold = .5; % 0.25 for session 2, 0.5 for session 1?
cfg_other.smooth = 0.05; % convolve with Gaussian of this SD
cfg_other.dcn = '<';
non_gamma_iv = DetectEvents(cfg_other, csc, ExpKeys);

non_gamma_iv.usr.evt_len = (non_gamma_iv.tend - non_gamma_iv.tstart)';

% remove long events
cfg_max_len = [];
cfg_max_len.operation = '<';
cfg_max_len.threshold = max(gamma_iv.low.tend - gamma_iv.low.tstart);
non_gamma_iv = SelectIV(cfg_max_len,non_gamma_iv,'evt_len');

%     remove events that are too short
cfg_max_len = [];
cfg_max_len.operation = '>';
cfg_max_len.threshold = min(gamma_iv.low.tend - gamma_iv.low.tstart);
evt = SelectIV(cfg_max_len,non_gamma_iv,'evt_len');

% select a random assortment of the non-gamma events to get the same number
% of events
keep_idx = datasample(1:length(non_gamma_iv.tend),length(gamma_iv.low.tend)); % get nSamples from gamma events chosen randomly from the number of samples in non gamma
non_gamma_iv = SelectIV([],non_gamma_iv,keep_idx);

% PlotTSDfromIV([], non_gamma_iv, csc)

%% try multiraster
% cfg = [];
% cfg.evt = gamma_iv.low; 
% cfg.lfp = csc;
% cfg.spkColor = 'jet';
% MultiRaster(cfg,S) 

%% get ACF for each gamma event

this_band = 'low';
%
cfg_acf = [];
cfg_acf.binsize = .001; 
cfg_acf.maxlag = 200;
all_acf = cell(1,length(S.t));
for iEvt = 1:500%length(gamma_iv.(this_band).tstart) % loop events
    disp(num2str(iEvt))
    this_event = restrict(S, gamma_iv.(this_band).tstart(iEvt), gamma_iv.(this_band).tend(iEvt));
    % create bins for this gamma event. 
    tbin_edges = gamma_iv.(this_band).tstart(iEvt):cfg_acf.binsize:gamma_iv.(this_band).tend(iEvt); % vector of time bin edges (for histogram)
    tbin_centers = tbin_edges(1:end-1)+cfg_acf.binsize/2; % vector of time bin centers (for plotting)

    % compute
    for iC = 1:length(S.t) % loop cells
        if ~isempty(this_event.t{iC})
            if length(this_event.t{iC})>1
                % bin the data for this cell. 
                spk_count = histc(S.t{iC},tbin_edges); % get spike counts for each bin
                spk_count = spk_count(1:end-1);
                
             [all_acf{iC}(iEvt,:) ,tvec_acf] = ComputeACF(cfg_acf, spk_count);
             
            end
        end
    end
end

%% same thing for NON-GAMMA epochs
cfg_acf = [];
cfg_acf.binsize = .001; 
cfg_acf.maxlag = 200;
all_acf_non = cell(1,length(S.t));
for iEvt = 1:500%length(non_gamma_iv.tstart) % loop events
    this_event = restrict(S, non_gamma_iv.tstart(iEvt), non_gamma_iv.tend(iEvt));
    % create bins for this gamma event. 
    tbin_edges = non_gamma_iv.tstart(iEvt):cfg_acf.binsize:non_gamma_iv.tend(iEvt); % vector of time bin edges (for histogram)
    tbin_centers = tbin_edges(1:end-1)+cfg_acf.binsize/2; % vector of time bin centers (for plotting)

    % compute
    for iC = 1:length(S.t) % loop cells
        if ~isempty(this_event.t{iC})
            if length(this_event.t{iC})>1
                % bin the data for this cell. 
                spk_count = histc(S.t{iC},tbin_edges); % get spike counts for each bin
                spk_count = spk_count(1:end-1);
                
             [all_acf_non{iC}(iEvt,:) ,tvec] = ComputeACF(cfg_acf, spk_count);
            end
        end
    end
end

%% remove cells that have very few data points
min_points = 10;
cells_to_keep =[];
for iC = 1:length(all_acf)
    if max(sum(all_acf{iC})) > min_points &&  max(sum(all_acf_non{iC})) > min_points
        cells_to_keep =[cells_to_keep iC];
    else
        fprintf('\nCell %.0f had < %.0f number of data points\n', iC,  min_points)
    end
end

%% plot all the gamma events vs non events.  

figure(222) 
x_plot = 1; 
for iC = 1:length(cells_to_keep)
    subplot(4,4,x_plot)

plot(tvec, nanmean(all_acf{cells_to_keep(iC)}),'r', tvec, nanmean(all_acf_non{cells_to_keep(iC)}), 'k')
title(num2str(cells_to_keep(iC)))
x_plot = x_plot+1;
end

%% plots?
figure(100)
plot_n = 0;
Fs = 1./cfg_acf.binsize;
S_gamma = restrict(S, gamma_iv.(this_band).tstart, gamma_iv.(this_band).tend);
S_non = restrict(S, non_gamma_iv.tstart, non_gamma_iv.tend);

%%
figure(100)
plot_n = 0;
for iC =cells_to_keep(1:ceil(length(cells_to_keep)/2))
    
        plot_n = plot_n+1;
        subplot(4,length(cells_to_keep), plot_n)
        hold on
        plot(tvec_acf, sum(all_acf_non{iC}), 'r')
        plot(tvec_acf, sum(all_acf{iC}), 'b')
        xlim([-0.05 0.05])
        ylim([0 inf])
        
        % plot the PSD of the ACORR
         subplot(4, length(cells_to_keep), plot_n+length(cells_to_keep))
         
         [P,F] = pwelch(all_acf{iC}(1:floor(length(all_acf{iC})/2)),2^7,2^6,[],Fs);
         
         plot(F,10*log10(P)); xlim([0 100]);
         title('acorr spectrum');
         
         % get the sta in gamma events
         this_spk_binned = zeros(size(csc.tvec));
         idx = nearest_idx3(S_gamma.t{iC},csc.tvec);
         this_spk_binned(idx) = 1;
         [xc_gamma,tvec] = xcorr(csc.data,this_spk_binned,1000);
         tvec = tvec.*median(diff(csc.tvec));
         % same for non gamma
                  this_spk_binned = zeros(size(csc.tvec));
         idx = nearest_idx3(S_non.t{iC},csc.tvec);
         this_spk_binned(idx) = 1;
         [xc_non,tvec] = xcorr(csc.data,this_spk_binned,1000);
                  
         subplot(4, length(cells_to_keep), plot_n+length(cells_to_keep)*2)

         plot(tvec,xc_non,'r', tvec, xc_gamma, 'b');
         
         title('raw STA');

end
tightfig
maximize
%%
figure(200)
plot_n = 0;
for iC =cells_to_keep(ceil(length(cells_to_keep)/2):end)
    
        plot_n = plot_n+1;
        subplot(4,ceil(length(cells_to_keep)/2), plot_n)
        hold on
        plot(tvec_acf, sum(all_acf_non{iC}), 'r')
        plot(tvec_acf, sum(all_acf{iC}), 'b')
        xlim([-0.05 0.05])
        ylim([0 inf])
        
        % plot the PSD of the XCORR
         subplot(4, ceil(length(cells_to_keep)/2), plot_n+ceil(length(cells_to_keep)/2))
         
         [P,F] = pwelch(all_acf{iC}(1:floor(length(all_acf{iC})/2)),2^10,2^9,[],Fs);
         
         plot(F,10*log10(P)); xlim([0 100]);
         title('acorr spectrum');
         
         % get the sta in gamma events
         this_spk_binned = zeros(size(csc.tvec));
         idx = nearest_idx3(S_gamma.t{iC},csc.tvec);
         this_spk_binned(idx) = 1;
         [xc_gamma,tvec] = xcorr(csc.data,this_spk_binned,1000);
         tvec = tvec.*median(diff(csc.tvec));
         % same for non gamma
                  this_spk_binned = zeros(size(csc.tvec));
         idx = nearest_idx3(S_non.t{iC},csc.tvec);
         this_spk_binned(idx) = 1;
         [xc_non,tvec] = xcorr(csc.data,this_spk_binned,1000);
                  
         subplot(4, ceil(length(cells_to_keep)/2), plot_n+(ceil(length(cells_to_keep)/2))*2)

         plot(tvec,xc_non,'r', tvec, xc_gamma, 'b');
         
         title('raw STA');

end
tightfig
maximize

%% get the psd for each STA across all events.  






