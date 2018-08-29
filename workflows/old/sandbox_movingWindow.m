%%
% remember to set path

%addpath(genpath('D:\My_Documents\GitHub\striatal-spike-rhythms\chronux_2_12\spectral_analysis'));

%fd = {'D:\data\adrlab\R117\R117-2007-06-20'};
fd = {'C:\data\adrlab\R117-2007-06-20'};
cd(fd{1});

LoadExpKeys;
cfg = []; cfg.fc = ExpKeys.goodGamma_vStr;
csc = LoadCSC(cfg); csc.data = csc.data-nanmean(csc.data); % may need to locdetrend
S = LoadSpikes([]);

%% Fourier transform of xcorr -- way faster!
cfg = [];
% no binsize because this is taken from csc.tvec
%cfg.maxlag = 10000;
cfg.STAwidth = 1000; % currenly in samples; Fs is about 2000
cfg.nShuf = 10;
cfg.wSize = 60; % moving window size (s)
cfg.wStep = 30;
%cfg.gauss_w = 1; cfg.gauss_sd = 0.002;

%gauss_window = cfg.gauss_w./cfg.binsize; % 1 second window
%gauss_SD = cfg.gauss_sd./cfg.binsize; % 0.02 seconds (20ms) SD
%gk = gausskernel(gauss_window,gauss_SD); 
%gk = gk.*cfg.binsize; % normalize by binsize

%tbin_edges = firstSpike(S):cfg.binsize:lastSpike(S);

for iC = 11;%length(S.t):-1:1
    
    fprintf('Cell %d/%d...\n',iC,length(S.t));
    
    spk_t = S.t{iC};
    
    % binarize spike train on csc.tvec timebase
    spk_binned = zeros(size(csc.tvec));
    idx = nearest_idx3(spk_t, csc.tvec);
    spk_binned(idx) = 1;
    cfg.binsize = median(diff(csc.tvec)); cfg.Fs = 1 ./ cfg.binsize;
    
    % could convolve spike train here for acorr computation later
    
    % define moving window start and end points
    wSizeBins = round(cfg.wSize ./ cfg.binsize); wStepBins = round(cfg.wStep ./ cfg.binsize);
    winStart_idx = 1:wStepBins:length(spk_binned) - wSizeBins;
    winEnd_idx = wSizeBins:wStepBins:length(spk_binned);
    
    if length(winStart_idx) ~= length(winEnd_idx)
        winStart_idx = winStart_idx(1:end - 1);
    end
    
    % loop over moving window time bins
    for iW = length(winStart_idx):-1:1
        
        this_spk_binned = spk_binned(winStart_idx(iW):winEnd_idx(iW));
        
        % xcorr
        [acf,acf_tvec] = ComputeACF(cfg,this_spk_binned);
        [P,F_acf] = pwelch(acf, 2^12, 2^11, [], cfg.Fs);
        
        P_all(iC,iW,:) = P; 
               
        % STA    
        [STA,tvec] = xcorr(csc.data, this_spk_binned, cfg.STAwidth);
        STA_tvec = tvec ./ cfg.Fs;
        
        STA_all(iC,iW,:) = STA;
        
        [STAp_all(iC,iW,:),F_STA] = pwelch(STA, length(STA), [], [], cfg.Fs);
                
        % shuffles
        clear this_shufP this_shufSTA this_shufSTAp;
        for iShuf = cfg.nShuf:-1:1
            spk_binned_shuf = this_spk_binned(randperm(length(this_spk_binned)));
            
            acf = ComputeACF(cfg,spk_binned_shuf);
            this_shufP(iShuf,:) = pwelch(acf, 2^12, 2^11, [], cfg.Fs);
            
            this_shufSTA(iShuf,:) = xcorr(csc.data, spk_binned_shuf, cfg.STAwidth);
            this_shufSTAp(iShuf,:) = pwelch(this_shufSTA(iShuf,:), length(this_shufSTA(iShuf,:)), [], [], cfg.Fs);
        end
        
        P_shufmean(iC,iW,:) = nanmean(this_shufP);
        P_shufSD(iC,iW,:) = nanstd(this_shufP);
        
        STA_shufmean(iC,iW,:) = nanmean(this_shufSTA); STAp_shufmean(iC,iW,:) = nanmean(this_shufSTAp);
        STA_shufSD(iC,iW,:) = nanstd(this_shufSTA); STAp_shufSD(iC,iW,:) = nanstd(this_shufSTAp);
        
    end

end

%% plot some basic things

imagesc(F_STA,1:103,10*log10(sq(STAp_all(11,:,:)))); xlim([0 100]);

% need z-score version

% how does this covary with LFP power?

%%
function [acf,tvec] = ComputeACF(cfg,spk_binned)
%[acf,tvec] = xcorr(spk_binned,spk_binned,cfg.maxlag);
[acf,tvec] = xcorr(spk_binned,spk_binned);
tvec = tvec.*cfg.binsize;
acf(ceil(length(acf)/2)) = 0;
end