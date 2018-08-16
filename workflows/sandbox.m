%%
% remember to set path

addpath(genpath('D:\My_Documents\GitHub\striatal-spike-rhythms\chronux_2_12\spectral_analysis'));

fd = {'D:\data\adrlab\R117\R117-2007-06-20'};
cd(fd{1});

LoadExpKeys;
cfg = []; cfg.fc = ExpKeys.goodGamma_vStr;
csc = LoadCSC(cfg);

%%

S = LoadSpikes([]);

t = 1:0.1:1000;
t = t + 0.01*randn(size(t));
t = t';

t = S.t{1};
%% chronux spike spectrum -- how can we control the number of frequency bins? seems unwieldy!
tic
params = []; 
params.Fs = 1000;
params.fpass = [0 100];
params.tapers = [3 5];
params.pad = 0;
[P,F,R] = mtspectrumpt(t,params);
subplot(221);
plot(F,10*log10(P));
title(toc)

%% Fourier transform of xcorr -- way faster!
cfg = [];
cfg.binsize = 0.001; cfg.Fs = 1./cfg.binsize;
cfg.nShuf = 1;

tbin_edges = firstSpike(S):cfg.binsize:lastSpike(S);

for iC = length(S.t):-1:1
    
    fprintf('Cell %d/%d...\n',iC,length(S.t));
    
    spk_t = S.t{iC};
    
    % binarize spike trains
    spk_binned = histc(spk_t,tbin_edges); spk_binned = spk_binned(1:end-1);
    
    % MATLAB xcorr
    [acf,tvec] = ComputeACF(cfg,spk_binned);
    figure;
    subplot(221);
    plot(tvec,acf); xlim([-0.5 0.5]);
    
    % transform
    [P,F] = pwelch(acf,2^12,2^11,[],cfg.Fs); P(1:5) = NaN;
    
    subplot(223);
    plot(F,10*log10(P)); xlim([0 100]);
    title('acorr spectrum');
    
    bigP(iC,:) = P;
    
    % shuffles
    clear this_shufP;
    for iShuf = cfg.nShuf:-1:1
        spk_binned_shuf = spk_binned(randperm(length(spk_binned)));
        [acf,tvec] = ComputeACF(cfg,spk_binned_shuf);
        this_shufP(iShuf,:) = pwelch(acf,2^12,2^11,[],cfg.Fs); this_shufP(iShuf,1:5) = NaN;
    end
    
    bigP_shufmean(iC,:) = nanmean(this_shufP);
    bigP_shufSD(iC,:) = nanstd(this_shufP);
    
    hold on;
    plot(F,10*log10(bigP_shufmean(iC,:)),'r');
    plot(F,10*log10(bigP_shufmean(iC,:)+2*bigP_shufSD(iC,:)),'r:');
    plot(F,10*log10(bigP_shufmean(iC,:)-2*bigP_shufSD(iC,:)),'r:');
    drawnow;
    
    % get STA by binning on CSC timebase
    this_spk_binned = zeros(size(csc.tvec));
    idx = nearest_idx3(spk_t,csc.tvec);
    this_spk_binned(idx) = 1;
    [xc,tvec] = xcorr(csc.data,this_spk_binned,1000);
    tvec = tvec.*median(diff(csc.tvec));
    
    subplot(222); 
    plot(tvec,xc);
    title('raw STA');
    
    % transform
    %[xc,tvec] = xcorr(csc.data,this_spk_binned);
    this_Fs = 1./median(diff(csc.tvec));
    [P,F] = pwelch(xc,length(xc),[],[],this_Fs);
    
    subplot(224);
    plot(F,10*log10(P)); xlim([0 100]);
    title('STA spectrum');
    
    drawnow;
end

%%
function [acf,tvec] = ComputeACF(cfg,spk_binned)
[acf,tvec] = xcorr(spk_binned,spk_binned);
tvec = tvec.*cfg.binsize;
acf(ceil(length(acf)/2)) = 0;
end