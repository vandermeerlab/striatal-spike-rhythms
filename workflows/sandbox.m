%%
% remember to set path

%addpath(genpath('D:\My_Documents\GitHub\striatal-spike-rhythms\chronux_2_12\spectral_analysis'));

%fd = {'D:\data\adrlab\R117\R117-2007-06-20'};
fd = {'C:\data\adrlab\R117-2007-06-20'};
cd(fd{1});

LoadExpKeys;
cfg = []; cfg.fc = ExpKeys.goodGamma_vStr;
csc = LoadCSC(cfg); csc.data = csc.data-nanmean(csc.data); % may need to locdetrend

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
%cfg.maxlag = 10000;
cfg.nShuf = 10;
%cfg.gauss_w = 1; cfg.gauss_sd = 0.002;

%gauss_window = cfg.gauss_w./cfg.binsize; % 1 second window
%gauss_SD = cfg.gauss_sd./cfg.binsize; % 0.02 seconds (20ms) SD
%gk = gausskernel(gauss_window,gauss_SD); 
%gk = gk.*cfg.binsize; % normalize by binsize

tbin_edges = firstSpike(S):cfg.binsize:lastSpike(S);

for iC = length(S.t):-1:1
    
    fprintf('Cell %d/%d...\n',iC,length(S.t));
    
    spk_t = S.t{iC};
    
    % binarize spike trains
    spk_binned = histc(spk_t,tbin_edges); spk_binned = spk_binned(1:end-1);
    % could convolve here
    %spk_binnedC = conv2(spk_binned,gk,'same');
    spk_binnedC = spk_binned;
    
    % MATLAB xcorr
    [acf,tvec] = ComputeACF(cfg,spk_binnedC);
    figure;
    subplot(221);
    plot(tvec,acf); xlim([-0.5 0.5]);
    title(sprintf('Cell %d',iC));
    
    % transform
    %[P,F] = pwelch(acf(1:floor(length(acf)/2)),2^11,2^10,[],cfg.Fs); P(1:5) = NaN;
    [P,F] = pwelch(acf,2^12,2^11,[],cfg.Fs); P(1:5) = NaN;
    
    subplot(223);
    plot(F,10*log10(P)); xlim([0 100]);
    title('acorr spectrum');
    
    bigP(iC,:) = P; F_acf = F;
    
    % shuffles
    clear this_shufP;
    for iShuf = cfg.nShuf:-1:1
        spk_binned_shuf = spk_binned(randperm(length(spk_binned)));
        %spk_binned_shuf = conv2(spk_binned_shuf,gk,'same');
        [acf,tvec] = ComputeACF(cfg,spk_binned_shuf);
        %this_shufP(iShuf,:) = pwelch(acf(1:floor(length(acf)/2)),2^11,2^10,[],cfg.Fs); this_shufP(iShuf,1:5) = NaN;
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
    
    STA_all(iC,:) = xc;
    
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
    
    STAp_all(iC,:) = P;
    
    % do STA shuffles
    clear this_shufSTA this_shufSTAp;
    for iShuf = cfg.nShuf:-1:1
        spk_binned_shuf = this_spk_binned(randperm(length(this_spk_binned)));
        this_shufSTA(iShuf,:) = xcorr(csc.data,spk_binned_shuf,1000);
        this_shufSTAp(iShuf,:) = pwelch(this_shufSTA(iShuf,:),length(this_shufSTA(iShuf,:)),[],[],this_Fs);
    end
    STA_shufmean(iC,:) = nanmean(this_shufSTA); STAp_shufmean(iC,:) = nanmean(this_shufSTAp);
    STA_shufSD(iC,:) = nanstd(this_shufSTA); STAp_shufSD(iC,:) = nanstd(this_shufSTAp);
    
    subplot(222); hold on;
    plot(tvec,STA_shufmean(iC,:),'r');
    plot(tvec,STA_shufmean(iC,:)+2*STA_shufSD(iC,:),'r:');
    plot(tvec,STA_shufmean(iC,:)-2*STA_shufSD(iC,:),'r:');
    
    subplot(224); hold on;
    plot(F,10*log10(STAp_shufmean(iC,:)),'r');
    plot(F,10*log10(STAp_shufmean(iC,:)+2*STAp_shufSD(iC,:)),'r:');
    plot(F,10*log10(STAp_shufmean(iC,:)-2*STAp_shufSD(iC,:)),'r:');

    drawnow;
end

%% compute some z-scores and plot
bigPz = (bigP-bigP_shufmean)./bigP_shufSD;
bigPzi = interp1(F_acf,bigPz',1:100);
subplot(221); imagesc(F_acf,1:27,bigPz); caxis([-10 10]);
subplot(222); imagesc(1:100,1:27,bigPzi'); caxis([-10 10]);

STApz = (STAp_all-STAp_shufmean)./STAp_shufSD;
STApzi = interp1(F,STApz',1:100);

subplot(223); imagesc(F,1:27,STApz); caxis([-10 10]);
subplot(224); imagesc(1:100,1:27,STApzi'); caxis([-10 10]);

for iF = 1:100
 
    [temp1,temp2] = corrcoef(STApzi(iF,:),bigPzi(iF,:));
    cc(iF) = temp1(1,2);
    pval(iF) = temp2(1,2);
    
end
figure; plot(cc); set(gca,'XTick',0:5:100); grid on;
xval = 1:100;
keep = pval < 0.05;
hold on; plot(xval(keep),cc(keep),'.b','MarkerSIze',20);
%%
function [acf,tvec] = ComputeACF(cfg,spk_binned)
%[acf,tvec] = xcorr(spk_binned,spk_binned,cfg.maxlag);
[acf,tvec] = xcorr(spk_binned,spk_binned);
tvec = tvec.*cfg.binsize;
acf(ceil(length(acf)/2)) = 0;
end