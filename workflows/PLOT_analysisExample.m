%% analysis example
clear
restoredefaultpath;
%addpath(genpath('D:\My_Documents\GitHub\striatal-spike-rhythms\shared'));
%addpath(genpath('D:\My_Documents\GitHub\striatal-spike-rhythms\chronux_2_12\spectral_analysis'));
addpath(genpath('C:\Users\mvdm\Documents\GitHub\striatal-spike-rhythms\shared'));
addpath(genpath('C:\Users\mvdm\Documents\GitHub\striatal-spike-rhythms\chronux_2_12\spectral_analysis'));

cd('C:\Data\R117\R117-2007-06-20');

LoadExpKeys;

if isfield(ExpKeys,'goodGamma_vStr')
    cfg = []; cfg.fc = ExpKeys.goodGamma_vStr;
elseif isfield(ExpKeys, 'goodGamma')
    cfg = []; cfg.fc = ExpKeys.goodGamma;
else
    error('Couldn''t find LFP field name.');
end

csc = LoadCSC(cfg); csc.data = csc.data-nanmean(csc.data); % could locdetrend to improve STA estimate

cfg = []; cfg.uint = '32';
S = LoadSpikes(cfg); S = SelectTS([], S, 1:2);

%% basic plotting
cols = {[0.5 0.5 0.7], [0.5 0.7 0.5], [0.7 0.5 0.5]};

t = [545.1 545.6];
cscR = restrict(csc, t(1), t(2));
SR = S; SR.t{1} = 1000;
SR.t{2} = [545.1852 545.3051 545.484]';

ax = subplot(5, 5, 1:3);

cfg_mr = [];
cfg_mr.openNewFig = 0;
cfg_mr.axislabel = 'off';
cfg_mr.lfp = cscR;
cfg_mr.SpikeHeight = 1;
MultiRaster(cfg_mr, SR); axis off;

ws = 0.1; yl = ylim;
for iS = 1:length(SR.t{2})
    
    rh(iS) = rectangle('Position', [SR.t{2}(iS)-ws/2 yl(1) ws yl(2)-yl(1)-3]);
    set(rh(iS), 'EdgeColor', cols{iS}, 'FaceColor', cols{iS});
end

MultiRaster(cfg_mr, SR); hold on; axis off;
for iS = 1:length(SR.t{2})
    lh(iS) = plot([SR.t{2}(iS) SR.t{2}(iS)], [yl(1) yl(2)-3], '--', 'Color', [0.3 0.3 0.3]);
end
%% raw STs
Fs = 1 ./ mode(diff(csc.tvec));
tvec = -ws/2: 1./Fs : ws/2; % time axis for STA
 
spk_t = SR.t{2};
for iSpk = length(spk_t):-1:1 % for each spike...
 
   sta_t = spk_t(iSpk) - ws/2;
   sta_idx = nearest_idx3(sta_t, cscR.tvec); % find index of leading window edge
   toAdd = cscR.data(sta_idx:sta_idx+length(tvec)-1); % grab LFP snippet for this window
 
   sta(iSpk,:) = toAdd'; 
 
   subplot(5, 5, 6 + (iSpk-1)*5);
   plot(tvec, toAdd', 'k'); set(gca,'Color', cols{iSpk}, 'XColor', cols{iSpk}, 'YColor', cols{iSpk}, 'XTick', [], 'YTick', [], 'FontSize', 12);
   ylim([-0.03 0.03]*10e-3); ylabel(sprintf('spike %d LFP', iSpk)); xlabel('time'); title(' ');
   
   hold on;
   lh(iS) = plot([0 0], [-0.03 0.03]*10e-3, '--', 'Color', [0 0 0]);
   
end

%% summed STA
sta(length(spk_t) + 1, :) = nanmean(sta);

subplot(5, 5, 21);
plot(tvec, sta(end, :), 'k');
ylim([-0.03 0.03]*10e-3);
hold on;
lh(iS) = plot([0 0], [-0.03 0.03]*10e-3, '--', 'Color', [0 0 0]);

set(gca,'Color', [0.7 0.7 0.7], 'XColor', [0.7 0.7 0.7], 'YColor', [0.7 0.7 0.7], 'XTick', [], 'YTick', [], 'FontSize', 12);
ylabel(sprintf('average LFP', iSpk)); xlabel('time'); title('STA');

%% pow & phase
nFFT = 1024; nP = size(sta, 2);
for iS = 1:size(sta, 1)
    
    this_fft = fft(sta(iS,:)-nanmean(sta(iS,:)), nFFT) / nP;
    FFTf = (Fs/2) * linspace(0, 1, (nFFT/2) + 1);
    
    pow(iS, :) = abs(this_fft(1:nFFT/2+1));
    phi(iS, :) = angle(this_fft(1:nFFT/2+1));

    subplot(5, 5, 7 + (iS-1)*5) 
    if iS ~= size(sta, 1) % plot phases for 3 STs
        plot(FFTf, phi(iS,:), '.', 'Color', cols{iSpk}); 
        set(gca, 'YLim', [-pi pi], 'YTick', [-pi:pi/2:pi], 'YTickLabel', {'-pi', '', '0', '', 'pi'});
       
    else
        plot(FFTf, pow(iS,:), 'k'); 
    end
    box off; set(gca, 'XLim', [0 100], 'XTick', [0:20:100], 'FontSize', 14, 'TickDir', 'out', 'LineWidth', 1);
    grid on;
    
    
end
