%% quick plot for single cells
for iC = 1:cc

    % session-specgram
    subplot(331); 
    imagesc(ALL.sessionTFR(iC).time, ALL.sessionTFR(iC).freq, sq(ALL.sessionTFR(iC).powspctrm)); axis xy; colorbar;
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out'); xlabel('time'); ylabel('frequency (Hz)'); title('LFP spectrogram');
    
    % spike spectrum -- COULD SMOOTH THIS SOME
    subplot(332);
    plot(ALL.spkSpec_freq, medfilt1(ALL.spkSpec(iC, :),25,'truncate'), 'k', 'LineWidth', 2);
    hold on;
    plot(ALL.spkSpec_freq, medfilt1(ALL.spkSpec_shufmean(iC, :),25) + medfilt1(ALL.spkSpec_shufSD(iC, :),25), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    plot(ALL.spkSpec_freq, medfilt1(ALL.spkSpec_shufmean(iC, :),25) - medfilt1(ALL.spkSpec_shufSD(iC, :),25), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out'); xlabel('frequency (Hz)'); ylabel('power'); title('spike spectrum');
    
    % STA (raw)
    subplot(333);
    plot(ALL.STAtime, ALL.STA(iC, :), 'k', 'LineWidth', 1);
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out','XLim',[-0.5 0.5],'XTick',-0.5:0.1:0.5); xlabel('time (s)'); ylabel('LFP'); title('STA');
    set(gca,'XTickLabel',{'-0.5','','','','','0','','','','','0.5'});
    grid on;
    
    subplot(334);
    plot(ALL.ppc_freq, ALL.ppc(iC, :), 'k', 'LineWidth', 2);
    hold on;
    plot(ALL.ppc_freq, ALL.ppc_shufmean(iC, :) + ALL.ppc_shufSD(iC, :), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    plot(ALL.ppc_freq, ALL.ppc_shufmean(iC, :) - ALL.ppc_shufSD(iC, :), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out'); xlabel('frequency (Hz)'); ylabel('ppc'); title('spike-field ppc');
    % try to recover this average from ppc-gram and spike count?
    
    subplot(335);
    imagesc(ALL.trPPCtime, ALL.trPPCfreq, sq(ALL.trPPC(iC,:,:))); axis xy; colorbar;
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out','XLim',[-5 5]); xlabel('time'); ylabel('frequency (Hz)'); title('raw PPC');
    
    subplot(336);
    plot(ALL.STAp_freq,ALL.STAp(iC,:), 'k', 'LineWidth', 2);
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out'); xlabel('frequency (Hz)'); ylabel('power'); title('STA power');
    
    subplot(337); bar(ALL.trPPCtime, sq(ALL.trPPCn(iC,1,:))); box off;
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out','XLim',[-5 5]); xlabel('time'); ylabel('count'); title('spike count histogram');
    
    subplot(338)
    imagesc(ALL.tr_time-5.5, ALL.tr_freq, sq(ALL.trP(iC,:,:))'); axis xy; colorbar;
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out','XLim',[-5 5]); xlabel('time'); ylabel('frequency (Hz)'); title('spike spectrum (raw)');
    % how to normalize this by spike count? need histogram on chronux timebase...
    
    this_trp = sq(ALL.trP(iC,:,:))';
    sc = repmat(ALL.trP_hist(iC,:), [size(this_trp,1) 1]);
    this_trp = this_trp ./ sc;
    
    subplot(339)
    imagesc(ALL.tr_time-5.5, ALL.tr_freq, this_trp); axis xy; colorbar;
    set(gca,'LineWidth',1,'FontSize',18,'TickDir','out','XLim',[-5 5]); xlabel('time'); ylabel('frequency (Hz)'); title('spike spectrum (n)');
    
    drawnow;
    pause; clf;

end


%% some averages
what = {'FSI', 'MSNonly', 'nonFSI', 'all'};

for iW = 1:length(what)
    
    switch what{iW} % select appropriate cell idxs
        case 'FSI'
            
        case 'MSNonly'
            
        case 'nonFSI'
            
        case 'all'
    
    end
    
    % assemble data
    
end

