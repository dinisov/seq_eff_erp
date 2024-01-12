function R = timeFrequencySpectrum(R, blocks)
%timeFrequencyAnalysis Summary of this function goes here
%   Detailed explanation goes here
    cwtERPs = [];
    
    %300 Hz seems like a sensible upper limit
    for s = 1:16
       [wt,f] = cwt(R.meanERPs(:,s),blocks(1).resampleFreq,'FrequencyLimits',[0 300]);
       cwtERPs(:,s,:) = wt; %#ok<AGROW>
    end
    
    magnitudeSEs = abs(cwtERPs);
    phaseSEs = angle(cwtERPs);
     
    R.magnitudeSEs = magnitudeSEs;
    R.phaseSEs = phaseSEs;
    R.f = f;
    
    % for block experiments, this is set in seq_eff_erp() as [0 InterBlockPeriod]
    time_bounds = blocks(1).window*1000;
    
    %make a tick every x steps
    y_ticks = 1:2:length(f);
    y_tick_labels = f(y_ticks);
    
    % spectrogram for AAAA minus AAAR
    figure; imagesc(squeeze(magnitudeSEs(:,16,:)) - squeeze(magnitudeSEs(:,8,:)),'xdata',time_bounds); colorbar;
    set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
    xlabel('time (ms)'); ylabel('Frequency (Hz)');
    title('AAAA minus AAAR');

    % spectrogram for RRRR minus RRRA
    figure; imagesc(squeeze(magnitudeSEs(:,9,:)) - squeeze(magnitudeSEs(:,1,:)),'xdata',time_bounds); colorbar;
    set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
    xlabel('time (ms)'); ylabel('Frequency (Hz)');
    title('RRRR minus RRRA');
end