function R = analyseSequentialEffectsTimeFrequency(blocks, aux_plots)

    %% sort data according to previous sequence

    n_back = 5; n_seq = 2^n_back; %TODO: make n_back variable
    
    n_blocks = length(blocks);
    
    blocks = sortSEs(blocks, n_back, aux_plots);
    
    window = blocks(1).window;%WARNING: we can only merge blocks with windows of same size

    %% sequential effects analysis across all blocks

    total_length = 0;

    for b = 1:n_blocks
        total_length = total_length + size(blocks(b).ERPS,3);
    end 
    
    % concatenate ERPs from different "experiments" (blocks)
    allERPs = zeros(length(window(1):window(2)), n_seq, total_length);
    allPHOTs = zeros(length(window(1):window(2)), n_seq, total_length);
    
    start_index = 1; goodTrials = [];

    % these are already separated by sequence so in order to group by block
    % it is only necessary to stack along third dimension
    for b = 1:n_blocks

        allERPs(:,:,start_index + 1:start_index + size(blocks(b).ERPS,3)) = blocks(b).ERPS;
        allPHOTs(:,:,start_index+ 1:start_index + size(blocks(b).ERPS,3)) = blocks(b).seqPHOT;

        start_index = start_index + size(blocks(b).ERPS,3);

        goodTrials = [goodTrials 1-blocks(b).badTrials]; %#ok<AGROW> 

    end
    
    %%

    % get rid of bad trials (trials with too long gaps between peaks)
    allERPs = allERPs(:,:,logical(goodTrials));
    allPHOTs = allPHOTs(:,:,logical(goodTrials));
    
    %% remove all ERPs with NaNs (from "shaving")
    allERPs = allERPs(:,:,~isnan(sum(squeeze(sum(allERPs,2)))));
    
    % onion plot showing all ERPs for this block
    all_erps = squeeze(sum(allERPs,2));
    if aux_plots
        figure;
        plot(all_erps,'b');
        hold on;
        plot(mean(all_erps,2), 'r');
    end
    
    % onion plot of the photodiode data
    if aux_plots
        figure;
        plot(squeeze(sum(allPHOTs,2)),'b');
        hold on
        plot(mean(squeeze(sum(allPHOTs,2)),2),'r');
    end
    
    % join ERPs corresponding to the same pattern (01001 and 10110 and so on)
    allERPs = allERPs + flip(allERPs,2);
    allPHOTs = allPHOTs + flip(allPHOTs,2);
    allERPs = allERPs(:,1:16,:);
    allPHOTs = allPHOTs(:,1:16,:);
    
    % in preparation for calculating the nan mean
    allERPs(allERPs == 0) = nan;
    allPHOTs(allPHOTs == 0) = nan;

    % mean across third dimension (ERP or PHOT stack)
    meanERPs = mean(allERPs, 3, 'omitnan');
    meanPHOTs = mean(allPHOTs, 3, 'omitnan');

    %number of ERPs for each of 16 sequence (equal to number of PHOTs)
    nERPs = sum(~isnan(allERPs(1,:,:)), 3);
    
    % reorder according to the literature
    meanERPs = meanERPs(:,seq_eff_order(n_back));
    meanPHOTs = meanPHOTs(:,seq_eff_order(n_back));
%     semERPs = semERPs(:,seq_eff_order(n_back));

%     figure;
%     plot(meanERPs)

    %standard errors of the mean for the maxima (use of linear indexing here)
%     semMax = semERPs(sub2ind(size(semERPs),ind_max_erp,1:16));
    
    % plot ERPs for each sequence separately in a 4x4 plot
    % for each sequence, highlight where the maxima (red) and minima (blue) are located 
%     figure;
% 
%     load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');
% 
%     %this is just to help turn horizontal sequences into vertical ones
%     ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]);
% 
%     for i = 1:16
%        subplot(4,4,i);
%        plot(normalize(meanERPs(:,i)));
%        hold on;
%        scatter(ind_max_erp(i), 0,40,'r','filled');
%        scatter(ind_min_erp(i), 0,40,'b','filled');
%        plot(normalize(meanPHOTs(:,i)));
%        plot([stim_onset(i) stim_onset(i)],ylim,'r');
%        title(binomial_x_labels_latex{i}(ind_horiz));
%     end
        
%     L = size(meanERPs,1);

%     fftERPs = fft(meanERPs,10*L);

    cwtERPs = [];
    
    for s = 1:16
%         figure; 
%         [x,f,t] = stft(meanERPs(:,s),3000,'FFTLength',512);
%         cwtERPs(:,s,:) = x; %#ok<AGROW>
%         cwtERPs = stft(meanERPs(:,s),3000);

       [wt,f] = cwt(meanERPs(:,s),3000,'FrequencyLimits',[0 300]);
       cwtERPs(:,s,:) = wt; %#ok<AGROW>

%        cwt(meanERPs(:,s),3000,'FrequencyLimits',[0 100],'VoicesPerOctave',48);  
%        saveas(gcf,['seq_' num2str(s) '.png']);
%        close;
    end
    
    magnitudeSEs = abs(cwtERPs);
    phaseSEs = angle(cwtERPs);
    
    % put all results into a neat structure
    R = struct;
    
    R.magnitudeSEs = magnitudeSEs;
    R.phaseSEs = phaseSEs;
    
    R.meanERPs = meanERPs;
    R.meanPHOTs = meanPHOTs;
    R.nERPs = nERPs;
    
%     R.t = t;
    R.f = f;
    
    time_bounds = [-blocks(1).time_before_peak, blocks(1).time_after_peak];
    
    %make a tick every x steps
    y_ticks = 1:2:length(f);
    y_tick_labels = f(y_ticks);
    
    % create 
    figure; imagesc(squeeze(magnitudeSEs(:,16,:))-squeeze(magnitudeSEs(:,8,:)),'xdata',time_bounds); colorbar;
    set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
    xlabel('time (ms)'); ylabel('Frequency (Hz)');
    title('AAAA_minus_AAAR');

    figure; imagesc(squeeze(magnitudeSEs(:,9,:))-squeeze(magnitudeSEs(:,1,:)),'xdata',time_bounds); colorbar;
    set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
    xlabel('time (ms)'); ylabel('Frequency (Hz)');
    title('RRRR_minus_RRRA');
    
    
end