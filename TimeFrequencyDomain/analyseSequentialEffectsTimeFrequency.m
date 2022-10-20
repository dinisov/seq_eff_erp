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
    
    % for plotting PHOT the right way round
%     if ~light_on_dark
%         allPHOTs = -allPHOTs;
%     end
    
    % onion plot of the photodiode data
    if aux_plots
        figure;
        plot(squeeze(sum(allPHOTs,2)),'b');
        hold on
        plot(mean(squeeze(sum(allPHOTs,2)),2),'r');
    end
    
    % in preparation for calculating the nan mean
    allERPs(allERPs == 0) = nan;
    allPHOTs(allPHOTs == 0) = nan;

    % mean across third dimension (ERP or PHOT stack)
    meanERPs = mean(allERPs, 3, 'omitnan');
    meanPHOTs = mean(allPHOTs, 3, 'omitnan');

    %number of ERPs for each of 32 sequence (also number of PHOTs)
    nERPs = sum(~isnan(allERPs(1,:,:)), 3);
    
    % in order to calculate weighted mean (broadcasting here)
    meanERPs = meanERPs .* nERPs;
    meanPHOTs = meanPHOTs .* nERPs;
    
    %group sequences 2 by 2 (00001 is the same as 11110 and so on)
    meanERPs = meanERPs + fliplr(meanERPs);
    meanPHOTs = meanPHOTs + fliplr(meanPHOTs);
    
    nERPs = nERPs + fliplr(nERPs);
    nERPs(:,n_seq/2 + 1:end) = [];
    
    meanERPs(:,n_seq/2 + 1:end) = []; 
    meanPHOTs(:,n_seq/2 + 1:end) = [];
    meanERPs = meanERPs ./ nERPs;
    meanPHOTs = meanPHOTs ./ nERPs;
    
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
    
end