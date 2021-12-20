function R = analyseSequentialEffects(blocks, aux_plots)

    %% sort data according to previous sequence

    n_back = 5; n_seq = 2^n_back; %TODO: make n_back variable
    
    n_blocks = length(blocks);
    
    for b = 1:n_blocks
        
        LOCS = blocks(b).LOCS;
        LFP = blocks(b).LFP;
        PHOT = blocks(b).PHOT;
%         rawPHOT = blocks(b).rawPHOT;
        randomSequence = blocks(b).randomSequence;
        time_before_peak = blocks(b).time_before_peak;
        time_after_peak = blocks(b).time_after_peak;
        resampleFreq = blocks(b).resampleFreq;
        badTrials = blocks(b).badTrials;
    
        window = floor([-time_before_peak*resampleFreq, time_after_peak*resampleFreq]);

        sequenceLength = length(randomSequence);

        ERPS = zeros(length(window(1):window(2)), n_seq, sequenceLength);
        seqPHOT = zeros(length(window(1):window(2)), n_seq, sequenceLength);
        
        for n = n_back:sequenceLength

            % decimal value of binary sequence of length n_back
            seq = bin2dec(num2str(randomSequence(n-n_back+1:n))) + 1;

            % stack ERPs and PHOTs along third dimension (first two dims are sequence and
            % time respectively)
            ERPS(:, seq, n) = LFP(LOCS(n) + window(1) : LOCS(n) + window(2));
            seqPHOT(:, seq, n) = normalize(PHOT(2-randomSequence(n), LOCS(n) + window(1) : LOCS(n) + window(2) ));

        end

        % matrix with all ERPs irrespective of sequence
        all_erps = squeeze(sum(ERPS,2));% squeeze removes first dimension

        %calculate mean and SEM for outlier calculations
        meanERP = mean(all_erps, 2);
        STDs = std(all_erps, [], 2);

        n_sd = 3;

        % remove ERPs beyond n_sd
        outliers = all_erps < (meanERP - n_sd*STDs) | all_erps > (meanERP + n_sd*STDs);

        good_erps = ~logical(sum(outliers));

        % remove ERP outliers
        ERPS = ERPS(:,:,good_erps);
        seqPHOT = seqPHOT(:,:,good_erps);
        badTrials = badTrials(good_erps);
        
        blocks(b).ERPS = ERPS;
        blocks(b).badTrials = badTrials;
        blocks(b).seqPHOT = seqPHOT;
        
        % onion plot showing all ERPs for this block
        all_erps = all_erps(:,good_erps);
        if aux_plots
            figure;
            plot(all_erps,'b');
            hold on;
            plot(mean(all_erps,2), 'r');
        end
    
    end

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

        allERPs(:,:,start_index:start_index + size(blocks(b).ERPS,3) - 1) = blocks(b).ERPS;
        allPHOTs(:,:,start_index:start_index + size(blocks(b).ERPS,3) - 1) = blocks(b).seqPHOT;

        start_index = start_index + size(blocks(b).ERPS,3);

        goodTrials = [goodTrials 1-badTrials];

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
    
    % SEM for each sequence (before grouping in pairs, so 32 sequences)
    sdERPs = std(allERPs,[],3,'omitnan');
    semERPs = sdERPs ./ sqrt(nERPs);
    
    % propagate the SEM for the pairs of sequences
    % according to sem_{(n_A*A + n_B*B)/(n_A+n_B)^2} = sqrt(n_A^2/(n_A+n_B)^2 sem_A^2 + n_B^2/(n_A+n_B)^2 sem_B^2)
    semERPs = sqrt(((nERPs.^2 .* semERPs.^2) + fliplr(nERPs.^2 .*semERPs.^2))./((nERPs + fliplr(nERPs)).^2));
    semERPs(:,n_seq/2 + 1:end) = [];
    
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
    semERPs = semERPs(:,seq_eff_order(n_back));
    
    figure;
    plot(meanERPs)
    
    % get the maxima and minima for all 16 sequences
    [max_erp, ind_max_erp] = max(meanERPs);
    [min_erp, ind_min_erp] = min(meanERPs);
    
    % get the maxima of the diff of the PHOT to mark stimulus onset
    [~, stim_onset] = max(diff(meanPHOTs));

    %standard errors of the mean for the maxima (use of linear indexing here)
    semMax = semERPs(sub2ind(size(semERPs),ind_max_erp,1:16));
    semMin = semERPs(sub2ind(size(semERPs),ind_min_erp,1:16));
    
    % plot ERPs for each sequence separately in a 4x4 plot
    % for each sequence, highlight where the maxima (red) and minima (blue) are located 
%     if aux_plots
        figure;
        y_limit = [min(meanERPs(:)) max(meanERPs(:))];

        load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');

        %this is just to help turn horizontal sequences into vertical ones
        ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]);

        for i = 1:16
           subplot(4,4,i);
           plot(meanERPs(:,i));
           ylim(y_limit);
           hold on;
           scatter(ind_max_erp(i), 0,40,'r','filled');
           scatter(ind_min_erp(i), 0,40,'b','filled');
           title(binomial_x_labels_latex{i}(ind_horiz));
        end
%     end
        
    % amplitude SEs and propagated SEM
    amplitudeSEs = max_erp - min_erp;
    semAmplSEs = sqrt(semMax.^2 + semMin.^2);
    
    % latency from stimulus onset to peak (not feasible to calculate error here)
    latencySEs = zeros(1,16);
    for i = 1:16
        latencySEs(i) = (ind_max_erp(i) - stim_onset(i))/blocks(1).resampleFreq;%this assumes for now that resampleFreq is always the same
    end
        
    % positive amplitude SEs and propagated SEM
    positiveAmplitudeSEs =  max_erp - meanERPs(1,:);
    semPosAmplSEs = sqrt(semMax.^2 + semERPs(1,:).^2);
    
    % negative amplitude and propagated SEM
    negativeAmplitudeSEs =  min_erp - meanERPs(1,:);
    semNegAmplSEs = sqrt(semMin.^2 + semERPs(1,:).^2);
    
    % put all results into a neat structure
    R = struct;
    
    R.amplitudeSEs = amplitudeSEs;
    R.semAmplSEs = semAmplSEs;
    
    R.latencySEs = latencySEs;
    
    R.positiveAmplitudeSEs = positiveAmplitudeSEs;
    R.semPosAmplSEs = semPosAmplSEs;
    
    R.negativeAmplitudeSEs = negativeAmplitudeSEs;
    R.semNegAmplSEs = semNegAmplSEs;
    
    R.meanERPs = meanERPs;
    R.meanPHOTs = meanPHOTs;
    R.nERPs = nERPs;
    
end