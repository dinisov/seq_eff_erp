function blocks = sortSEs(blocks, n_back)
%this function sorts ERPs according to the past sequence of events 

    n_seq = 2^n_back;

    for b = 1:length(blocks)
        
        LOCS = blocks(b).LOCS;
        LFP = blocks(b).LFP;
        PHOT = blocks(b).PHOT;
%         rawPHOT = blocks(b).rawPHOT;
        randomSequence = blocks(b).randomSequence;
        resampleFreq = blocks(b).resampleFreq;
        badTrials = blocks(b).badTrials;
        if isfield(blocks,'focusPeaks')
            focusPeaks = blocks.focusPeaks;
        end
        
        window = floor([-blocks(b).time_before_peak*resampleFreq, blocks(b).time_after_peak*resampleFreq]);

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
        
        %% PCA
%         [coeff,score] = pca(all_erps.','centered',false);
%         
%         figure;
%         for i=1:64
%            subplot(8,8,i);
%            plot(coeff(:,i));
%         end
        
        %%
        
        %calculate mean and SEM for outlier calculations
        meanERP = mean(all_erps, 2);
        STDs = std(all_erps, [], 2);

        n_sd = 3;

        % remove ERPs beyond n_sd (broadcasting here)
        outliers = all_erps < (meanERP - n_sd*STDs) | all_erps > (meanERP + n_sd*STDs);

        disp(['Data lost due to outliers: ' num2str(nnz(sum(outliers))/length(outliers)*100) '%']);

        good_erps = ~logical(sum(outliers));

        % remove ERP outliers
        ERPS = ERPS(:,:,good_erps);
        seqPHOT = seqPHOT(:,:,good_erps);
        badTrials = badTrials(good_erps);
        
        if isfield(blocks,'focusPeaks')
            focusPeaks = focusPeaks(good_erps);
            blocks(b).focusPeaks = focusPeaks;
        end
        
        blocks(b).ERPS = ERPS;
        blocks(b).badTrials = badTrials;
        blocks(b).seqPHOT = seqPHOT;
        blocks(b).window = window;
    
    end
    
end