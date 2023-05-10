function R = analyseSequentialEffects(blocks, aux_plots)

    %% sort data according to previous sequence

    n_back = 5;  %TODO: make n_back variable
    
    blocks = sortSEs(blocks, n_back, aux_plots);
    
    window = blocks(1).window;%WARNING: we can only merge blocks with windows of same size
    
    %% group blocks

    [allERPs, allPHOTs, goodTrials] = groupBlocks(blocks,window,n_back);

    %% get rid of bad trials (trials with too long gaps between peaks)
    allERPs = allERPs(:,:,logical(goodTrials));
    allPHOTs = allPHOTs(:,:,logical(goodTrials));
    
    %% remove all ERPs with NaNs (from "shaving")
    allERPs = allERPs(:,:,~isnan(sum(squeeze(sum(allERPs,2)))));
    
    %% onion plots
    
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
        title('')
    end
    
    %% plot mean ERPs for both stimuli separately
    
    plotSeparateERPs(allERPs);
    
    %AAAARRRRAARARRARa
%     plotIsomers(allERPs,[1 0 1 0 1 1 1 1 1 0 1 1 0 0 0 1 1 0 1 0]);
    
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
    
    % SEM for each sequence (this outputs a t X seq matrix of SDs)
    sdERPs = std(allERPs,[],3,'omitnan');
    semERPs = sdERPs ./ sqrt(nERPs);% broadcasting
    
    % reorder according to the literature
    meanERPs = meanERPs(:,seq_eff_order(n_back));
    meanPHOTs = meanPHOTs(:,seq_eff_order(n_back));
    semERPs = semERPs(:,seq_eff_order(n_back));
    nERPs = nERPs(seq_eff_order(n_back));
    
    if aux_plots
        figure;
        plot(meanERPs)
    end
    
    % get the maxima and minima for all 16 sequences
    [max_erp, ind_max_erp] = max(meanERPs);
    [min_erp, ind_min_erp] = min(meanERPs);
    
    % get the maxima of the diff of the PHOT to mark stimulus onset
%     [~, stim_onset] = max(diff(meanPHOTs));

    %SEM for the maxima and minima(use of linear indexing here)
    semMax = semERPs(sub2ind(size(semERPs),ind_max_erp,1:16));
    semMin = semERPs(sub2ind(size(semERPs),ind_min_erp,1:16));
    
    %% ANOVA
    allERPs = allERPs(:,seq_eff_order(n_back),:);
    
    % full data for maxima and minima across all stimuli
    aux = squeeze(reshape(allERPs,[1,size(allERPs,1)*size(allERPs,2),size(allERPs,3)])).';
    
    % data tables for performing an ANOVA to check for the effect of
    % sequence
    dataMaxima = aux(:,sub2ind(size(semERPs),ind_max_erp,1:16));
    dataMinima = aux(:,sub2ind(size(semERPs),ind_min_erp,1:16));
    dataAmplitude = dataMaxima-dataMinima;
    
    anova1(dataAmplitude);
    
    %%
    
    % plot ERPs for each sequence separately in a 4x4 plot
    % for each sequence, highlight where the maxima (red) and minima (blue) are located 
    if aux_plots
        
        figure;

        load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');

        %this is just to help turn horizontal sequences into vertical ones
        ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]);

        for i = 1:16
           subplot(4,4,i);
           plot(normalize(meanERPs(:,i)));
           hold on;
           scatter(ind_max_erp(i), 0,40,'r','filled');
           scatter(ind_min_erp(i), 0,40,'b','filled');
           plot(normalize(meanPHOTs(:,i)));
%            plot([window(1) window(1)],ylim,'r');
           title(binomial_x_labels_latex{i}(ind_horiz));
        end
    end
        
    % amplitude SEs and propagated SEM
    amplitudeSEs = max_erp - min_erp;
    semAmplSEs = sqrt(semMax.^2 + semMin.^2);
    
    % latency from stimulus onset to peak (not feasible to calculate error here)
    % this assumes for now that resampleFreq is the same for all blocks
    % the entire window is defined by the stimulus onset so all we need to
    % do is use the left side of the window to correct the time to
    % peak/trough
    latencyToPeakSEs = zeros(1,16);
    latencyToTroughSEs = zeros(1,16);
    for i = 1:16
        latencyToPeakSEs(i) = (ind_max_erp(i)-window(1))/blocks(1).resampleFreq;
        latencyToTroughSEs(i) = (ind_min_erp(i)-window(1))/blocks(1).resampleFreq;
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
    
    R.latencyToPeakSEs = latencyToPeakSEs;
    R.latencyToTroughSEs = latencyToTroughSEs;
    
    R.positiveAmplitudeSEs = positiveAmplitudeSEs;
    R.semPosAmplSEs = semPosAmplSEs;
    
    R.negativeAmplitudeSEs = negativeAmplitudeSEs;
    R.semNegAmplSEs = semNegAmplSEs;
    
    R.meanERPs = meanERPs;
    R.meanPHOTs = meanPHOTs;
    R.nERPs = nERPs;
    
    R.allERPs = allERPs;
    
end