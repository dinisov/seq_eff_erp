function R = analyseSequentialEffects(blocks, aux_plots)

    %% sort data according to previous sequence

    n_back = 5;  %TODO: make n_back variable
    
    % sort according to sequence
    blocks = sortSEs(blocks, n_back);
    
    window = floor(blocks(1).window * blocks(1).resampleFreq);
    
    % group blocks (WARNING: only possible if window is the same for all blocks)
    [allERPs, allPHOTs, goodTrials, focusPeaks] = groupBlocks(blocks,n_back);

    %TODO: find way to handle errors in block experiments
    if isfield(blocks,'focusPeaks')
        % get rid of all trials except the last one in each train
        allERPs = allERPs(:,:,logical(focusPeaks));
        allPHOTs = allPHOTs(:,:,logical(focusPeaks));
    else
        % get rid of bad trials (trials with too long gaps between peaks)
        allERPs = allERPs(:,:,logical(goodTrials));
        allPHOTs = allPHOTs(:,:,logical(goodTrials));
    end
    
    % remove all ERPs with NaNs (from "shaving")
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
    
    %% plot isomers
    
    %AAAARRRRAARARRARa
%     plotIsomers(allERPs,[1 0 1 0 1 1 1 1 1 0 1 1 0 0 0 1 1 0 1 0], window, n_back, blocks(1).resampleFreq);
    plotIsomers(allERPs,[], window, n_back, blocks(1).resampleFreq);
    
    %% calculate SEs
    
    % join ERPs corresponding to the same pattern (01001 and 10110 and so on)
    allERPs = allERPs + flip(allERPs,2);
    allPHOTs = allPHOTs + flip(allPHOTs,2);
    allERPs = allERPs(:,1:16,:);
    allPHOTs = allPHOTs(:,1:16,:);
    
    % reorder according to the literature
    allERPs = allERPs(:,seq_eff_order(n_back),:);
    allPHOTs = allPHOTs(:,seq_eff_order(n_back),:);
    
    R = calculateSEs(allERPs,allPHOTs,aux_plots,window,blocks(1).resampleFreq);
    
    R = timeFrequencySpectrum(R, blocks);
    
    % add window to results structure for convenience
    R.window = blocks(1).window;
    
    %% ANOVA
    SEAnova(R);
    
end