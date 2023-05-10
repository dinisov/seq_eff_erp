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
%     plotIsomers(allERPs,[1 0 1 0 1 1 1 1 1 0 1 1 0 0 0 1 1 0 1 0],n_back);
    
    % join ERPs corresponding to the same pattern (01001 and 10110 and so on)
    allERPs = allERPs + flip(allERPs,2);
    allPHOTs = allPHOTs + flip(allPHOTs,2);
    allERPs = allERPs(:,1:16,:);
    allPHOTs = allPHOTs(:,1:16,:);
    
    %mean photodiode traces
    allPHOTs(allPHOTs == 0) = nan;
    meanPHOTs = mean(allPHOTs, 3, 'omitnan');
    meanPHOTs = meanPHOTs(:,seq_eff_order(n_back));
    
    R = calculateSEs(allERPs,aux_plots,window,n_back,blocks(1).resampleFreq);
    
    R.meanPHOTs = meanPHOTs;
    
    % plot ERPs for each sequence separately in a 4x4 plot
    % for each sequence, highlight where the maxima (red) and minima (blue) are located 
%     if aux_plots
%         
%         figure;
% 
%         load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');
% 
%         %this is just to help turn horizontal sequences into vertical ones
%         ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]);
% 
%         for i = 1:16
%            subplot(4,4,i);
%            plot(normalize(meanERPs(:,i)));
%            hold on;
%            scatter(ind_max_erp(i), 0,40,'r','filled');
%            scatter(ind_min_erp(i), 0,40,'b','filled');
%            plot(normalize(meanPHOTs(:,i)));
% %            plot([window(1) window(1)],ylim,'r');
%            title(binomial_x_labels_latex{i}(ind_horiz));
%         end
%     end
    
end