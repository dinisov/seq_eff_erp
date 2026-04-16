function R = analyseSequentialEffects(blocks, aux_plots, plotSelector, reOrder, n_back, options)

    %% sort data according to previous sequence

    arguments
        blocks struct
        aux_plots double
        plotSelector double
        reOrder double
        n_back double
        options.suppressANOVA double = 0
        options.plotIndividualFlies double = 1
        options.scramLevel double = 0
    end

    %Extra
    suppressANOVA = options.suppressANOVA;
    plotIndividualFlies = options.plotIndividualFlies;
    scramLevel = options.scramLevel;

    %n_back = 5;  %TODO: make n_back variable
    
    % sort according to sequence
    blocks = sortSEs(blocks, n_back);

    %QA
    if size(blocks,2) > 1
        ['## Alert: blocks structure larger than expected ##']
        crash = yes
        %If this occurs, perhaps more than 1 block for this fly, or similar
            %This is bad because multiple following operations hardcoded use first element of blocks
    end
    
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

    %Scramble at sequence level if requested
    if scramLevel == 2
        scramInds = randperm( size(allERPs,2) );
        allERPs = allERPs(:,scramInds,:);
        allPHOTs = allPHOTs(:,scramInds,:);
        %NOTE: THIS DOES NOT DISTINGUISH BETWEEN ISOMERS
            %i.e. Scrambling will happily put LLLL as RRRRR etc 
            %It is probably possible to scramble *within* isomers, but more effort is needed
        disp(['Sequences scrambled (Scram 2)'])
    end

    
    % remove all ERPs with NaNs (from "shaving")
        %Disabled based on NaN-present nature
            %This can only be done safely on account of NaN pre-checks in sortSEs
    %%allERPs = allERPs(:,:,~isnan(sum(squeeze(sum(allERPs,2)))));
    
    %% onion plots
    
    % onion plot showing all ERPs for this block
    %all_erps = squeeze(sum(allERPs,2));
    all_erps = squeeze(nansum(allERPs,2)); %NaN version
    if aux_plots
        figure;
        plot(all_erps,'b');
        hold on;
        %plot(mean(all_erps,2), 'r');
        plot(nanmean(all_erps,2), 'r');
    end
    
    % onion plot of the photodiode data
    if aux_plots
        figure;
        plot(squeeze(nansum(allPHOTs,2)),'b');
        hold on
        %plot(mean(squeeze(sum(allPHOTs,2)),2),'r');
        plot(nanmean(squeeze(nansum(allPHOTs,2)),2),'r');
        title('')
    end
    
    %% plot mean ERPs for both stimuli separately
    
    %plotSeparateERPs(allERPs);
    if plotIndividualFlies
        plotSeparateERPs(allERPs);
    end
    
    %% plot isomers
    
    %AAAARRRRAARARRARa
%     plotIsomers(allERPs,[1 0 1 0 1 1 1 1 1 0 1 1 0 0 0 1 1 0 1 0], window, n_back, blocks(1).resampleFreq);
    %plotIsomers(allERPs,[], window, n_back, blocks(1).resampleFreq);
    %plotIsomers(allERPs,[], window, n_back, blocks(1).resampleFreq,...
    %    blocks(1).transectTime, blocks(1).avTransectWindow, plotSelector, reOrder,...
    %    [blocks(1).date,' (#',num2str(blocks(1).fly),') B',blocks(1).block]);
    [ISOMER] = plotIsomers(allERPs,[], window, n_back, blocks(1).resampleFreq,...
        blocks(1).transectTime, blocks(1).avTransectWindow, plotSelector, reOrder,...
        [blocks(1).date,' (#',num2str(blocks(1).fly),') B',blocks(1).block],...
        plotIndividualFlies, allPHOTs);
    
    %% calculate SEs

    allERPs( isnan(allERPs) ) = 0; %This nonsensical step required because flipping below requires 0s
        %That is, adding 0 to a value -> value, but adding NaN to a value -> NaN, so we briefly convert NaNs to 0s 
    
    % join ERPs corresponding to the same pattern (01001 and 10110 and so on)
    allERPs = allERPs + flip(allERPs,2);
    allPHOTs = allPHOTs + flip(allPHOTs,2);
    %allERPs = allERPs(:,1:16,:);
    %allPHOTs = allPHOTs(:,1:16,:);
    allERPs = allERPs(:,1:0.5*2^n_back,:);
    allPHOTs = allPHOTs(:,1:0.5*2^n_back,:);
    
    % reorder according to the literature
    allERPs = allERPs(:,seq_eff_order(n_back),:);
    allPHOTs = allPHOTs(:,seq_eff_order(n_back),:);

    %MOVED HERE FROM calculateSEs TO SIMPLIFY ZERO-SHIFTING MATTERS
    allERPs(allERPs == 0) = nan;

    %Shift all above zero, if requested (And data actually dips below zero)
        %Pushed into individual functions
    %if zeroShiftMode == 1 && nanmin(allERPs,[],'all') < 0
    %    allERPs = allERPs - nanmin(allERPs,[],'all');    
    %    disp(['Zero shift applied'])
    %end
    
    %R = calculateSEs(allERPs,allPHOTs,aux_plots,window,blocks(1).resampleFreq);
    R = calculateSEs(allERPs,allPHOTs,aux_plots,window,blocks(1).resampleFreq, blocks(1).transectTime, blocks(1).avTransectWindow, plotSelector, n_back);
    %Inject isomer data
    R.ISOMER = ISOMER;
    R = timeFrequencySpectrum(R, blocks, n_back);
    
    % add window to results structure for convenience
    R.window = blocks(1).window;

    R.avTransectWindow = blocks(1).avTransectWindow; %Okay to force first block?

    R.fly = blocks(1).fly;
    R.date = blocks(1).date;
    R.block = blocks(1).block;

    %R.SEQS = blocks(1).SEQS; %Disabled for the moment until flipped/compacted/etc (See above)

    
    %% ANOVA
    if ~suppressANOVA
        SEAnova(R, reOrder, n_back);
    end
    
end