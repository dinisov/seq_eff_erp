function blocks = sortSEs(blocks, n_back, options)
%this function sorts ERPs according to the past sequence of events 
arguments
    blocks struct
    n_back double
    options.arrowMode double = 0 %Whether to implement 'Arrow of Time' control (Collection of data from immediately before stimuli)
        %arrowMode = 1 - randomSequence read backwards, data collected from 'behind' and is backwards (e.g. n1:5 is read as 5:1, and window is 0:1 etc)
        %arrowMode = 1.5 - As above, but LFP/PHOT/Times are reflipped to return them to their 'true' orientation
    options.firstLastPlot double = 0 %Whether to do a testatory first/last plot
end

    arrowMode = options.arrowMode;
    if arrowMode
        disp(['-> Arrow mode requested <-'])
    end
    firstLastPlot = options.firstLastPlot;

    %n_seq = 2^nBackActual;

    for b = 1:length(blocks)

        transProbDesign = blocks(b).transProbDesign;
        if ~transProbDesign
            n_seq = 2^n_back; %Moved slightly inefficiently here so as to be part of per-block type checking
            nBackActual = n_back; %Use active form for multiple type support
        else
            disp(['Block to be sorted for transition probability nature'])
            n_stimuli = blocks(b).transProbAncillary.nStimuli;
            nBackActual = 2; %Hardcoded, for now
            disp([num2str(n_stimuli),' stimuli x ',num2str(nBackActual),' transition nBack = ',num2str(n_stimuli^nBackActual),' "sequences"'])
            %Old, calculated-here system
            %{
            transMatrix = []; %TO DO: SWITCH TO LOADING THIS FROM A CENTRAL FILE?
            for r = 1:n_stimuli
                for c = 1:n_stimuli
                    transMatrix{r,c} = [r,c];
                end
            end
            transMatrix = cell2mat( reshape(transMatrix, n_stimuli^nBackActual, 1) );
                %Thus, canonically, the sequences will be ordered as 1->1, 2->1, 3->1, 4->1, 1->2, etc
            %}
            %New, load-from-central system
            try
                temp = load(['transMatrix_',num2str(n_stimuli),'stim_',num2str(nBackActual),'back.mat']);
                transMatrix = temp.transMatrix;
                transLabels = temp.transMatrixLegend; %Save legend for later use with calculateSEs; Note that this is the full form, like transMatrix
                histLabels = temp.transMatrixHistLabels; %Save isomer-friendly version as well
            catch
                ['-# Could not load transition probability legend/matrix #-']
                crash = yes
            end
        end
        
        LOCS = blocks(b).LOCS;
        LFP = blocks(b).LFP;
        PHOT = blocks(b).PHOT;
        TIMES = blocks(b).times;
%         rawPHOT = blocks(b).rawPHOT;
        randomSequence = blocks(b).randomSequence;
        resampleFreq = blocks(b).resampleFreq;
        badTrials = blocks(b).badTrials;
        
        if isfield(blocks,'focusPeaks')
            focusPeaks = blocks.focusPeaks;
        end
        
        window = floor(resampleFreq*blocks(b).window);

        sequenceLength = length(randomSequence);
        %Note: No explicit check that this will match up with transProb nBack?

        if ~transProbDesign 
            %ERPS = zeros(length(window(1):window(2)), n_seq, sequenceLength);
            %seqPHOT = zeros(length(window(1):window(2)), n_seq, sequenceLength);
            ERPS = nan(length(window(1):window(2)), n_seq, sequenceLength); %Switch to NaN because zero issues
            seqPHOT = nan(length(window(1):window(2)), n_seq, sequenceLength);
            seqTIME = nan(length(window(1):window(2)), n_seq, sequenceLength);
            SEQS = cell(1,n_seq);
        else
            ERPS = nan(length(window(1):window(2)), n_stimuli^nBackActual, sequenceLength); %Switch to NaN because zero issues
            seqPHOT = nan(length(window(1):window(2)), n_stimuli^nBackActual, sequenceLength);
            seqTIME = nan(length(window(1):window(2)), n_stimuli^nBackActual, sequenceLength);
            STIMS = cell(1,n_stimuli^nBackActual);
        end

        %Pre-check for 'true' NaNs
        if any(isnan(LFP))
            disp(['-# Caution: True NaN/s present in LFP data #-'])
            %crash = yes %Maybe overkill; Will have to see if this ever occurs
        end
        %Prepare arrow data if applicable
        if arrowMode
            arrowSequence = fliplr(randomSequence); %CHECK IF ACTUALLY NECESSARY
            arrowLOCS = fliplr( numel(LFP) - LOCS); %e.g. If LOCS are [7,8,9] for a LFP len 10, this converts them to [1,2,3]
            arrowLFP = fliplr(LFP);
            arrowPHOT = fliplr(PHOT);
            arrowTIMES = fliplr(TIMES);
        end

        %Prepare some plot stuff
        if firstLastPlot
            h = figure;
            set(gcf,'Name', ['Fly ',num2str(blocks(b).fly), ' Block ',blocks(b).block,' First/Last plot'])
            c = 1;
        end
        
        for n = nBackActual:sequenceLength

            %Pre-check to make sure not attempting to acquire LFP data from after end  of experiment (i.e. If stimuli still occurring at end)
            if ~arrowMode && ( ( LOCS(n) + window(2) > size(LFP,2) ) || ( LOCS(n) + window(2) > size(PHOT,2) )  ) %Borrow below indicising
                disp(['-# Attempted acquisition of event #',num2str(n),' would exceed LFP data; Skipping #-'])
                continue %NOTE: MAY HAVE ISSUES BY LEAVING NaNs IN ERPS?
                %Note: Choosing deliberately to not modify randomSequence, but if it is to ever be returned, this will be a potential issue
            elseif arrowMode && ( (arrowLOCS(n) + window(2) > size(LFP,2)) || (arrowLOCS(n) + window(2) > size(PHOT,2)) ) %Ditto
                disp(['-# Attempted acquisition of event #',num2str(n),' would exceed LFP data; Skipping #-'])
                continue %NOTE: MAY HAVE ISSUES BY LEAVING NaNs IN ERPS?
            end

            % decimal value of binary sequence of length nBackActual
            %%seq = bin2dec(num2str(randomSequence(n-nBackActual+1:n))) + 1;
            if ~arrowMode
                if ~transProbDesign
                    seq = bin2dec(num2str(randomSequence(n-nBackActual+1:n))) + 1; %This should evaluate to a number between 1 and 32 eg
                    SEQS{1,seq} = randomSequence(n-nBackActual+1:n); %Store this for posterity
                else
                    [yen,seq] = ismember( randomSequence(n-nBackActual+1:n) , transMatrix, 'rows' );
                    %QA
                    if yen == 0
                        ['## Transition sequence not found in matrix! ##']
                        crash = yes
                    end
                    STIMS{1,seq} = randomSequence(n-nBackActual+1:n);
                end
            else
                if ~transProbDesign
                    seq = bin2dec(num2str(arrowSequence(n-nBackActual+1:n))) + 1;
                    SEQS{1,seq} = arrowSequence(n-nBackActual+1:n);
                else
                    [yen,seq] = ismember( arrowSequence(n-nBackActual+1:n) , transMatrix, 'rows' );
                    %QA
                    if yen == 0
                        ['## (Arrowed) Transition sequence not found in matrix! ##']
                        crash = yes
                    end
                    STIMS{1,seq} = arrowSequence(n-nBackActual+1:n);
                end
            end

            % stack ERPs and PHOTs along third dimension (first two dims are sequence and
            % time respectively)
            %%ERPS(:, seq, n) = LFP(LOCS(n) + window(1) : LOCS(n) + window(2));
            %%seqPHOT(:, seq, n) = normalize(PHOT(2-randomSequence(n), LOCS(n) + window(1) : LOCS(n) + window(2) ));
            %%seqTIME(:, seq, n) = TIMES(LOCS(n) + window(1) : LOCS(n) + window(2));

            if ~arrowMode
                theseLOCInds = LOCS(n) + window(1) : LOCS(n) + window(2);
                ERPS(:, seq, n) = LFP( theseLOCInds );
                if ~transProbDesign
                    seqPHOT(:, seq, n) = normalize(PHOT(2-randomSequence(n), theseLOCInds )); %seqPHOT is a rather sneaky use of randomSequence to pull applicable PHOT channel
                else
                    seqPHOT(:, seq, n) = normalize(nansum( PHOT(:, theseLOCInds ),1) ); %For transition probabilities I currently cannot be bothered to preserve the individual phot channel information
                end
                seqTIME(:, seq, n) = TIMES( theseLOCInds );
            else
                theseLOCInds = arrowLOCS(n) + window(1) : arrowLOCS(n) + window(2);
                if arrowMode == 1
                    ERPS(:, seq, n) = arrowLFP( theseLOCInds );
                    if ~transProbDesign
                        seqPHOT(:, seq, n) = normalize(arrowPHOT(2-arrowSequence(n), theseLOCInds ));
                        %Note: This selects the applicable phot row for this event, 
                        % but if the window spans >1 event, the next event may be not shown, 
                        % because it occurred in the other channel
                    else
                        seqPHOT(:, seq, n) = normalize(nansum( PHOT(:, theseLOCInds ),1) ); %This may misrepresent the data?
                    end
                    seqTIME(:, seq, n) = arrowTIMES( theseLOCInds );
                elseif arrowMode == 1.5
                    ERPS(:, seq, n) = fliplr( arrowLFP( theseLOCInds ) );
                    if ~transProbDesign
                        seqPHOT(:, seq, n) = fliplr( normalize(arrowPHOT(2-arrowSequence(n), theseLOCInds )) );
                    else
                        seqPHOT(:, seq, n) = normalize(nansum( PHOT(:, theseLOCInds ),1) ); %Again: Not tested for arrowMode
                    end
                    seqTIME(:, seq, n) = fliplr( arrowTIMES( theseLOCInds ) );
                end
            end

            %SEQS{1,seq} = randomSequence(n-nBackActual+1:n); %Store this for posterity

            %Testatory plots, if requested
            if firstLastPlot
                %h = figure;
                %c = 1;
                figure(h)
                if n == nBackActual || n == sequenceLength
                    subplot(2,1,c)
                    ploti = 1:resampleFreq/100:numel(LFP); %Subsample so plot isn't as supermassive
                    %Plot base LFP data
                    plot( ploti , LFP(ploti) ) %Only plot OG data?
                    hold on
                    photi = [1:size(PHOT,2)]; %No subsampling?
                    for i = 1:size(PHOT,1)
                        plot( photi, normalize( PHOT(i,photi), 'range', [nanmin(LFP)+0.1*range(LFP), nanmin(LFP)+0.2*range(LFP)] ) )
                    end
                    %plot()
                    %Plot associated LOCs
                    if arrowMode
                        theseLOCIndsActual = numel(LFP) - theseLOCInds; %Note: No fliplr applied here
                        LOCSActual = numel(LFP) - arrowLOCS;
                    else
                        theseLOCIndsActual = theseLOCInds;
                        LOCSActual = LOCS;
                    end
                    %Window extent
                    line([theseLOCIndsActual(1),theseLOCIndsActual(end)],...
                        [LFP(theseLOCIndsActual(1)),LFP(theseLOCIndsActual(end))],'LineStyle','--', 'Color','k')
                    %Window start
                    line([theseLOCIndsActual(1),theseLOCIndsActual(1)],[nanmin(LFP),nanmax(LFP)],'Color','g')
                    %Window end
                    line([theseLOCIndsActual(end),theseLOCIndsActual(end)],[nanmin(LFP),nanmax(LFP)],'Color','r')

                    scatter( LOCSActual, LFP(LOCSActual) ) %All LOCS
                    scatter( LOCSActual( n-nBackActual+1:n ) , repmat( nanmax(LFP)*0.8 , 1 , nBackActual ) )
                    for i = n-nBackActual+1:n
                        text([LOCSActual(i)],[nanmax(LFP)*0.78],num2str(i))
                    end

                    %Plot ERP
                    %plot( sort(theseLOCIndsActual), ERPS(:, seq, n) ) %Borrow from above
                        %Note: The sort applied to the inds is intended to force unflipped arrow data (If applicable) to appear unflipped
                            %(Since plotting with original flipped indices will conceal this)
                     plot(theseLOCIndsActual, ERPS(:, seq, n) ) %Disabled sort
                        %Without sort, the ERP should be read from green to red to understand its shape
                        %Secondary note: Since no subsampling, may look different to 'original' base LFP data, even if same orientation

                    xlim([ nanmin(theseLOCIndsActual)-nBackActual*range(theseLOCIndsActual)*1.1 ,...
                        nanmax(theseLOCIndsActual)+nBackActual*range(theseLOCIndsActual)*1.1 ])

                    titleStr = ['n ',num2str(n-nBackActual+1),' : ',num2str(n), ' capture'];
                    if arrowMode 
                        titleStr =  [titleStr,' [arrowMode ',num2str(arrowMode),']'];
                    end
                    title(titleStr)

                    c = c + 1;
                end
            end

        end

        % matrix with all ERPs irrespective of sequence
        %all_erps = squeeze(sum(ERPS,2));% squeeze removes first dimension
        all_erps = squeeze(nansum(ERPS,2));% squeeze removes first dimension; Need to use nansum now
        
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
        %QA
        if nansum(isnan(meanERP)) ~= 0
            ['## ALERT: NAN CONTAMINATION IN meanERP ##']
            crash = yes         
        end

        n_sd = 4;

        % remove ERPs beyond n_sd (broadcasting here)
        outliers = all_erps < (meanERP - n_sd*STDs) | all_erps > (meanERP + n_sd*STDs);

        disp(['Data lost due to outliers: ' num2str(nnz(sum(outliers))/length(outliers)*100) '%']);

        good_erps = ~logical(sum(outliers));

        % remove ERP outliers
        ERPS = ERPS(:,:,good_erps);
        seqPHOT = seqPHOT(:,:,good_erps);
        seqTIME = seqTIME(:,:,good_erps);
        badTrials = badTrials(good_erps);
        
        if isfield(blocks,'focusPeaks')
            focusPeaks = focusPeaks(good_erps);
            blocks(b).focusPeaks = focusPeaks;
        end
        
        blocks(b).ERPS = ERPS;
        blocks(b).badTrials = badTrials;
        blocks(b).seqPHOT = seqPHOT;
        blocks(b).seqTIME = seqTIME;
        if ~transProbDesign
            blocks(b).SEQS = SEQS;
        else
            blocks(b).STIM = STIMS;
            blocks(b).transMatrix = transMatrix;
            blocks(b).transProbAncillary.nBackActual = nBackActual;
            blocks(b).transProbAncillary.transLabels = transLabels;
            blocks(b).transProbAncillary.histLabels = histLabels;
        end
    
    end
    
end