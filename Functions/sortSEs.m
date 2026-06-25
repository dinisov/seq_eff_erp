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

    n_seq = 2^n_back;

    for b = 1:length(blocks)
        
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

        %ERPS = zeros(length(window(1):window(2)), n_seq, sequenceLength);
        %seqPHOT = zeros(length(window(1):window(2)), n_seq, sequenceLength);
        ERPS = nan(length(window(1):window(2)), n_seq, sequenceLength); %Switch to NaN because zero issues
        seqPHOT = nan(length(window(1):window(2)), n_seq, sequenceLength);
        seqTIME = nan(length(window(1):window(2)), n_seq, sequenceLength);
        SEQS = cell(1,n_seq);

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
        
        for n = n_back:sequenceLength

            % decimal value of binary sequence of length n_back
            %%seq = bin2dec(num2str(randomSequence(n-n_back+1:n))) + 1;
            if ~arrowMode
                seq = bin2dec(num2str(randomSequence(n-n_back+1:n))) + 1;
                SEQS{1,seq} = randomSequence(n-n_back+1:n); %Store this for posterity
            else
                seq = bin2dec(num2str(arrowSequence(n-n_back+1:n))) + 1;
                SEQS{1,seq} = arrowSequence(n-n_back+1:n);
            end

            % stack ERPs and PHOTs along third dimension (first two dims are sequence and
            % time respectively)
            %%ERPS(:, seq, n) = LFP(LOCS(n) + window(1) : LOCS(n) + window(2));
            %%seqPHOT(:, seq, n) = normalize(PHOT(2-randomSequence(n), LOCS(n) + window(1) : LOCS(n) + window(2) ));
            %%seqTIME(:, seq, n) = TIMES(LOCS(n) + window(1) : LOCS(n) + window(2));

            if ~arrowMode
                theseLOCInds = LOCS(n) + window(1) : LOCS(n) + window(2);
                ERPS(:, seq, n) = LFP( theseLOCInds );
                seqPHOT(:, seq, n) = normalize(PHOT(2-randomSequence(n), theseLOCInds ));
                seqTIME(:, seq, n) = TIMES( theseLOCInds );
            else
                theseLOCInds = arrowLOCS(n) + window(1) : arrowLOCS(n) + window(2);
                if arrowMode == 1
                    ERPS(:, seq, n) = arrowLFP( theseLOCInds );
                    seqPHOT(:, seq, n) = normalize(arrowPHOT(2-arrowSequence(n), theseLOCInds ));
                    %Note: This selects the applicable phot row for this event, 
                    % but if the window spans >1 event, the next event may be not shown, 
                    % because it occurred in the other channel
                    seqTIME(:, seq, n) = arrowTIMES( theseLOCInds );
                elseif arrowMode == 1.5
                    ERPS(:, seq, n) = fliplr( arrowLFP( theseLOCInds ) );
                    seqPHOT(:, seq, n) = fliplr( normalize(arrowPHOT(2-arrowSequence(n), theseLOCInds )) );
                    seqTIME(:, seq, n) = fliplr( arrowTIMES( theseLOCInds ) );
                end
            end

            %SEQS{1,seq} = randomSequence(n-n_back+1:n); %Store this for posterity

            %Testatory plots, if requested
            if firstLastPlot
                %h = figure;
                %c = 1;
                figure(h)
                if n == n_back || n == sequenceLength
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
                    scatter( LOCSActual( n-n_back+1:n ) , repmat( nanmax(LFP)*0.8 , 1 , n_back ) )
                    for i = n-n_back+1:n
                        text([LOCSActual(i)],[nanmax(LFP)*0.78],num2str(i))
                    end

                    %Plot ERP
                    %plot( sort(theseLOCIndsActual), ERPS(:, seq, n) ) %Borrow from above
                        %Note: The sort applied to the inds is intended to force unflipped arrow data (If applicable) to appear unflipped
                            %(Since plotting with original flipped indices will conceal this)
                     plot(theseLOCIndsActual, ERPS(:, seq, n) ) %Disabled sort
                        %Without sort, the ERP should be read from green to red to understand its shape
                        %Secondary note: Since no subsampling, may look different to 'original' base LFP data, even if same orientation

                    xlim([ nanmin(theseLOCIndsActual)-n_back*range(theseLOCIndsActual)*1.1 ,...
                        nanmax(theseLOCIndsActual)+n_back*range(theseLOCIndsActual)*1.1 ])

                    titleStr = ['n ',num2str(n-n_back+1),' : ',num2str(n), ' capture'];
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
        blocks(b).SEQS = SEQS;
    
    end
    
end