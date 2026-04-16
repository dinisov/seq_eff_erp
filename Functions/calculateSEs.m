function R = calculateSEs(allERPs,allPHOTs,aux_plots,window, resampleFreq, transectTime, avTransectWindow, plotSelector, n_back)
%calculateSEs Takes a matrix of ERPs separated by sequence and calculates SEs
%   Detailed explanation goes here

    nActive = 0.5*2^n_back; %Replaces hardcoded 16 in multiple places

    R = struct;
    
    if ~isempty(allPHOTs)
        %mean photodiode traces
        allPHOTs(allPHOTs == 0) = nan;
        meanPHOTs = mean(allPHOTs, 3, 'omitnan');
        R.meanPHOTs = meanPHOTs;
    end
    
    % in preparation for calculating the nan mean
    %allERPs(allERPs == 0) = nan;
    meanERPs = mean(allERPs, 3, 'omitnan');

    %number of ERPs for each of 16 sequence (equal to number of PHOTs)
    nERPs = sum(~isnan(allERPs(1,:,:)), 3);
    
    % SEM for each sequence (this outputs a t X seq matrix of SDs)
    sdERPs = std(allERPs,[],3,'omitnan');
    semERPs = sdERPs ./ sqrt(nERPs);% broadcasting
    
    if aux_plots
        figure;
        plot(meanERPs);
    end
    
    % get the maxima and minima for all 16 sequences
    [max_erp, ind_max_erp] = max(meanERPs);
    [min_erp, ind_min_erp] = min(meanERPs);
    
    % get the maxima of the diff of the PHOT to mark stimulus onset
%     [~, stim_onset] = max(diff(meanPHOTs));

    %SEM for the maxima and minima(use of linear indexing here)
    %semMax = semERPs(sub2ind(size(semERPs),ind_max_erp,1:16));
    %semMin = semERPs(sub2ind(size(semERPs),ind_min_erp,1:16));
    semMax = semERPs(sub2ind(size(semERPs),ind_max_erp,1:nActive));
    semMin = semERPs(sub2ind(size(semERPs),ind_min_erp,1:nActive));
    
    %%
        
    % amplitude SEs and propagated SEM
    amplitudeSEs = max_erp - min_erp;
    semAmplSEs = sqrt(semMax.^2 + semMin.^2);
    
    % latency from stimulus onset to peak (not feasible to calculate error here)
    % this assumes for now that resampleFreq is the same for all blocks
    % the entire window is defined by the stimulus onset so all we need to
    % do is use the left side of the window to correct the time to
    % peak/trough
    latencyToPeakSEs = zeros(1,nActive);
    latencyToTroughSEs = zeros(1,nActive);
    for i = 1:nActive
        latencyToPeakSEs(i) = (ind_max_erp(i)-window(1))/resampleFreq;
        latencyToTroughSEs(i) = (ind_min_erp(i)-window(1))/resampleFreq;
    end
        
    % positive amplitude SEs and propagated SEM
    positiveAmplitudeSEs =  max_erp - meanERPs(1,:); %Note: If max amp occurs at comparison index (1 here) then this will yield a [true] 0
    semPosAmplSEs = sqrt(semMax.^2 + semERPs(1,:).^2);
    
    % negative amplitude and propagated SEM
    negativeAmplitudeSEs =  min_erp - meanERPs(1,:); %Ditto, for min amp (probably)
    semNegAmplSEs = sqrt(semMin.^2 + semERPs(1,:).^2);

    %Transect (New)
    %disp(['Mean ERP at transect (',num2str(transectTime),' of ',num2str(meanERPs,2),' timepoints)'])
    %disp(['Mean ERP at transect t ',num2str(transectTime),' of ',num2str(size(meanERPs,1)),' total timepoints'])
    %size(allERPs)
    %size(meanERPs)
    transectSEs = abs(meanERPs(transectTime,:)); %Collect mean ERP at transectTime
    %size(semERPs)
    %semTransectSEs = semERPs(sub2ind(size(semERPs),transectTime,1:16)); %Collect singular SEM from same time
    semTransectSEs = semERPs(transectTime,1:nActive); %No need sub2ind when 1 timepoint
        %Note: Only 80% confident in correctness of this SEM line
    %Transect across window
    if isempty(avTransectWindow)
        %avTransectWindowSEs = [];
        %semAvTransectWindowSEs = [];
        if plotSelector(7) %Don't bother mentioning this if not actually looking at transect window
            disp(['No transect window specified; Using provided transect time (',num2str(transectTime),') +- 10'])
        end
        avTransectWindow = [transectTime-10,transectTime+10];
        %Quick QA
        if avTransectWindow(1) < 1
            avTransectWindow(1) = 1;
            if plotSelector(7)
                ['-# Caution: Transect first coordinate calculated to be smaller than 0; Rectifying #-']
            end
        end
        if avTransectWindow(2) > size(meanERPs,1)
            avTransectWindow(2) = size(meanERPs,1);
            if plotSelector(7)
                ['-# Caution: Transect second coordinate calculated to be larger than time; Rectifying #-']
            end
        end
    end
    %if ~isempty(avTransectWindow)
    avTransectWindowSEs = abs( nanmean( meanERPs( avTransectWindow(1):avTransectWindow(2) ,:),1) ); %Collect mean ERP across transect window
    %semAvTransectWindowSEs = nanstd( allERPs( avTransectWindow(1):avTransectWindow(2), :, : ), [], [1,3] ) ./ ...
    %    sqrt( nansum( ~isnan( allERPs( avTransectWindow(1):avTransectWindow(2), :, : ) ), [1,3] ) );
        %Calculate SEM manually from underlying, divide by number of events per sequence
    semAvTransectWindowSEs = mean( nanstd( allERPs( avTransectWindow(1):avTransectWindow(2), :, : ), [], [3] ) ./ ...
        sqrt( nansum( ~isnan( allERPs( avTransectWindow(1):avTransectWindow(2), :, : ) ), [3] ) ) , 1 );
        %Calculate SEM per each slice of transect, then average
    %end

    %testatory qa
        %Note: This QA defunct, given current assumption that zeroes are in fact 'true'
    %if any(amplitudeSEs == 0) | any(positiveAmplitudeSEs == 0) | any(negativeAmplitudeSEs == 0) | any(transectSEs == 0) | any(avTransectWindowSEs == 0)
    %    ['## error: true zeroes present in data; mean will fail ##']
    %    crasj = yes
    %end


    % put all results into a neat structure  
    R.PROFILE.amplitude = amplitudeSEs;
    R.ERROR.amplitude = semAmplSEs;
    
    R.PROFILE.positiveAmplitude = positiveAmplitudeSEs;
    R.ERROR.positiveAmplitude = semPosAmplSEs;
    
    R.PROFILE.negativeAmplitude = negativeAmplitudeSEs;
    R.ERROR.negativeAmplitude = semNegAmplSEs;
    
    R.PROFILE.latencyToPeak = latencyToPeakSEs;
    R.ERROR.latencyToPeak = [];%just to make plotting a bit easier downstream
    
    R.PROFILE.latencyToTrough = latencyToTroughSEs;
    R.ERROR.latencyToTrough = [];

    R.PROFILE.transect = transectSEs;
    R.ERROR.transect = semTransectSEs;
    R.transectTime = transectTime;

    R.PROFILE.avTransectWindow = avTransectWindowSEs;
    R.ERROR.avTransectWindow = semAvTransectWindowSEs;
    R.avTransectWindow = avTransectWindow;
    
    R.meanERPs = meanERPs;
    R.nERPs = nERPs;
    
    R.ind_max_erp = ind_max_erp;
    R.ind_min_erp = ind_min_erp;
    
    R.allERPs = allERPs;
    R.semERPs = semERPs;

    % plot ERPs for each sequence separately in a 4x4 plot
    % for each sequence, highlight where the maxima (red) and minima (blue) are located 
    if aux_plots && ~isempty(allPHOTs)
        
        figure;

        switch n_back
            case 5
                load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');
                labels = binomial_x_labels_latex;
            otherwise
                loadName = [num2str(n_back),'-back_legend.mat'];
                eval(['load ',loadName])
                %labels = anynomial_x_labels_latex; %Old style; Native ordering
                labels = anynomial_x_labels_latex_canonical; %Matches what is applied by seq_eff_order in analyseSequentialEffects
        end

        %this is just to help turn horizontal sequences into vertical ones
        %ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]);
        ind_horiz = sub2ind(size(labels{1}),1:n_back-1,[ones(1,n_back-2) 5]);
        %labels{i}(ind_horiz)

        for i = 1:nActive
           %subplot(4,4,i);
           subplot(ceil(sqrt(0.5*2^n_back)),ceil(sqrt(0.5*2^n_back)),i);
           plot(normalize(meanERPs(:,i)));
           hold on;
           scatter(ind_max_erp(i), 0,40,'r','filled');
           scatter(ind_min_erp(i), 0,40,'b','filled');
           plot(normalize(meanPHOTs(:,i)));
%            plot([window(1) window(1)],ylim,'r');
           %title(binomial_x_labels_latex{i}(ind_horiz));
           title(labels{i}(ind_horiz));
        end
    end
    
    %Matt extra plots
    if aux_plots
        figure
        hold on
        for seq = 1:size(meanERPs,2)
            plot(meanERPs(:,seq))
        end
        yLim = get(gca,'YLim');

        if plotSelector(3)
        temp = nanmedian(ind_min_erp,'all');
        temp2 = nanmedian(min_erp,'all');
        line([temp,temp], [temp2*0.9,temp2*1.1],'LineStyle','--','Color','m')
        text(temp,temp2*0.9,'Min','Color','m')
        end

        if plotSelector(2)
        temp = nanmedian(ind_max_erp,'all');
        temp2 = nanmedian(max_erp,'all');
        line([temp,temp], [temp2*0.9,temp2*1.1],'LineStyle','--','Color','b')
        text(temp,temp2*0.9,'Max','Color','b')
        end

        if plotSelector(1)
        temp = nanmean( [nanmedian(ind_min_erp,'all'),nanmedian(ind_max_erp,'all')] );
        line([temp,temp],[nanmedian(min_erp,'all'),nanmedian(max_erp,'all')],'Color','k')
        text(temp,nanmedian(max_erp,'all'),'Amp.','Color','b')
        end

        if plotSelector(6)
        line([transectTime,transectTime], [yLim(1),yLim(2)],'LineStyle','--','Color','c')
        text(transectTime,yLim(2),'Transect','Color','c')
        end

        if ~isempty(avTransectWindow) && plotSelector(7)
            line([avTransectWindow(1),avTransectWindow(1)], [yLim(1),yLim(2)],'LineStyle','--','Color','g')
            line([avTransectWindow(2),avTransectWindow(2)], [yLim(1),yLim(2)],'LineStyle','--','Color','g')
            text(transectTime,nanmean(yLim),'Av. transect window','Color','g')
        end

    end

end

