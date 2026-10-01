function SEAnova(R, reOrder, n_back, transProbAncillary)
%SEAnova Performs an ANOVA on SE data

if isempty(transProbAncillary)
    transProbDesign = 0;
else
    transProbDesign = 1;
end

    if ~transProbDesign

        nBackActual = n_back;

        %Matthew system for turning Dinis X labels into something conveniently usable
        switch nBackActual
            case 5
                load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');
                labels = binomial_x_labels_latex;
            otherwise
                loadName = [num2str(nBackActual),'-back_legend.mat'];
                eval(['load ',loadName])
                %labels = anynomial_x_labels_latex; %Old style; Native ordering
                labels = anynomial_x_labels_latex_canonical; %Matches what is applied by seq_eff_order in analyseSequentialEffects
        end
        %this is just to help turn horizontal sequences into vertical ones
        %ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]); %Hardcoded n-back of 5
        ind_horiz = sub2ind(size(labels{1}),1:nBackActual-1,[ones(1,nBackActual-2) 5]); %Dynamic
        exLabels = [];
        for s = 1:size(labels,2)
            %exLabels{s} = binomial_x_labels_latex{s}(ind_horiz);
            exLabels{s} = labels{s}(ind_horiz);
        end
        reOrderActual = reOrder;

    else

        nBackActual = transProbAncillary.nBackActual;

        exLabels = [];
        for i = 1:size(transProbAncillary.transLabels,1)
            exLabels{i} = transProbAncillary.transLabels(i,:);
        end   
        reOrderActual = transProbAncillary.reOrderActual;

    end

    allERPs = R.allERPs;
    semERPs = R.semERPs;
        
    if ~transProbDesign

        % full data for maxima and minima across all stimuli
        aux = squeeze(reshape(allERPs,[1,size(allERPs,1)*size(allERPs,2),size(allERPs,3)])).';
            %Collapses allERPs from T x Seq x Event(Sparse) to Event(Sparse) x TSeq[?]
                %Thus, each event is represented on one row of aux, where column is T[?]

        % data tables for performing an ANOVA to check for the effect of
        % sequence
        %dataMaxima = aux(:,sub2ind(size(semERPs),R.ind_max_erp,1:16));
        %dataMinima = aux(:,sub2ind(size(semERPs),R.ind_min_erp,1:16));
        dataMaxima = aux(:,sub2ind(size(semERPs),R.ind_max_erp,1:0.5*2^nBackActual));
        dataMinima = aux(:,sub2ind(size(semERPs),R.ind_min_erp,1:0.5*2^nBackActual));
        dataAmplitude = dataMaxima-dataMinima;

    else
        %Old, aux-based system that fails for trans prob
        %aux = squeeze(reshape(allERPs,[1,size(allERPs,1)*size(allERPs,2),size(allERPs,3)])).';   
        %dataMaxima = aux(:,sub2ind(size(semERPs),R.ind_max_erp, size(allERPs,2) )); %Replace nBack exponential with size query
        %dataMinima = aux(:,sub2ind(size(semERPs),R.ind_min_erp, size(allERPs,2) ));
        %dataAmplitude = dataMaxima-dataMinima;


        %dataMaxima = allERPs( R.ind_max_erp , : , : ); %
        dataMaxima = squeeze(allERPs([R.ind_max_erp],:,:));
        dataMinima = squeeze(allERPs([R.ind_min_erp],:,:));
            %Yields an nStimuli*nBackActual x nStimuli*nBackActual x Event # matrix
                %'Duplication' of 'sequence' identity is because sparse in column *and* layer

        dataAmplitude = dataMaxima - dataMinima; %Subtract at a per-event level (Why sparse is tolerated)
        dataAmplitude = squeeze( nansum( dataAmplitude, 2) )'; %Collapse to an Event x Seq matrix as done for SEs
        
        ['## Alert: Math currently incorrect; Please suppress ANOVAs for now for trans probs ##']
        crash = yes


    end
    
    %[p, ANOVATAB,STATS]= anova1(dataAmplitude);
    [p, ANOVATAB,STATS]= anova1(dataAmplitude(:,reOrderActual));
    disp(['Amplitude p value: ',num2str(p)])
    %title(['Amplitude ANOVA (p=',num2str(p),')'])
    %title([R.date,' B',R.block,' Amplitude ANOVA (p=',num2str(p),')',newline,newline])
    titleStr = [R.date,' B',R.block,' Amplitude ANOVA (p=',num2str(p),')'];
    for surp = 1:nBackActual - 3
        titleStr = [titleStr,newline]; %Add an appropriate number of newlines
    end
    title(titleStr)
    %Utterly overengineered system to place SEs labelling on top axis
    xticklabels(reOrderActual) %Critical if reordering used
    ax1 = gca;
    ax2 = axes('Position', ax1.Position, ...
                   'XAxisLocation', 'top', ...
                   'Color', 'none', ...
                   'YTick', [], ...
                   'Box', 'off');
    linkaxes([ax1, ax2], 'x');
    xticks([1:size(dataAmplitude,2)])
    xticklabels(exLabels(reOrderActual))
    xtickangle(270)


    dataTransect = squeeze(allERPs(R.transectTime,:,:))';

    %[p, ANOVATAB,STATS]= anova1(dataTransect);
    [p, ANOVATAB,STATS]= anova1(dataTransect(:,reOrderActual));
    disp(['Transect p value: ',num2str(p)])
    %title(['Transect ANOVA (p=',num2str(p),')'])
    %figure
    %multcompare(STATS)
    %title([R.date,' B',R.block,' Transect ANOVA (p=',num2str(p),')',newline,newline])
    titleStr = [R.date,' B',R.block,' Transect ANOVA (p=',num2str(p),')'];
    for surp = 1:nBackActual - 3
        titleStr = [titleStr,newline]; %Add an appropriate number of newlines
    end
    title(titleStr)
    
    %Utterly overengineered system to place SEs labelling on top axis
    xticklabels(reOrderActual) %Critical if reordering used
    ax1 = gca;
    ax2 = axes('Position', ax1.Position, ...
                   'XAxisLocation', 'top', ...
                   'Color', 'none', ...
                   'YTick', [], ...
                   'Box', 'off');
    linkaxes([ax1, ax2], 'x');
    xticks([1:size(dataTransect,2)])
    xticklabels(exLabels(reOrderActual))
    xtickangle(270)

end

