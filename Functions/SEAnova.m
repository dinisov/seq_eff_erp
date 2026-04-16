function SEAnova(R, reOrder, n_back)
%SEAnova Performs an ANOVA on SE data

    %Matthew system for turning Dinis X labels into something conveniently usable
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
    %ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]); %Hardcoded n-back of 5
    ind_horiz = sub2ind(size(labels{1}),1:n_back-1,[ones(1,n_back-2) 5]); %Dynamic
    exLabels = [];
    for s = 1:size(labels,2)
        %exLabels{s} = binomial_x_labels_latex{s}(ind_horiz);
        exLabels{s} = labels{s}(ind_horiz);
    end

    allERPs = R.allERPs;
    semERPs = R.semERPs;
    
    % full data for maxima and minima across all stimuli
    aux = squeeze(reshape(allERPs,[1,size(allERPs,1)*size(allERPs,2),size(allERPs,3)])).';
    
    % data tables for performing an ANOVA to check for the effect of
    % sequence
    %dataMaxima = aux(:,sub2ind(size(semERPs),R.ind_max_erp,1:16));
    %dataMinima = aux(:,sub2ind(size(semERPs),R.ind_min_erp,1:16));
    dataMaxima = aux(:,sub2ind(size(semERPs),R.ind_max_erp,1:0.5*2^n_back));
    dataMinima = aux(:,sub2ind(size(semERPs),R.ind_min_erp,1:0.5*2^n_back));
    dataAmplitude = dataMaxima-dataMinima;
    
    %[p, ANOVATAB,STATS]= anova1(dataAmplitude);
    [p, ANOVATAB,STATS]= anova1(dataAmplitude(:,reOrder));
    disp(['Amplitude p value: ',num2str(p)])
    %title(['Amplitude ANOVA (p=',num2str(p),')'])
    %title([R.date,' B',R.block,' Amplitude ANOVA (p=',num2str(p),')',newline,newline])
    titleStr = [R.date,' B',R.block,' Amplitude ANOVA (p=',num2str(p),')'];
    for surp = 1:n_back - 3
        titleStr = [titleStr,newline]; %Add an appropriate number of newlines
    end
    title(titleStr)
    %Utterly overengineered system to place SEs labelling on top axis
    xticklabels(reOrder) %Critical if reordering used
    ax1 = gca;
    ax2 = axes('Position', ax1.Position, ...
                   'XAxisLocation', 'top', ...
                   'Color', 'none', ...
                   'YTick', [], ...
                   'Box', 'off');
    linkaxes([ax1, ax2], 'x');
    xticks([1:size(dataAmplitude,2)])
    xticklabels(exLabels(reOrder))
    xtickangle(270)


    dataTransect = squeeze(allERPs(R.transectTime,:,:))';

    %[p, ANOVATAB,STATS]= anova1(dataTransect);
    [p, ANOVATAB,STATS]= anova1(dataTransect(:,reOrder));
    disp(['Transect p value: ',num2str(p)])
    %title(['Transect ANOVA (p=',num2str(p),')'])
    %figure
    %multcompare(STATS)
    %title([R.date,' B',R.block,' Transect ANOVA (p=',num2str(p),')',newline,newline])
    titleStr = [R.date,' B',R.block,' Transect ANOVA (p=',num2str(p),')'];
    for surp = 1:n_back - 3
        titleStr = [titleStr,newline]; %Add an appropriate number of newlines
    end
    title(titleStr)
    
    %Utterly overengineered system to place SEs labelling on top axis
    xticklabels(reOrder) %Critical if reordering used
    ax1 = gca;
    ax2 = axes('Position', ax1.Position, ...
                   'XAxisLocation', 'top', ...
                   'Color', 'none', ...
                   'YTick', [], ...
                   'Box', 'off');
    linkaxes([ax1, ax2], 'x');
    xticks([1:size(dataTransect,2)])
    xticklabels(exLabels(reOrder))
    xtickangle(270)

end

