function ERROR = calculateError(FLIES, chosenFlies, reOrder, n_back, options)

arguments
    FLIES struct
    chosenFlies double
    reOrder double
    n_back double
    options.suppressANOVA double = 0
end

%Extra
suppressANOVA = options.suppressANOVA;

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

%amplitudeSEs = zeros(16,length(chosenFlies));
%positiveAmplitudeSEs = zeros(16,length(chosenFlies));
%negativeAmplitudeSEs = zeros(16,length(chosenFlies));
%transectSEs = zeros(16,length(chosenFlies));
%avTransectWindowSEs = zeros(16,length(chosenFlies));
amplitudeSEs = zeros(0.5*2^n_back,length(chosenFlies));
positiveAmplitudeSEs = zeros(0.5*2^n_back,length(chosenFlies));
negativeAmplitudeSEs = zeros(0.5*2^n_back,length(chosenFlies));
transectSEs = zeros(0.5*2^n_back,length(chosenFlies));
avTransectWindowSEs = zeros(0.5*2^n_back,length(chosenFlies));

%collect profiles in matrix
for fly = 1:length(chosenFlies)        
    %sequential effects results
    amplitudeSEs(:,fly) = FLIES(chosenFlies(fly)).PROFILE.amplitude.';
    positiveAmplitudeSEs(:,fly) = FLIES(chosenFlies(fly)).PROFILE.positiveAmplitude.';
    negativeAmplitudeSEs(:,fly) = FLIES(chosenFlies(fly)).PROFILE.negativeAmplitude.';
    transectSEs(:,fly) = FLIES(chosenFlies(fly)).PROFILE.transect.';
    avTransectWindowSEs(:,fly) = FLIES(chosenFlies(fly)).PROFILE.avTransectWindow.';
end

%testatory qa
    %Defunct for same reasons as stated in calculateSEs (i.e. Zeroes are likely true, for posAmp / negAmp at least)
%if any(amplitudeSEs == 0) |any(positiveAmplitudeSEs == 0) |any(negativeAmplitudeSEs == 0) |any(transectSEs == 0) |any(avTransectWindowSEs == 0)
%    ['## error: true zeroes present in data; mean will fail ##']
%    crasj = yes
%end

%subtract mean
amplitudeSEs = amplitudeSEs-mean(amplitudeSEs);
positiveAmplitudeSEs = positiveAmplitudeSEs-mean(positiveAmplitudeSEs);
negativeAmplitudeSEs = negativeAmplitudeSEs-mean(negativeAmplitudeSEs);
transectSEs = transectSEs-mean(transectSEs);
avTransectWindowSEs = avTransectWindowSEs-mean(avTransectWindowSEs);

%onion plot of profiles
%figure; create_seq_eff_plot(amplitudeSEs,[]);
figure; create_seq_eff_plot(amplitudeSEs,[],'reOrder',reOrder, 'n_back',n_back,'histlength',n_back-1);
title(['Cross-fly amps profiles'])

if ~suppressANOVA
    %do ANOVA on subtracted mean data (for amplitude)
    %[p,ANOVATAB,STATS] = anova1(amplitudeSEs.');
    [p,ANOVATAB,STATS] = anova1(amplitudeSEs([reOrder],:).');
    %title(['Cross-fly amplitude ANOVA (p=',num2str(p),')',newline,newline])
    titleStr = ['Cross-fly Amplitude ANOVA (p=',num2str(p),')'];
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
    xticks([1:size(transectSEs,1)])
    xticklabels(exLabels(reOrder)) %Ditto criticality
    xtickangle(270)
    
    %Repeat for transect
    %figure; create_seq_eff_plot(transectSEs,[]);
    figure; create_seq_eff_plot(transectSEs,[],'reOrder',reOrder,'n_back',n_back,'histlength',n_back-1);
    title(['Cross-fly transect profiles'])
    %[p,ANOVATAB,STATS] = anova1(transectSEs.');
    [p,ANOVATAB,STATS] = anova1(transectSEs([reOrder],:).');
    %title(['Cross-fly transect ANOVA (p=',num2str(p),')',newline,newline])
    titleStr = ['Cross-fly Transect ANOVA (p=',num2str(p),')'];
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
    xticks([1:size(transectSEs,1)])
    xticklabels(exLabels(reOrder))
    xtickangle(270)
end


semAmplSEs = std(amplitudeSEs,[],2)/sqrt(length(chosenFlies));
semPosAmplSEs = std(positiveAmplitudeSEs,[],2)/sqrt(length(chosenFlies));
semNegAmplSEs = std(negativeAmplitudeSEs,[],2)/sqrt(length(chosenFlies));
semTransSEs = std(transectSEs,[],2)/sqrt(length(chosenFlies));
semAvTransWinSEs = std(avTransectWindowSEs,[],2)/sqrt(length(chosenFlies));

%Combined seq dep plot/s
figure; create_seq_eff_plot(mean(amplitudeSEs,2),[],'errors',semAmplSEs, 'reOrder',reOrder,'n_back',n_back,'histlength',n_back-1);
title('Cross-fly combined amplitude profile')
figure; create_seq_eff_plot(mean(transectSEs,2),[],'errors',semTransSEs,'reOrder',reOrder,'n_back',n_back,'histlength',n_back-1);
title('Cross-fly combined transect profile')

%propagated errors
ERROR.amplitude = semAmplSEs.';
ERROR.positiveAmplitude = semPosAmplSEs.';
ERROR.negativeAmplitude = semNegAmplSEs.';
ERROR.latencyToPeak = [];
ERROR.latencyToTrough = [];
ERROR.transect = semTransSEs.';
ERROR.avTransectWindow = semAvTransWinSEs.';