%% Extra transect analyses
function extraTransectAnalyses(FLIES, chosenFlies, resultsDirectory,options)

arguments
    FLIES struct
    chosenFlies double
    resultsDirectory char
    options.alphaVal double = 0.05 % p
    options.correctForNTimepoints double = 1 %Whether to perform conservation division of P-value for how many ANOVAs will be performed (in addition to Bonferroni within ANOVAs)
    options.showPValPlots double = 0
    options.patchMethod char = 'average'
    options.n_back double = 5
    options.reOrder double = [1:16]
    options.plotIndividualFlies double = 0;
end

alphaVal = options.alphaVal;
correctForNTimepoints = options.correctForNTimepoints;
showPValPlots = options.showPValPlots;
patchMethod = options.patchMethod;
n_back = options.n_back;
reOrder = options.reOrder;
plotIndividualFlies = options.plotIndividualFlies;

%%

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

resultsDirectory = [resultsDirectory '\Transect\'];
if exist(resultsDirectory) ~= 7
    mkdir(resultsDirectory)
    disp(['Transect results directory had to be made'])
    disp(resultsDirectory)
end

%%

%All Timepoint ANOVA
%Per fly
for fly = 1:size(chosenFlies,2)
    thisFly = chosenFlies(fly);
    %disp(['Now analysing fly #',num2str(thisFly),' for extra transect analyses'])
    thisFlyData = FLIES(thisFly);

    %QA for n-back
    if size(thisFlyData.nERPs,2) ~= 0.5*2^n_back
        ['## Alert: Potentially incorrect n_back being used for extra transect analyses ##']
        crash = yes
        %Mostly a concern for figure labels/etc
    end

    %Testatory figure
    %{
    figure
    hold on
    for seq = 1:size(thisFlyData.meanERPs,2)
        plot(thisFlyData.meanERPs(:,seq))
    end
    title([thisFlyData.date,' - B',num2str(thisFlyData.block),' mean ERPs'])
    %}

    %Decide on P value
    if ~correctForNTimepoints
        effectiveAlpha = alphaVal;
    else
        effectiveAlpha = alphaVal / size( thisFlyData.allERPs, 1 );
    end

    pVals = [];
    for timep = 1:size( thisFlyData.allERPs, 1 )
        [pVals(timep),ANOVATAB,STATS] = anova1( squeeze( thisFlyData.allERPs(timep,:,:) ).', [], 'off' ); %Need to use all ERPs data for ANOVA because n
    end

    %sigPatches = bwlabel( pVals < alphaVal );
    sigPatches = bwlabel( pVals < effectiveAlpha );
    numSig = nanmax(sigPatches);

    %And another
    if showPValPlots && plotIndividualFlies
        figure
        subplot(3,1,1)
        hold on
        for seq = 1:size(thisFlyData.meanERPs,2)
            plot(thisFlyData.meanERPs(:,seq))
        end
        title([thisFlyData.date,' - B',num2str(thisFlyData.block),' mean ERPs'])
        subplot(3,1,2)
        plot(pVals)
        hold on
        %line([0,size(pVals,2)],[effectiveAlpha,effectiveAlpha],'LineStyle',':','Color','k')
        title(['(',num2str(thisFly),') ',thisFlyData.date,' - B',num2str(thisFlyData.block),' inter-seq. transect ANOVA P-values'])
        subplot(3,1,3)
        plot(pVals)
        hold on
        line([0,size(pVals,2)],[effectiveAlpha,effectiveAlpha],'LineStyle',':','Color','k')
        ylim([0,effectiveAlpha*1.5])
        title(['P-values zoom, wrt p<',num2str(effectiveAlpha)])
        saveas(gcf,[ resultsDirectory '/' 'PVal' '_fly' num2str(thisFly) '.png']);
    end

    %Profile/s for sig over time
    sigProfiles = {};
    if numSig > 0 && plotIndividualFlies
        for sigpa = 1:nanmax(sigPatches)
            switch patchMethod
                case 'average'
                    thesePInd = find( sigPatches == sigpa ); %Both indices and logicals work here
                    sigDataAll = squeeze( thisFlyData.allERPs( thesePInd , : , : ) ); %Time x Seq x Instance; Absed
                        %Note that at any one 'instance', there will only be one non-NaN Time x Seq, due to sparse nature
                    %sigData = nansum( sigDataAll, 3 ); %Time x Seq [Instances summed]; Note: Not technically incorrect, but probably not correct
                    %sigData = nanmean( sigDataAll, 3 ); %Time x Seq [Instances averaged]
                    sigData = abs( nanmean( sigDataAll, 3 ) ); %Time x Seq [Instances averaged]
                    meanData = nanmean( sigData, 1 );
                        %Note: Profiles 'averaged' by virtue of underlying numbers, not profiles themselves
                    underN = size(sigData,1);
                    meanSEM = nanstd( sigData, [], 1 ) / sqrt( underN );
                        %Note: SEM here represents error across time (e.g. 25 frames)
                case 'lowestP'
                    [theseP, thesePInd] = min( pVals( sigPatches == sigpa ) );
                    thesePInd = thesePInd + find(sigPatches == sigpa,1,'first') - 1; %Necessary because thesePInd originally in patch-specific reference frame
                    sigData = abs( squeeze( thisFlyData.allERPs( thesePInd , : , : ) ) )'; %Instance x Seq; Absed
                        %Sparse as above, now manifesting as one element per row (rather than layer)
                    %meanData = nanmean(sigData,1);
                    meanData = nanmean(sigData,1); %Apply abs
                    underN = nansum( ~isnan( sigData ) , 1 );
                    %meanSEM = nanstd( sigData, [], 1 ) ./ sqrt( underN );
                    meanSEM = nanstd( sigData, [], 1 ) ./ sqrt( underN );
                        %Note: Different to above, SEM here represents # of actual instances of seq (e.g. 210)

            end

            %Individual figure
            %{
            figure
            errorbar(blerg,blergSEM)
            %hold on
            %indivPos = floor(size(blorg,1)/2);
            %errorbar(blorg( indivPos, : ),blergSEM) %Plot just one transect
            xticks([1:size(blerg,2)])
            xticklabels(exLabels)
            xtickangle(270)
            title(['Sig. patch #',num2str(sigpa),' averaged ERP profile'])
            %}

            sigProfiles{sigpa}{1} = meanData; %Mean profile
            sigProfiles{sigpa}{2} = meanSEM; %SEM of above
            sigProfiles{sigpa}{3} = [ unique( [min(thesePInd),max(thesePInd)] ) ]; %Lowest/Highest index for sampled area (Unitary if lowestP)

        end
    %end

        %Display of sig across time
        figure
        subplot(2,numSig,[1:numSig]) %subplot across first row
        hold on
        for seq = 1:size(thisFlyData.meanERPs,2)
            plot(thisFlyData.meanERPs(:,seq))
        end
        yLim = get(gca,'YLim');
        safeAlt = max( yLim );
        if numSig > 0
            for sigpa = 1:numSig
                startEnd = [ find(sigPatches == sigpa,1,'first'), find(sigPatches == sigpa,1,'last') ];
                line([startEnd],[safeAlt,safeAlt],'LineWidth',2,'Color','k')
                %text(nanmean(startEnd),safeAlt-20,num2str(sigpa),'Color','k')
                text(nanmean(startEnd),nanmax(yLim)*0.9,num2str(sigpa),'Color','k')
            end
        end
        if safeAlt < 0
            ylim([min(yLim),safeAlt*0.9])
        else
            ylim([min(yLim),safeAlt*1.1]) %Too lazy to perform sign-invariant y limit calcs
        end
        title(['(#',num2str(thisFly),') ',thisFlyData.date,' - B',num2str(thisFlyData.block),' sig. (p<',num2str(effectiveAlpha),') across time'])
        %Quick QA
        if max(get(gca,'YLim')) < safeAlt
            ['-# Caution: Sig. displayed outside figure Y limits! #-']
            ylim([min(yLim),safeAlt])
        end
        %Sig profiles
        %if numSig > 0
        for sigpa = 1:numSig
            subplot(2,numSig,numSig+sigpa)
            %errorbar(sigProfiles{sigpa}{1},sigProfiles{sigpa}{2})
            errorbar(sigProfiles{sigpa}{1}(reOrder),sigProfiles{sigpa}{2}(reOrder))
            xticks([1:size(meanData,2)]) %Note: Not inherently reordered
            %xticklabels(exLabels)
            xticklabels(exLabels(reOrder))
            xtickangle(270)
            xlim([0,size(meanData,2)+1])
            title(['Sig. patch #',num2str(sigpa),' ',patchMethod,' (Ind:',num2str(sigProfiles{sigpa}{3}),') profile (n per seq = ~',num2str(floor( nanmean(underN)) ),')'])
        end
        %end
        saveas(gcf,[ resultsDirectory '/' 'TransectSigWProf' '_fly' num2str(thisFly) '.png']);
    end

end

end