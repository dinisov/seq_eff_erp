%% Extra transect analyses
function [CROSSISOMER] = extraIsomerAnalyses(FLIES, chosenFlies, resultsDirectory,options)

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
    options.plotSelector double = [];
    options.plotIndividualFlies double = 1;
    options.timeStep double = 20; %What timestep to use for the isomer x time plot
    options.doAnimatedPlot double = 0; %Whether to do a once-off animated plot that rolls across all timepoints (Note: Only done if at least one transect-type analysis requested in plotSelector)
    options.animFrameRate double = 10; %Tentative attempted framerate for animated plot to display/save at
    options.saveFigVid double = 0; %Whether to save the animation to the figure folder
    options.extraFigSubplots double = 0; %Whether to add ERP subplots to animated plot (1) or rearrange (2)
end

alphaVal = options.alphaVal;
correctForNTimepoints = options.correctForNTimepoints;
showPValPlots = options.showPValPlots;
patchMethod = options.patchMethod;
n_back = options.n_back;
reOrder = options.reOrder;
plotSelector = options.plotSelector;
plotIndividualFlies = options.plotIndividualFlies;
timeStep = options.timeStep;
doAnimatedPlot = options.doAnimatedPlot;
animFrameRate = options.animFrameRate;
saveFigVid = options.saveFigVid;
extraFigSubplots = options.extraFigSubplots;

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
    %Pre-data QA
    if ~isfield(thisFlyData,'ISOMER')
        ['-# Alert: Fly #',num2str(thisFly),' lacks stored ISOMER data #-']
        continue
    end
    numIsomers = size(fieldnames(thisFlyData.ISOMER),1); %Futureproofing?
    if numIsomers > 2
        ['## Caution: Code case not written for >2 isomers yet ##']
        crash = yes %cbf right now writing true generalised (and combinatorial) code for multiple isomer comparison
    end
    R1 = thisFlyData.ISOMER.R1;
    R2 = thisFlyData.ISOMER.R2;

    %QA for n-back
    if size(thisFlyData.nERPs,2) ~= 0.5*2^n_back
        ['## Alert: Potentially incorrect n_back being used for extra isomer analyses ##']
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

    rVals = [];
    pVals = [];
    for timep = 1:size( thisFlyData.allERPs, 1 )
        %[pVals(timep),ANOVATAB,STATS] = anova1( squeeze( thisFlyData.allERPs(timep,:,:) ).', [], 'off' ); %Need to use all ERPs data for ANOVA because n
        [R,P] = corrcoef(R1.meanERPs(timep,:),R2.meanERPs(timep,:)); %Note: Will crash if allERPs and meanERPs different size, but bigger issues if so
        rVals(timep) = R(2,1); %ONLY VALID AS LONG AS ONLY TWO ISOMERS BEING COMPARED
        pVals(timep) = P(2,1);
    end

    %sigPatches = bwlabel( pVals < alphaVal );
    sigPatches = bwlabel( pVals < effectiveAlpha );
    numSig = nanmax(sigPatches);

    %And another
    if showPValPlots && plotIndividualFlies
        figure
        subplot(3,1,1)
        hold on
        %for seq = 1:size(thisFlyData.meanERPs,2)
        %    plot(thisFlyData.meanERPs(:,seq))
        %end
        plot(rVals,'Color','m')
        title([thisFlyData.date,' - B',num2str(thisFlyData.block),' R vals'])
        subplot(3,1,2)
        plot(pVals,'Color','m')
        hold on
        %line([0,size(pVals,2)],[effectiveAlpha,effectiveAlpha],'LineStyle',':','Color','k')
        title(['(',num2str(thisFly),') ',thisFlyData.date,' - B',num2str(thisFlyData.block),' inter-profile corr. P-values'])
        subplot(3,1,3)
        plot(pVals,'Color','m')
        hold on
        line([0,size(pVals,2)],[effectiveAlpha,effectiveAlpha],'LineStyle',':','Color','k')
        ylim([0,effectiveAlpha*1.5])
        title(['P-values zoom, wrt p<',num2str(effectiveAlpha)])
        saveas(gcf,[ resultsDirectory '/' 'RCoeffPVal' '_fly' num2str(thisFly) '.png']);
    end

    %Profile/s for sig over time
    sigProfiles = {};
    if numSig > 0
        for sigpa = 1:nanmax(sigPatches)
            meanData = [];
            switch patchMethod
                case 'average'
                    thesePInd = find( sigPatches == sigpa ); %Both indices and logicals work here
                    %sigDataAll = squeeze( thisFlyData.allERPs( thesePInd , : , : ) ); %Time x Seq x Instance; Absed
                    %    %Note that at any one 'instance', there will only be one non-NaN Time x Seq, due to sparse nature
                    %sigData = abs( nanmean( sigDataAll, 3 ) ); %Time x Seq [Instances averaged]
                    meanData{1} = nanmean(R1.meanERPs(thesePInd,:),1); 
                    meanData{2} = nanmean(R2.meanERPs(thesePInd,:),1); 
                    %meanData = nanmean( sigData, 1 );
                    %    %Note: Profiles 'averaged' by virtue of underlying numbers, not profiles themselves
                    %underN = size(sigData,1);
                    %meanSEM = nanstd( sigData, [], 1 ) / sqrt( underN );
                        %Note: SEM here represents error across time (e.g. 25 frames)
                    %As below wrt not using allERPs
                    rData = nanmean(rVals(thesePInd));
                    pData = nanmean(pVals(thesePInd));
                case 'lowestP'
                    [theseP, thesePInd] = min( pVals( sigPatches == sigpa ) );
                    thesePInd = thesePInd + find(sigPatches == sigpa,1,'first') - 1; %Necessary because thesePInd originally in patch-specific reference frame
                    %sigData = abs( squeeze( thisFlyData.allERPs( thesePInd , : , : ) ) )'; %Instance x Seq; Absed
                    %    %Sparse as above, now manifesting as one element per row (rather than layer)
                    %meanData = nanmean(sigData,1); %Apply abs
                    meanData{1} = nanmean(R1.meanERPs(thesePInd,:),1); 
                    meanData{2} = nanmean(R2.meanERPs(thesePInd,:),1); 
                    %underN = nansum( ~isnan( sigData ) , 1 );
                    %meanSEM = nanstd( sigData, [], 1 ) ./ sqrt( underN );
                        %Note: Different to above, SEM here represents # of actual instances of seq (e.g. 210)
                    %Note: Since errorless meanERP is what is fed to corrcoef, this is used in this plot too
                        %With a little divergence, the actual underlying error could be shown
                    rData = nanmean(rVals(thesePInd));
                    pData = nanmean(pVals(thesePInd));
            end

            %Individual figure
            %{
            figure
            %errorbar(blerg,blergSEM)
            plot(meanData{1})
            hold on
            plot(meanData{2})
            %indivPos = floor(size(blorg,1)/2);
            %errorbar(blorg( indivPos, : ),blergSEM) %Plot just one transect
            xticks([1:size(meanData{1},2)])
            xticklabels(exLabels)
            xtickangle(270)
            title(['R sig. patch #',num2str(sigpa),' isomers (Not reordered)'])
            %}

            sigProfiles{sigpa}{1} = meanData; %Mean profile
            %sigProfiles{sigpa}{2} = meanSEM; %SEM of above
            sigProfiles{sigpa}{3} = [ unique( [min(thesePInd),max(thesePInd)] ) ]; %Lowest/Highest index for sampled area (Unitary if lowestP)
            sigProfiles{sigpa}{4} = rData;
            sigProfiles{sigpa}{5} = pData;
        end
    %end

        if plotIndividualFlies
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
                %errorbar(sigProfiles{sigpa}{1}(reOrder),sigProfiles{sigpa}{2}(reOrder))
                plot(sigProfiles{sigpa}{1}{1}(reOrder),'Color','b')
                hold on
                plot(sigProfiles{sigpa}{1}{2}(reOrder),'Color','r')
                xticks([1:size(sigProfiles{sigpa}{1}{1},2)]) 
                %xticklabels(exLabels)
                xticklabels(exLabels(reOrder))
                xtickangle(270)
                xlim([0,size(sigProfiles{sigpa}{1}{1},2)+1])
                title(['Sig. patch #',num2str(sigpa),' ',patchMethod,' (Ind:',num2str(sigProfiles{sigpa}{3}),') profile (Mean R value: ',num2str(sigProfiles{sigpa}{4}),' [p=',num2str(sigProfiles{sigpa}{5}),'])'])
                legend({['R1'],['R2']})
            end
            %end
            saveas(gcf,[ resultsDirectory '/' 'RCoeffSigWIsom' '_fly' num2str(thisFly) '.png']);
        end

    end

    %masses

    if plotIndividualFlies
        %Plot of all isomers across time
        isomTime = struct;
        a = 1;
        for timep = 1:timeStep:size( thisFlyData.allERPs, 1 )
            isomTime(a).timep = timep;
            isomTime(a).data = [R1.meanERPs(timep,:);R2.meanERPs(timep,:)];
            isomTime(a).mean = [nanmean(squeeze(R1.allERPs(timep,:,:)),2)';...
                                nanmean(squeeze(R2.allERPs(timep,:,:)),2)'];
            %isomTime(a).SEM = [nanstd(squeeze(R1.allERPs(timep,:,:)),[],2)' ./ nansum( ~isnan( squeeze(R1.allERPs(timep,:,:)) ), 2 )';... %A REMINDER OF MISTAKES
            %                   nanstd(squeeze(R2.allERPs(timep,:,:)),[],2)' ./ nansum( ~isnan( squeeze(R2.allERPs(timep,:,:)) ), 2 )'];
            isomTime(a).SEM = [nanstd(squeeze(R1.allERPs(timep,:,:)),[],2)' ./ sqrt( nansum( ~isnan( squeeze(R1.allERPs(timep,:,:)) ), 2 ) )';...
                               nanstd(squeeze(R2.allERPs(timep,:,:)),[],2)' ./ sqrt( nansum( ~isnan( squeeze(R2.allERPs(timep,:,:)) ), 2 ) )'];
            a = a + 1;
        end

        figure
        for subl = 1:size(isomTime,2)
            subplot(1,size(isomTime,2),subl)
            %plot(isomTime(subl).data(1,reOrder),'Color','b')
            errorbar( isomTime(subl).mean(1,reOrder), isomTime(subl).SEM(1,reOrder), 'Color', 'b' )
            hold on
            %plot(isomTime(subl).data(2,reOrder),'Color','r')
            errorbar( isomTime(subl).mean(2,reOrder), isomTime(subl).SEM(2,reOrder), 'Color', 'r' )
            xticks([1:size(isomTime(subl).data,2)])
            xticklabels([])
            title(['T:',num2str(isomTime(subl).timep)])
        end
        legend({'R1','R2'})
        set(gcf,'Name',['IsomxTime' '_fly' num2str(thisFly)])
        saveas(gcf,[ resultsDirectory '/' 'IsomxTime' '_fly' num2str(thisFly) '.png']);
        
    end


end

%Cross-fly isomers (if applicable)
CROSSISOMER = [];
if size(chosenFlies,2) > 1
    if isempty(plotSelector)
        ['-# Cross-fly isomer analysis would be attempted, but no plotSelector specified #-']
        return
    end
    %Quick pre-grab of plot types
    plotTypes = fieldnames( FLIES(chosenFlies(1)).ISOMER.R1.PROFILE); %If this fails because R1 does not exist then you are an idiot for trying to run this
        %Note: Only 80% confident that fieldnames follows the known ordering, but this saves effort on hardcoding each plot type
    isomerNames = fieldnames( FLIES(chosenFlies(1)).ISOMER);
    isomerCount = size(isomerNames,1);
    isomerColours = jet(isomerCount);
    %Quick check of timepoints in ERPs, since it is important to be consistent
    temp = [];
    for fly = chosenFlies
        temp = [temp,size(FLIES(fly).meanERPs,1)];
    end
    if numel(unique(temp)) > 1
        ['## Alert: Dissimilar timepoints detected in ERPs across fly ##']
        [num2str(temp)]
    end
    nTimepoints = unique(temp);

    for plotI = find(plotSelector == 1)
        %unstop
        disp(['Calculating cross-fly isomer profile for ',plotTypes{plotI}])

        crossData = cell(isomerCount,1);
        if isequal( plotTypes{plotI} , 'avTransectWindow' ) | isequal( plotTypes{plotI} , 'transect' ) %Only do if transecting
            crossDataFull = cell(isomerCount,1); %Same as above, but contains full time information, and obtained differently
            transecting = 1;
        else
            transecting = 0;
        end
        crossMean = nan(0.5*2^n_back, isomerCount);
        crossSEM = nan(0.5*2^n_back, isomerCount);
        for isoI = 1:isomerCount
            crossData{isoI} = nan(0.5*2^n_back, size(chosenFlies,2)); %e.g. 16 seqs x 5 flies
            crossDataFull{isoI} = nan(nTimepoints, 0.5*2^n_back, size(chosenFlies,2)); %Note different architecture, for ease of data intake

            %Collect flies data
            for flyI = 1:size(chosenFlies,2)
                thisFly = chosenFlies(flyI);
                temp = FLIES(thisFly).ISOMER.(isomerNames{isoI}).PROFILE.(plotTypes{plotI});
                %                               e.g. R1                    e.g. transect
                %QA
                if any(temp == 0) && ( contains(plotTypes{plotI},'transect') ) %For positiveAmplitude (and maybe negativeAmplitude) this is 'okay', but not (presumably) for transect
                    disp(['-# Caution: Dataset #',num2str(flyI),' (Fly #',num2str(chosenFlies(flyI)),') contained 0s in data #-'])
                end
                %temp(temp == 0) = NaN;
                    %Note: Currently unsure how true zeroes are being managed elsewhere in code
                crossData{isoI}(:,flyI) = temp;
                if transecting
                    crossDataFull{isoI}(:, :, flyI) = FLIES(thisFly).ISOMER.(isomerNames{isoI}).meanERPs;
                end
            end

            %Perform meaning, as done by Dinis in calculateError/etc
            %crossData{isoI} = crossData{isoI} - mean(crossData{isoI});
            crossData{isoI} = crossData{isoI} - nanmean(crossData{isoI}); 
                %The intended behaviour here is that all profiles (1x16) are balanced around their own mean
            %Repeat for crossDataFull if applicable
            if transecting
                crossDataFull{isoI} = crossDataFull{isoI} - nanmean( crossDataFull{isoI}, 2 ); %Note: Dimension for meaning derived empirically
                %Testatory animated plot
                %{
                figure
                for t = 1:size(crossDataFull{isoI}, 1 )
                clf
                blirg = squeeze( crossDataFull{isoI}(t,:,:) );
                for fly = 1:size(blirg,2)
                plot( blirg(:,fly) )
                hold on
                end
                title('post auto meaning, t = ',num2str(t))
                pause(0.1)
                end
                %}
            end
            
            crossMean(:,isoI) = nanmean( crossData{isoI},2 );
            crossSEM(:,isoI) = nanstd( crossData{isoI}, [], 2 ) / sqrt( size(crossData{isoI},2) );

        end

        %body

        %Plot SEs with Dinis functions
            %Based on implementation in plotIsomers
        figure; 
        %create_seq_eff_plot([thisData1.' thisData2.'],[],'errors',[thisError1.' thisError2.'],...
        create_seq_eff_plot([crossMean],[],'errors',[crossSEM],...
            'reOrder', reOrder,'n_back',n_back,'histlength',n_back-1);
        %h = legend({'First Isomer','Second Isomer'});
        h = legend(isomerNames);
        %title([flyID,' isomers - ',thisName,char(10),'Iso #1 n: ',num2str(nERPs1(reOrder)),char(10),'Iso #2 n: ',num2str(nERPs2(reOrder))], 'FontSize', 8 )
        title(['Cross-fly isomers - ',plotTypes{plotI},' - (N=',num2str(size(chosenFlies,2)),')'])
        set(h,'FontSize',6);
        
        h = findobj(gca,'Type','ErrorBar');
        set(h(1),'color','r');

        %Testatory double-check
        %{
        figure
        hold on
        for isoI = 1:size(crossMean,2)
            errorbar([1:size(crossMean,1)],crossMean(reOrder,isoI),crossSEM(reOrder,isoI))
        end
        xticks([1:size(crossMean,1)])
        xticklabels(exLabels(reOrder))
        xtickangle(270)
        legend(isomerNames);
        title(['Cross-fly Isomers testatory'])
        %}

        CROSSISOMER.(plotTypes{plotI}).crossData = crossData;
        CROSSISOMER.(plotTypes{plotI}).crossMean = crossMean;
        CROSSISOMER.(plotTypes{plotI}).crossSEM = crossSEM;


        %Follow-on timepoint isomer plot (if applicable)
        if transecting %Note: Depending on plotSelector complement, this process may be done twice
            %Borrow from above single-fly method
            isoMax = struct;
            %isoIndiv = struct; %Beginnings of an extension of isomax to invdividual flies
            %a = 1;
            for timep = 1:nTimepoints %Unlike above do for *all* timepoints (Mostly to support animated plot)
                isoMax(timep).timep = timep;
                %isomTime(a).data = [R1.meanERPs(timep,:);R2.meanERPs(timep,:)];
                temp = [];
                tempMean = [];
                tempSEM = [];
                for isoI = 1:size(crossDataFull,1)
                    temp = squeeze( crossDataFull{isoI}(timep,:,:) ); %More cautious method of obtaining data used due to increased complexity of architecture
                    tempMean = [tempMean; nanmean(temp,2)'];
                    tempSEM = [tempSEM; nanstd(temp,[],2)' ./ sqrt(size(temp,2))];
                end
                %isomTime(a).mean = [nanmean(squeeze(R1.allERPs(timep,:,:)),2)';...
                %                    nanmean(squeeze(R2.allERPs(timep,:,:)),2)'];
                %isomTime(a).SEM = [nanstd(squeeze(R1.allERPs(timep,:,:)),[],2)' ./ nansum( ~isnan( squeeze(R1.allERPs(timep,:,:)) ), 2 )';...
                %                   nanstd(squeeze(R2.allERPs(timep,:,:)),[],2)' ./ nansum( ~isnan( squeeze(R2.allERPs(timep,:,:)) ), 2 )'];
                isoMax(timep).mean = tempMean; 
                isoMax(timep).SEM = tempSEM;

                %isoIndiv
                %a = a + 1;
            end
            %Cut down to useful for below plot
            isomTime = struct;
            a = 1;
            for timep = 1:timeStep:nTimepoints
                isomTime(a).timep = timep;
                isomTime(a).mean = isoMax(timep).mean;
                isomTime(a).SEM = isoMax(timep).SEM;
                a = a + 1;
            end

            %And plot
            if isoI > 2
                ['## Alert: Subsequent plot not equipped to deal with more than 2 isomers ##']
            end
            figure
            for subl = 1:size(isomTime,2)
                subplot(1,size(isomTime,2),subl)
                %plot(isomTime(subl).data(1,reOrder),'Color','b')
                errorbar( isomTime(subl).mean(1,reOrder), isomTime(subl).SEM(1,reOrder), 'Color', 'b' ) %Note: Hardcoded 2 isomers displayed
                hold on
                %plot(isomTime(subl).data(2,reOrder),'Color','r')
                errorbar( isomTime(subl).mean(2,reOrder), isomTime(subl).SEM(2,reOrder), 'Color', 'r' )
                %xticks([1:size(isomTime(subl).data,2)])
                xticklabels([])
                title(['T:',num2str(isomTime(subl).timep)])
            end
            legend({'R1','R2'})
            set(gcf,'Name',['Cross-fly IsomxTime'])
            saveas(gcf,[ resultsDirectory '/' 'Cross-fly IsomxTime.png']);


            if doAnimatedPlot

                animFrameDelay = 1/animFrameRate; %Convert to useful

                %Pre-calculate some things (for all flies) if doing extra subplots
                if extraFigSubplots %Fun fact: Even a value of 2 counts as 'true' here
                    superERPMeanIsom = cell(1,isomerCount);
                    superERPMeanSeqIsom = cell(1,isomerCount);
                    superPhotMeanIsom = cell(1,isomerCount);
                    for isoI = 1:isomerCount
                        temp = [];
                        temp2 = [];
                        a = 1;
                        for i = chosenFlies
                            temp(:,:,a) = FLIES(i).ISOMER.(isomerNames{isoI}).meanERPs; %T x Seq x Fly
                            temp2(:,:,a) = FLIES(i).ISOMER.(isomerNames{isoI}).meanPHOTs; %T x Seq x Fly; Theoretically might not exist
                            a = a + 1;
                        end
                        %superERPMean = nanmean(temp,[2,3]); %Mean
                        %superERPMeanSeq = nanmean(temp,3);
                        superERPMeanIsom{isoI} = nanmean(temp,[2,3]); %Mean
                        superERPMeanSeqIsom{isoI} = nanmean(temp,3); %Mean, sep by seq
                        superPhotMeanIsom{isoI} = nanmean(temp2,[2,3]); %Mean of phot
                    end
                end

                %Similar, for individuals

                if saveFigVid
                    vidName = [ resultsDirectory , filesep, 'isomAnim.mp4'];
                    disp(['Saving figure video to ',vidName])
                    %figxtory = getframe(hfCollab);
                    figVidWrite = VideoWriter(vidName, 'MPEG-4')
                    figVidWrite.FrameRate = animFrameRate; %Might complain if non-codec value?
                    open(figVidWrite)                    
                end

                hfCollab = figure;
                %isomerNamesSeq = {};
                if extraFigSubplots == 1
                    set(hfCollab,'units','normalized','outerposition',[0 0 0.4 1])
                %elseif extraFigSubplots == 2
                %    set(hfCollab,'units','normalized','outerposition',[0 0 1 0.6])
                end
                for t = 1:size(isoMax, 2 )
                    clf
                    %if extraFigSubplots ~= 1 %Default
                    %    clf
                    %else
                    if extraFigSubplots == 1
                        subplot(2,2,[1:2]) %Functionally same as clf, for one panel?
                    elseif extraFigSubplots == 2
                        subplot(1,2,2) %Place on right side
                    end
                    %Basal
                    hold on
                    for isoI = 1:size(isoMax(t).mean,1)
                        errorbar( isoMax(t).mean(isoI,reOrder), isoMax(t).SEM(isoI,reOrder) )
                    end
                    xticks([1:size(isoMax(t).mean,2)])
                    xticklabels(exLabels(reOrder))
                    xtickangle(270)
                    xlim([0,size(isoMax(t).mean,2)+1])
                    legend(isomerNames)
                    title('Cross-fly isomers at t = ',num2str(t))
                    %Extra
                    if extraFigSubplots == 1 || extraFigSubplots == 2
                        %Mean isomer ERP/s
                        if extraFigSubplots == 1
                            %---
                            subplot(2,2,3)
                            hold on
                            %Plot isomer means (across seq)
                            for isoI = 1:isomerCount
                                plot(superERPMeanIsom{isoI})
                            end
                            legend(isomerNames,'Location','north')
                            legend('AutoUpdate','off')
                            %Plot phot means (across seq)
                            for isoI = 1:isomerCount
                                temp = range( superERPMeanIsom{isoI} );
                                plot(normalize( superPhotMeanIsom{isoI} , 'Range', [-0.5*temp, 0.5*temp] ), 'Color',[1,0.5,0]);
                            end
                            hold off
                            %---
                        elseif extraFigSubplots == 2
                            %---
                            subplot(1,2,1)
                            hold on
                            %Pre-plot some isomers, to cheese out the legend
                            plot(superERPMeanSeqIsom{1}(:,1),'Color',isomerColours(1,:))
                            plot(superERPMeanSeqIsom{2}(:,1),'Color',isomerColours(2,:))
                            legend(isomerNames,'Location','north')
                            legend('AutoUpdate','off')
                            %Plot isomers sep by seq, coloured by isomer
                            for isoI = 1:isomerCount
                                for seq = 1:size(superERPMeanSeqIsom{isoI},2)
                                    plot(superERPMeanSeqIsom{isoI}(:,seq),'Color',isomerColours(isoI,:))
                                    %alpha(0.15)
                                end
                            end
                            %Plot phot means (across seq)
                            for isoI = 1:isomerCount
                                temp = range( superERPMeanIsom{isoI} );
                                plot(normalize( superPhotMeanIsom{isoI} , 'Range', [-0.5*temp, 0.5*temp] ), 'Color',[1,0.5,0]);
                            end
                            hold off
                            %---
                        end

                        yLim = get(gca,'YLim'); %Might be slow?
                        line([isoMax(t).timep,isoMax(t).timep],[yLim(1),yLim(2)],'LineWidth',2,'LineStyle',':','Color','g')
                        xlim([0,size(superERPMeanIsom{isoI},1)])
                        xlabel('Frame') %Note: Only correct as long as isoMax iterates at 1:1

                        %Separated sequence isomer ERPs rolling
                        if extraFigSubplots == 1
                            subplot(2,2,4)
                            hold on
                            for isoI = 1:isomerCount
                                %plot(superERPMeanIsom{isoI})
                                coords = [isoMax(t).timep-0.5*timeStep:isoMax(t).timep+0.5*timeStep];
                                coordsActual = floor(coords);
                                coordsActual( coordsActual < 1 ) = []; coordsActual( coordsActual > size(superERPMeanSeqIsom{isoI},1) ) = [];
                                for seq = 1:size(superERPMeanSeqIsom{isoI},2)
                                    %plot(superERPMeanSeqIsom{isoI}(:,seq),'Color',isomerColours(isoI,:))
                                    plot(coordsActual,superERPMeanSeqIsom{isoI}(coordsActual,seq),'Color',isomerColours(isoI,:))
                                end
                            end
                            xlim([coords(1),coords(end)])
                            yLim = get(gca,'YLim');
                            line([isoMax(t).timep,isoMax(t).timep],[yLim(1),yLim(2)],'LineWidth',2,'LineStyle',':','Color','g')
                            hold off
                        end

                    end
                    drawnow
                    hold off
                    %pause(0.1)
                    if saveFigVid
                        figxtory = getframe(hfCollab);
                        writeVideo(figVidWrite,figxtory); %Might be faster to save frames and do this outside, but more memory use
                        clear figxtory
                    end 
                    if ~extraFigSubplots %If extra stuff is being plotted will probably be pretty slow anyway
                        pause(animFrameDelay)
                    end
                    %clf
                end

                if saveFigVid == 1
                    close(figVidWrite)
                    disp(['-- Video saved successfully --'])
                end


            end



        
        end


    end  


else
    disp(['(Not performing cross-fly isomer analyses because only one fly analysed)'])
end

end