function plotSeparateERPs(allERPs,transProbAncillary)

if isempty(transProbAncillary) %Use as a proxy to infer no transition probability type
    % calculate ERP for stimuli 0 and 1
    %allERPs0 = allERPs(:,1:2:31,:); allERPs0 = sum(allERPs0,2); %Hardcoded 5back
    %allERPs1 = allERPs(:,2:2:32,:); allERPs1 = sum(allERPs1,2);
    allERPs0 = allERPs(:,1:2:size(allERPs,2)-1,:); allERPs0 = nansum(allERPs0,2); %Note: Not entirely sure data interleaved at >5back
    allERPs1 = allERPs(:,2:2:size(allERPs,2),:); allERPs1 = nansum(allERPs1,2);
    
    allERPs0(allERPs0 == 0) = nan; allERPs1(allERPs1 == 0) = nan;

    %Shift all above zero, if requested (And data actually dips below zero)
    %if zeroShiftMode == 1 && nanmin(allERPs,[],'all') < 0
    %    allERPs = allERPs - nanmin(allERPs,[],'all');    
    %    disp(['Zero shift applied'])
    %end
    
    % SEM for each stimulus
    nERPs0 = sum(~isnan(allERPs0(1,1,:)), 3);
    sdERPs0 = std(allERPs0,[],3,'omitnan');
    semERPs0 = sdERPs0 ./ sqrt(nERPs0);% broadcasting
    
    nERPs1 = sum(~isnan(allERPs1(1,1,:)), 3);
    sdERPs1 = std(allERPs1,[],3,'omitnan');
    semERPs1 = sdERPs1 ./ sqrt(nERPs1);% broadcasting
    
    % plot mean ERP for both stimuli
    figure; 
    h1 = plot(mean(allERPs0,3,'omitnan'),'b'); hold on; 
    h2 = plot(mean(allERPs1,3,'omitnan'),'r');
    
    plot(mean(allERPs0,3,'omitnan')+semERPs0,'b:');
    plot(mean(allERPs0,3,'omitnan')-semERPs0,'b:');
    
    plot(mean(allERPs1,3,'omitnan')+semERPs1,'r:');
    plot(mean(allERPs1,3,'omitnan')-semERPs1,'r:');
    
    legend([h1 h2],{'Stimulus 1','Stimulus 2'});

else %Transition probability

    nStimuli = transProbAncillary.nStimuli;
    nBack = transProbAncillary.nBackActual;

    colMap = jet(nStimuli);

    %Do everything online in one figure
        %cowabunga
    figure
    hold on
    legList = {};
    for stimi = 1:nStimuli

        thisERPs = allERPs(:, (1:nStimuli) + (stimi-1)*nStimuli ,:); %Should acquire 1:4, 5:8, etc, based on canonical order
        thisERPs = nansum(thisERPs,2);
        thisERPs(thisERPs == 0) = NaN;

        nERPs = sum(~isnan(thisERPs(1,1,:)), 3);
        sdERPs = std(thisERPs,[],3,'omitnan');
        semERPs = sdERPs ./ sqrt(nERPs);

        %Plot
        plot(mean(thisERPs,3,'omitnan'),'Color',colMap(stimi,:));
        plot(mean(thisERPs,3,'omitnan')+semERPs,'Color',colMap(stimi,:),'LineStyle',':');
        plot(mean(thisERPs,3,'omitnan')-semERPs,'Color',colMap(stimi,:),'LineStyle',':');
        
        legList = [legList,{['Stim. ',num2str(stimi),' mean']},{['Stim. ',num2str(stimi),' SEM']},{['Stim. ',num2str(stimi),' SEM']}];

    end
    legend(legList)


end