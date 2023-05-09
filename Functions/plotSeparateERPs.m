function plotSeparateERPs(allERPs)

    % calculate ERP for stimuli 0 and 1
    allERPs0 = allERPs(:,1:2:31,:); allERPs0 = sum(allERPs0,2);
    allERPs1 = allERPs(:,2:2:32,:); allERPs1 = sum(allERPs1,2);
    
    allERPs0(allERPs0 == 0) = nan; allERPs1(allERPs1 == 0) = nan;
    
    % SEM for each stimulus
    nERPs0 = sum(~isnan(allERPs0(1,1,:)), 3);
    sdERPs0 = std(allERPs0,[],3,'omitnan');
    semERPs0 = sdERPs0 ./ sqrt(nERPs0);% broadcasting
    
    nERPs1 = sum(~isnan(allERPs1(1,1,:)), 3);
    sdERPs1 = std(allERPs1,[],3,'omitnan');
    semERPs1 = sdERPs1 ./ sqrt(nERPs1);% broadcasting
    
    % plot mean ERP for both stimuli
    figure; 
    plot(mean(allERPs0,3,'omitnan'),'r'); hold on; 
    plot(mean(allERPs0,3,'omitnan')+semERPs0,'r:');
    plot(mean(allERPs0,3,'omitnan')-semERPs0,'r:');
    
    plot(mean(allERPs1,3,'omitnan'),'b');
    plot(mean(allERPs1,3,'omitnan')+semERPs1,'b:');
    plot(mean(allERPs1,3,'omitnan')-semERPs1,'b:');
    
    legend({'Stimulus 1','Stimulus 0'});