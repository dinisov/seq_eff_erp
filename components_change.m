close all; clear;

load results

struct_names = {'freq1dot25hz','freq6dot25hz','freq12dot5hz','freq25hz'};

scores = struct;

%concatenate data
for f = 1:4
   
    data = results.(struct_names{f});
    
    slrp = []; lrpr = []; weird = [];
    
    for j = 1: length(data)
       
        slrp = [slrp abs(data(j).fit_params(1))]; %#ok<AGROW>
        lrpr = [lrpr abs(data(j).fit_params(2))]; %#ok<AGROW>
        weird = [weird abs(data(j).fit_params(3))]; %#ok<AGROW>
        
    end
    
    scores(f).slrp = slrp;
    scores(f).lrpr = lrpr;
    scores(f).weird = weird;
    
end

mean_slrp = zeros(1,4);
sem_slrp = zeros(1,4);
mean_lrpr = zeros(1,4);
sem_lrpr = zeros(1,4);
mean_weird = zeros(1,4);
sem_weird = zeros(1,4);

for f = 1:4
    
    mean_slrp(f) = mean(scores(f).slrp);
    sem_slrp(f) = std(scores(f).slrp)/sqrt(length(scores(f).slrp));
    mean_lrpr(f) = mean(scores(f).lrpr);
    sem_lrpr(f) = std(scores(f).lrpr)/sqrt(length(scores(f).lrpr));
    mean_weird(f) = mean(scores(f).weird);
    sem_weird(f) = std(scores(f).weird)/sqrt(length(scores(f).weird));
    
end

errorbar(mean_slrp,2*sem_slrp,'linewidth',2); hold on;
errorbar(mean_lrpr,2*sem_lrpr,'linewidth',2);
errorbar(mean_weird,2*sem_weird,'linewidth',2);
xlim([0 5]); legend({'SLRP','LRPR','WEIRD'});
set(gca,'xtick',[1 2 3 4],'xticklabel',[1.25 6.25 12.5 25]);
ylabel('abs(Component score)'); xlabel('Frequency');