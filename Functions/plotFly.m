function plotFly(data,plotSelector,resultsDirectory)

seTypes = {'amplitudeSEs','positiveAmplitudeSEs','negativeAmplitudeSEs','latencyToPeakSEs','latencyToTroughSEs'};

errorNames = {'semAmplSEs','semPosAmplSEs','semNegAmplSEs'};

for p = find(plotSelector)

    if p < 4% no errors for latencies
        figure('Name',seTypes{p},'NumberTitle','off');
        create_seq_eff_plot(data.(seTypes{p}),[],'errors',data.(errorNames{p}));
    
    else
       
        figure('Name',seTypes{p},'NumberTitle','off');
        create_seq_eff_plot(data.(seTypes{p}),[]);
        
    end

    saveas(gcf,[ resultsDirectory seTypes{p} '.png']);
    
end

end