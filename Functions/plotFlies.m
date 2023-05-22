function plotFlies(data,chosenFlies,plotSelector,resultsDirectory)

seTypes = {'amplitudeSEs','positiveAmplitudeSEs','negativeAmplitudeSEs','latencyToPeakSEs','latencyToTroughSEs'};

errorNames = {'semAmplSEs','semPosAmplSEs','semNegAmplSEs'};

directoryNames = {'Amplitude','Positive amplitude','Negative amplitude','Latency','Latency'};

if length(data) == 1

    for p = find(plotSelector)
        
        if isfield(data,'FIT')
            modelFit = data.FIT.model_fit_amplitude; 
            modelScores = data.FIT.fit_params(1:3); 
        else
            modelFit = []; 
            modelScores = []; 
        end

        if p == 1
            figure('Name',seTypes{p},'NumberTitle','off');
            create_seq_eff_plot(data.(seTypes{p}).',modelFit,'errors',data.(errorNames{p}).','scores',modelScores);
        elseif (p > 1) && (p < 4)
            figure('Name',seTypes{p},'NumberTitle','off');
            create_seq_eff_plot(data.(seTypes{p}).',[],'errors',data.(errorNames{p}).');
        else
            figure('Name',seTypes{p},'NumberTitle','off');
            create_seq_eff_plot(data.(seTypes{p}).',[]);
        end

        saveas(gcf,[ resultsDirectory seTypes{p} '.png']);

    end

else
    
    for fly = chosenFlies
        
        if isfield(data(fly),'FIT')
            modelFit = data(fly).FIT.model_fit_amplitude; 
            modelScores = data(fly).FIT.fit_params(1:3); 
        else
            modelFit = []; 
            modelScores = []; 
        end
   
        for p = find(plotSelector)
            
            if p == 1
                figure('Name',seTypes{p},'NumberTitle','off');
                create_seq_eff_plot(data(fly).(seTypes{p}).',modelFit,'errors',data(fly).(errorNames{p}).','scores',modelScores);
            elseif (p > 1) && (p < 4)
                figure('Name',seTypes{p},'NumberTitle','off');
                create_seq_eff_plot(data(fly).(seTypes{p}).',[],'errors',data(fly).(errorNames{p}).');
            else% no errors for latencies
                figure('Name',seTypes{p},'NumberTitle','off');
                create_seq_eff_plot(data(fly).(seTypes{p}).',[]);
            end

            saveas(gcf,[ resultsDirectory directoryNames{p} '/' seTypes{p} '_fly' num2str(fly) '.png']);

        end
    
    end

end