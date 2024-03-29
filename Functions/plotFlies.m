function plotFlies(data,chosenFlies,plotSelector,resultsDirectory)

seTypes = {'amplitude','positiveAmplitude','negativeAmplitude','latencyToPeak','latencyToTrough'};

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
            create_seq_eff_plot(data.PROFILE.(seTypes{p}).',modelFit,'errors',data.ERROR.(seTypes{p}).','scores',modelScores);
        else
            figure('Name',seTypes{p},'NumberTitle','off');
            create_seq_eff_plot(data.PROFILE.(seTypes{p}).',[],'errors',data.ERROR.(seTypes{p}).');
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
                figure('Name',[seTypes{p} '_fly_' num2str(fly)],'NumberTitle','off');
                create_seq_eff_plot(data(fly).PROFILE.(seTypes{p}).',modelFit,'errors',data(fly).ERROR.(seTypes{p}).','scores',modelScores);
            else
                figure('Name',[seTypes{p} '_fly_' num2str(fly)],'NumberTitle','off');
                create_seq_eff_plot(data(fly).PROFILE.(seTypes{p}).',[],'errors',data(fly).ERROR.(seTypes{p}).');
            end

            saveas(gcf,[ resultsDirectory directoryNames{p} '/' seTypes{p} '_fly' num2str(fly) '.png']);

        end
    
    end

end