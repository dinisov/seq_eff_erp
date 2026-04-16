function plotFlies(data,chosenFlies,plotSelector,resultsDirectory, reOrder, n_back, scramLevel)

seTypes = {'amplitude','positiveAmplitude','negativeAmplitude','latencyToPeak','latencyToTrough','transect','avTransectWindow'};

directoryNames = {'Amplitude','Positive amplitude','Negative amplitude','Latency','Latency','Transect','AvTransectWindow'};

if length(data) == 1

    for p = find(plotSelector)
        
        if isfield(data,'FIT')
            modelFit = data.FIT.model_fit_amplitude; 
            modelScores = data.FIT.fit_params(1:3); 
        else
            modelFit = []; %This is actually mostly redundant now, since fit only applied if valid (i.e. Requested and n_back)
            modelScores = []; 
        end

        if p == 1 && isfield(data,'FIT')
            figure('Name',seTypes{p},'NumberTitle','off');
            %create_seq_eff_plot(data.PROFILE.(seTypes{p}).',modelFit,'errors',data.ERROR.(seTypes{p}).','scores',modelScores);
            create_seq_eff_plot(data.PROFILE.(seTypes{p}).',modelFit(reOrder),'errors',data.ERROR.(seTypes{p}).','scores',modelScores,...
                'reOrder',reOrder,'n_back',n_back,'histlength',n_back-1);
        else
            figure('Name',seTypes{p},'NumberTitle','off');
            %create_seq_eff_plot(data.PROFILE.(seTypes{p}).',[],'errors',data.ERROR.(seTypes{p}).');
            create_seq_eff_plot(data.PROFILE.(seTypes{p}).',[],'errors',data.ERROR.(seTypes{p}).',...
                'reOrder',reOrder,'n_back',n_back,'histlength',n_back-1);
        end
        if scramLevel == 1
            title(['RAW SEQUENCE SCRAMBLED'])
        elseif scramLevel == 2
            title(['SEQUENCES SCRAMBLED'])
        end

        %Make folder if not existing
        if exist( [ resultsDirectory directoryNames{p} ] ) == 0 %i.e. Not existing as folder
            mkdir( [ resultsDirectory directoryNames{p} ] )
            disp(['Results directory for ',directoryNames{p},' had to be made'])
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
            
            if p == 1 && isfield(data,'FIT')
                figure('Name',[seTypes{p} '_fly_' num2str(fly)],'NumberTitle','off');
                %create_seq_eff_plot(data(fly).PROFILE.(seTypes{p}).',modelFit,'errors',data(fly).ERROR.(seTypes{p}).','scores',modelScores);
                create_seq_eff_plot(data(fly).PROFILE.(seTypes{p}).',modelFit,'errors',data(fly).ERROR.(seTypes{p}).','scores',modelScores,...
                    'reOrder',reOrder,'n_back',n_back,'histlength',n_back-1);
            else
                figure('Name',[seTypes{p} '_fly_' num2str(fly)],'NumberTitle','off');
                %create_seq_eff_plot(data(fly).PROFILE.(seTypes{p}).',[],'errors',data(fly).ERROR.(seTypes{p}).');
                create_seq_eff_plot(data(fly).PROFILE.(seTypes{p}).',[],'errors',data(fly).ERROR.(seTypes{p}).',...
                    'reOrder',reOrder,'n_back',n_back,'histlength',n_back-1);
            end

            %Make folder if not existing
            if exist( [ resultsDirectory directoryNames{p} ] ) == 0 %i.e. Not existing as folder
                mkdir( [ resultsDirectory directoryNames{p} ] )
                disp(['Results directory for ',directoryNames{p},' had to be made'])
            end

            saveas(gcf,[ resultsDirectory directoryNames{p} '/' seTypes{p} '_fly' num2str(fly) '.png']);

        end
    
    end

end