function plotFlies(data,chosenFlies,plotSelector,resultsDirectory, reOrder, n_back, scramLevel)

seTypes = {'amplitude','positiveAmplitude','negativeAmplitude','latencyToPeak','latencyToTrough','transect','avTransectWindow'};

directoryNames = {'Amplitude','Positive amplitude','Negative amplitude','Latency','Latency','Transect','AvTransectWindow'};

addTitle = 1; %Whether to add a little title above the plot relating to the fly and block number

if length(data) == 1 %Note: May not actually be reachable nowadays (Due to data being length maximum fly number)

    for p = find(plotSelector)
        
        if isfield(data,'FIT')
            modelFit = data.FIT.model_fit_amplitude; 
            modelScores = data.FIT.fit_params(1:3); 
        else
            modelFit = []; %This is actually mostly redundant now, since fit only applied if valid (i.e. Requested and n_back)
            modelScores = []; 
        end

        if ~data.transProbDesign %Normal

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

        else %Transition probabilities [Do not support model fitting]

            figure('Name',seTypes{p},'NumberTitle','off');
            create_seq_eff_plot(data.PROFILE.(seTypes{p}).',[],'errors',data.ERROR.(seTypes{p}).',...
                'reOrder', reOrder,'n_back',-data.transProbAncillary.nStimuli,'histlength',data.transProbAncillary.nBackActual,...
                    'nStimuli',data.transProbAncillary.nStimuli); 
                check parameter correctness

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
            
            if ~data(fly).transProbDesign

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

            else

                %(Remember: No models)
                figure('Name',[seTypes{p} '_fly_' num2str(fly)],'NumberTitle','off');
                create_seq_eff_plot(data(fly).PROFILE.(seTypes{p}).',[],'errors',data(fly).ERROR.(seTypes{p}).',...
                    'reOrder', data(fly).transProbAncillary.reOrderActual,...
                    'n_back',-data(fly).transProbAncillary.nStimuli,...
                    'histlength',data(fly).transProbAncillary.nBackActual,... %Unlike in plotIsomers, no -1 here, because not isomers
                    'nStimuli',data(fly).transProbAncillary.nStimuli, ...
                    'overrideLabels',data(fly).transProbAncillary.transLabels);
                    %Reminder: For this type of analysis, n_back is purely to denote trans prob nature, histlength should be 'actual' nBack (-> hl^2 e.g. 2^2),
                    %and nStimuli is nStimuli, used only for label collection?
                    %As set currently, by feeding nStimuli twice, we prepare correctly for the full range of transition combinations, with no stacking
                        %(Graphs relating to stacking are covered by plotIsomers)
                if addTitle
                    title(['Fly ',num2str(fly),' - Block ',data(fly).block,' - ',seTypes{p}])
                end

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