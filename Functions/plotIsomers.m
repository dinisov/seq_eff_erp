function [ISOMER] = plotIsomers(allERPs, isomer, window, n_back, resampleFreq,...
    transectTime, avTransectWindow, plotSelector, reOrder, flyID, plotIndividualFlies, allPHOTs, allTIMEs,...
    transProbDesign, transProbAncillary)
%plotIsomers Plot isomers
%   If isomer input is empty it calculates "natural" isomers (all ending with 0/1)
%   "isomer" defines a set of 16 sequences, the other isomer being 1-isomer

    if ~isempty(isomer) %Not normal case?

        if n_back ~= 5
            ['## Alert: Isomer case not coded for non-5 n_back ##']
            crash = yes
            %Too lazy to code right now; Will do if this case actually occurs
        end

        isomer1 = isomer; isomer2 = 1-isomer;

        index1 = zeros(1,(2^n_back)); index2 = zeros(1,(2^n_back));

        %convert to index
        for n = 5:length(isomer)

           idx1 = bin2dec(num2str(isomer1(n-n_back+1:n))) + 1;
           idx2 = bin2dec(num2str(isomer2(n-n_back+1:n))) + 1;

           index1(idx1) = idx1;
           index2(idx2) = idx2;

        end

        index1 = index1+fliplr(index1); index1(17:end) = []; index1 = index1(seq_eff_order(n_back));
        index2 = index2+fliplr(index2); index2(17:end) = []; index2 = index2(seq_eff_order(n_back));

        %uncomment to check if your isomer is sound (should be 1:32)
    % 	disp(sort(union(index1,index2)));

        %grab all the ERPs for each isomer, already in SE standard order
        allERPs1 = allERPs(:,index1,:); 
        allERPs2 = allERPs(:,index2,:);
        %New addition: Ditto for phot
        if ~isempty(allPHOTs)
            allPHOTs1 = allPHOTs(:,index1,:); %Note: Assumes same architecture as allERPs
            allPHOTs2 = allPHOTs(:,index2,:); %Also, hardcoded only two isomers?
        else
            allPHOTs1 = [];
            allPHOTs2 = [];
        end
        %New addition: And for times
        if ~isempty(allTIMEs)
            allTIMEs1 = allTIMEs(:,index1,:); %Note: Assumes same architecture as allERPs
            allTIMEs2 = allTIMEs(:,index2,:); %Also, hardcoded only two isomers?
        else
            allTIMEs1 = [];
            allTIMEs2 = [];
        end

        %R1 = calculateSEs(allERPs1,[],0,window, resampleFreq); %Missing arguments for current gen implementation?
        %R2 = calculateSEs(allERPs2,[],0,window, resampleFreq);
        R1 = calculateSEs(allERPs1,allPHOTs1,allTIMEs1,0,window, resampleFreq); %Support for phot even for isomers
        R2 = calculateSEs(allERPs2,allPHOTs2,allTIMEs2,0,window, resampleFreq);

        %figure; create_seq_eff_plot(R1.amplitudeSEs.',R2.amplitudeSEs.');
        figure; create_seq_eff_plot(R1.amplitudeSEs.',R2.amplitudeSEs.', 'reOrder',reOrder,'n_back',n_back);
        h = legend({'First Isomer','Second Isomer'});
        title([flyID,' isomers'])
        set(h,'FontSize',6);

        check isomer architecture wrt below usage
        ISOMER = isomer;

        %----------------------------------------------------------------------------------------------
    
    else %Nowadays normal case?
        
        if ~transProbDesign %Normal

            index1 = zeros(1,(2^n_back)); 
            index2 = zeros(1,(2^n_back));
            
            %Hardcoded 5-back
            %{
            for i = 1:2:31
                index1(i) = i; 
                index2(i+1) = i+1;
            end
            %}
            %Dynamic
            for i = 1:2:(2^n_back)-1
                index1(i) = i; 
                index2(i+1) = i+1;
            end
            
            %Hardcode
            %index1 = index1+fliplr(index1); index1(17:end) = []; index1 = index1(seq_eff_order(n_back));
            %index2 = index2+fliplr(index2); index2(17:end) = []; index2 = index2(seq_eff_order(n_back));
            %Dynamic
            index1 = index1+fliplr(index1); 
            index1((size(index1,2)/2)+1:end) = []; 
            index1 = index1(seq_eff_order(n_back)); %Will crash if 1/2 not whole, but that is fine
            index2 = index2+fliplr(index2); 
            index2((size(index2,2)/2)+1:end) = []; 
            index2 = index2(seq_eff_order(n_back));
       
            %disp(sort(union(index1,index2)));
            nunion = nanmax(union(index1,index2));
            disp(['n-back: ',num2str(n_back),', total number of unique sequences: ',num2str(nunion),' (',num2str(0.5*nunion),' functional)'])
            %disp(sort(union(index1,index2)));
            
            %grab all the ERPs for each isomer, already in SE standard order
            allERPs1 = allERPs(:,index1,:); 
            allERPs2 = allERPs(:,index2,:);
            %New addition: Ditto for phot
            if ~isempty(allPHOTs)
                allPHOTs1 = allPHOTs(:,index1,:); %Note: Assumes same architecture as allERPs
                allPHOTs2 = allPHOTs(:,index2,:); %Also, hardcoded only two isomers?
            else
                allPHOTs1 = [];
                allPHOTs2 = [];
            end
            %And for times
            if ~isempty(allTIMEs)
                allTIMEs1 = allTIMEs(:,index1,:); %Note: Assumes same architecture as allERPs
                allTIMEs2 = allTIMEs(:,index2,:); %Also, hardcoded only two isomers?
            else
                allTIMEs1 = [];
                allTIMEs2 = [];
            end
    
            nERPs1 = nansum( ~isnan(allERPs1(1,:,:)) , 3 );
            nERPs2 = nansum( ~isnan(allERPs2(1,:,:)) , 3 );
    
            %R1 = calculateSEs(allERPs1,[],0,window, resampleFreq);
            %R2 = calculateSEs(allERPs2,[],0,window, resampleFreq);
            R1 = calculateSEs(allERPs1,allPHOTs1,allTIMEs1,0,window, resampleFreq, transectTime, avTransectWindow, plotSelector, n_back, [], []);
            R2 = calculateSEs(allERPs2,allPHOTs2,allTIMEs2,0,window, resampleFreq, transectTime, avTransectWindow, plotSelector, n_back, [], []);
    
            %Amplitude isomers (Hardcoded)
            %{
            %figure; create_seq_eff_plot([R1.PROFILE.amplitude.' R2.PROFILE.amplitude.'],[],'errors',[R1.ERROR.amplitude.' R2.ERROR.amplitude.']);
            figure; create_seq_eff_plot([R1.PROFILE.amplitude.' R2.PROFILE.amplitude.'],[],'errors',[R1.ERROR.amplitude.' R2.ERROR.amplitude.'],...
                'reOrder', reOrder,'n_back',n_back,'histlength',n_back-1);
            h = legend({'First Isomer','Second Isomer'});
            title([flyID,' isomers'])
            set(h,'FontSize',6);
            
            h = findobj(gca,'Type','ErrorBar');
            set(h(1),'color','r');
            %}
            %Generalised isomer plotting
            %whichPlots = find(plotSelector == 1);
            if plotIndividualFlies
                for i = find(plotSelector == 1)
                    if i == 1
                        thisData1 = R1.PROFILE.amplitude;
                        thisData2 = R2.PROFILE.amplitude;
                        thisError1 = R1.ERROR.amplitude;
                        thisError2 = R2.ERROR.amplitude;
                        thisName = 'amplitude';
                    elseif i == 2
                        thisData1 = R1.PROFILE.positiveAmplitude;
                        thisData2 = R2.PROFILE.positiveAmplitude;
                        thisError1 = R1.ERROR.positiveAmplitude;
                        thisError2 = R2.ERROR.positiveAmplitude;
                        thisName = 'positiveAmplitude';
                    elseif i == 3
                        thisData1 = R1.PROFILE.negativeAmplitude;
                        thisData2 = R2.PROFILE.negativeAmplitude;
                        thisError1 = R1.ERROR.negativeAmplitude;
                        thisError2 = R2.ERROR.negativeAmplitude;
                        thisName = 'negativeAmplitude';
                    elseif i == 4
                        thisData1 = R1.PROFILE.latencyToPeak;
                        thisData2 = R2.PROFILE.latencyToPeak;
                        thisError1 = R1.ERROR.latencyToPeak;
                        thisError2 = R2.ERROR.latencyToPeak;
                        thisName = 'latencyToPeak';
                    elseif i == 5
                        thisData1 = R1.PROFILE.latencyToTrough;
                        thisData2 = R2.PROFILE.latencyToTrough;
                        thisError1 = R1.ERROR.latencyToTrough;
                        thisError2 = R2.ERROR.latencyToTrough;
                        thisName = 'latencyToTrough';
                    elseif i == 6
                        thisData1 = R1.PROFILE.transect;
                        thisData2 = R2.PROFILE.transect;
                        thisError1 = R1.ERROR.transect;
                        thisError2 = R2.ERROR.transect;
                        thisName = 'transect';
                    elseif i == 7
                        thisData1 = R1.PROFILE.avTransectWindow;
                        thisData2 = R2.PROFILE.avTransectWindow;
                        thisError1 = R1.ERROR.avTransectWindow;
                        thisError2 = R2.ERROR.avTransectWindow;
                        thisName = 'avTransectWindow';
                    else
                        ['## Plot case unknown to plotIsomers ##']
                        crash = yes
                    end
                    figure; create_seq_eff_plot([thisData1.' thisData2.'],[],'errors',[thisError1.' thisError2.'],...
                    'reOrder', reOrder,'n_back',n_back,'histlength',n_back-1);
                    h = legend({'First Isomer','Second Isomer'});
                    %title([flyID,' isomers'])
                    %title([flyID,' isomers - ',thisName])
                    title([flyID,' isomers - ',thisName,char(10),'Iso #1 n: ',num2str(nERPs1(reOrder)),char(10),'Iso #2 n: ',num2str(nERPs2(reOrder))], 'FontSize', 8 )
                    set(h,'FontSize',6);
                    
                    try
                        h = findobj(gca,'Type','ErrorBar');
                        set(h(1),'color','r');
                    catch
                        ['-# Error bar find failure #-']
                    end
                    
                end
        
                %Create plot of n
                if exist('nERPs1') %Only do if calculated
                    figure
                    create_seq_eff_plot([nERPs1.' nERPs2.'],[],...
                        'reOrder', reOrder,'n_back',n_back,'histlength',n_back-1);
                    h = legend({'First Isomer','Second Isomer'});
                    ylabel('Event count')
                    title([flyID,' isomer n plot',char(10),'Iso #1 n: ',num2str(nERPs1(reOrder)),char(10),'Iso #2 n: ',num2str(nERPs2(reOrder))], 'FontSize', 8 )
                    set(h,'FontSize',6);
                end
            end
    
            %Save isomer data
            ISOMER = struct;
            ISOMER.R1 = R1;
            ISOMER.R2 = R2;

        else %Transition probabilities

            nBackActual = transProbAncillary.nBackActual;
            nStimuli = transProbAncillary.nStimuli;

            indexes = {};
            for stimi = 1:nStimuli
                indexes{stimi} = [(1:nStimuli) + (stimi-1)*nStimuli]; %Up for debate, but currently we define a transition probability 'isomer' as 'ends on a particular stim'
            end

            disp([num2str(size(indexes,2)),' transition probability isomer indices prepared'])

            allERPsIsomer = {};
            nERPsIsomer = {};
            %if ~isempty(allPHOTs)
            allPHOTsIsomer = {};
            %end
            allTIMEsIsomer = {};
            for stimi = 1:size(indexes,2)
                allERPsIsomer{stimi} = allERPs(:,indexes{stimi},:);
                nERPsIsomer{stimi} = nansum( ~isnan(allERPsIsomer{stimi}(1,:,:)) , 3 );
                if ~isempty(allPHOTs)
                    allPHOTsIsomer{stimi} = allPHOTs(:,indexes{stimi},:);
                end
                if ~isempty(allTIMEs)
                    allTIMEsIsomer{stimi} = allTIMEs(:,indexes{stimi},:);
                end
            end

            Rs = {};
            for stimi = 1:size(indexes,2)
                Rs{stimi} = calculateSEs( allERPsIsomer{stimi},allPHOTsIsomer{stimi},allTIMEsIsomer{stimi},0,window, resampleFreq, transectTime, avTransectWindow, plotSelector, ...
                    [], nStimuli, [] );
                    %A new argument is added to calculateSEs that allows for an overridden nBack, which is useful here
                        %Note that we don't supply custom labels, which means we can't do aux plots here
            end
                %Reminder: Each cell of Rs is an isomer (equivalent to R1, R2, etc), which here are defined as "Transition ending in Stim #<stimi>"
                    %Columns within PROFILE etc are individual 'sequences' (e.g. 1->1, 2->1, 3->1, 4->1, etc)

            %Now plot (maybe)
            metricIndex = {'amplitude','positiveAmplitude','negativeAmplitude','latencyToPeak','latencyToTrough','transect','avTransectWindow'}; %I am lazy
            if plotIndividualFlies
                for i = find(plotSelector == 1)
                    thisMetric = metricIndex{i};
                    thisData = [];
                    thisError = [];
                    legList = {};
                    for stimi = 1:nStimuli
                        thisData = [thisData, Rs{stimi}.PROFILE.(thisMetric)' ];
                        thisError = [thisError, Rs{stimi}.ERROR.(thisMetric)' ];
                        %These lines subserve both the hardcoded R1/R2 specification above, and the collation previously only performed at the moment of plot creation
                        legList = [legList,{['Isomer ',num2str(stimi)]}];
                    end
                    %Rows of thisData are 'sequences' within each isomer, Columns are isomers
                        %(Because ')

                    figure; 
                    create_seq_eff_plot([thisData],[],'errors',[thisError],...
                    'reOrder', reOrder,'n_back',-nStimuli,'histlength',nBackActual-1,...
                    'nStimuli',nStimuli, 'overrideLabels',transProbAncillary.histLabels); %Note that here histlength is -1 because isomers
                        %Note: The indicisation of stimuli and sequences has been empirically confirmed, but not rigorously tested
                    %h = legend({'First Isomer','Second Isomer'});
                    legList = {};
                    %nERPs = [];
                    for stimi = 1:nStimuli
                        legList = [legList,{['Stim. n',num2str(stimi)]}];
                        %nERPs = [nERPs;nERPsIsomer{stimi}];
                    end
                    h = legend(legList);
                    legend('AutoUpdate','off')
                    titleStr = [flyID,' isomers - ',thisMetric];
                    for stimi = 1:nStimuli
                        titleStr  = [titleStr,char(10),'n',num2str(stimi),': ',num2str(nERPsIsomer{stimi})]; %Not 100% confident about this indicisation
                            %Reminder that columns of nERPsIsomer are *probably* sequence (and rows are the stim ending)
                    end
                    title(titleStr, 'FontSize', 8 )
                    set(h,'FontSize',6);
                    
                    %h = findobj(gca,'Type','ErrorBar');
                    %set(h(1),'color','r'); %This just makes things confusing for multiple stims

                    
                end
            end

            ISOMER = struct;
            for stimi = 1:nStimuli
                ISOMER.(['R',num2str(stimi)]) = Rs{stimi};
            end

        end %transProbDesign end
        
    end %isempty end

end %function end

