function [ISOMER] = plotIsomers(allERPs, isomer, window, n_back, resampleFreq, transectTime, avTransectWindow, plotSelector, reOrder, flyID, plotIndividualFlies, allPHOTs)
%plotIsomers Plot isomers
%   If isomer input is empty it calculates "natural" isomers (all ending with 0/1)
%   "isomer" defines a set of 16 sequences, the other isomer being 1-isomer

    if ~isempty(isomer)

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

        %R1 = calculateSEs(allERPs1,[],0,window, resampleFreq); %Missing arguments for current gen implementation?
        %R2 = calculateSEs(allERPs2,[],0,window, resampleFreq);
        R1 = calculateSEs(allERPs1,allPHOTs1,0,window, resampleFreq); %Support for phot even for isomers
        R2 = calculateSEs(allERPs2,allPHOTs2,0,window, resampleFreq);

        %figure; create_seq_eff_plot(R1.amplitudeSEs.',R2.amplitudeSEs.');
        figure; create_seq_eff_plot(R1.amplitudeSEs.',R2.amplitudeSEs.', 'reOrder',reOrder,'n_back',n_back);
        h = legend({'First Isomer','Second Isomer'});
        title([flyID,' isomers'])
        set(h,'FontSize',6);

        check isomer architecture wrt below usage
        ISOMER = isomer;

        %----------------------------------------------------------------------------------------------
    
    else
        
        index1 = zeros(1,(2^n_back)); index2 = zeros(1,(2^n_back));
        
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
        index1 = index1+fliplr(index1); index1((size(index1,2)/2)+1:end) = []; index1 = index1(seq_eff_order(n_back)); %Will crash if 1/2 not whole, but that is fine
        index2 = index2+fliplr(index2); index2((size(index2,2)/2)+1:end) = []; index2 = index2(seq_eff_order(n_back));
       
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

        nERPs1 = nansum( ~isnan(allERPs1(1,:,:)) , 3 );
        nERPs2 = nansum( ~isnan(allERPs2(1,:,:)) , 3 );

        %R1 = calculateSEs(allERPs1,[],0,window, resampleFreq);
        %R2 = calculateSEs(allERPs2,[],0,window, resampleFreq);
        R1 = calculateSEs(allERPs1,allPHOTs1,0,window, resampleFreq, transectTime, avTransectWindow, plotSelector, n_back);
        R2 = calculateSEs(allERPs2,allPHOTs2,0,window, resampleFreq, transectTime, avTransectWindow, plotSelector, n_back);

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
                
                h = findobj(gca,'Type','ErrorBar');
                set(h(1),'color','r');
                
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
        
    end

end

