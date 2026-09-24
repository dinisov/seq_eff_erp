function blocks = calculatePeaks(blocks, aux_plots, options)
%calculatePeaks Calculate photodiode peaks
%   Detailed explanation goes here
arguments
    blocks struct
    aux_plots double
    options.photMovMaxWindow double = [20,20]
    options.n_back double = 5; %Adds support for negative n_back (Transition probabilities)
    options.transMergeWindow double = 0.1; %Maximum ISI size that potentially simultaneous events are allowed to be
end
photMovMaxWindow = options.photMovMaxWindow;
n_back = options.n_back;
transMergeWindow = options.transMergeWindow;

%Specify index for transition probabilities, if relevant
if n_back <= 0
    photTransIndex = [1,0,0;...
                      0,1,0;...
                      1,1,0;...
                      0,0,1]; %This *needs* to be a 1:1 match with the extraDots specification in the PTB script (Presumably for photMode 2)
end


for b = 1:length(blocks)
    
    peakThreshold = blocks(b).peakThreshold; %Note: Implicitly this is the threshold for PHOT2
    if isfield(blocks,'peakThresholdOne') %Acquire 'new' threshold for PHOT1
        peakThresholdOne = blocks(b).peakThresholdOne;
    else
        peakThresholdOne = .1; %Old default value
    end
    ISI = blocks(b).ISI;
    resampleFreq = blocks(b).resampleFreq;
    
    % one photodiode (different size peaks)
    if n_back > 0 %Normal
        if blocks(b).PHOTType == 1
    
            PHOT = -blocks(b).PHOT(3,:)/max(blocks(b).PHOT(3,:));
            %QA
            if ~any(~isnan(PHOT)) %Boolean may be assembled incorrectly
                ['## No non-NaN PHOT data found! ##']
                crash = yes
            end
    
            if aux_plots
                figure; plot(PHOT);
            end
            %PHOT = movmax(PHOT,[20 20]);
            PHOT = movmax(PHOT,photMovMaxWindow);
            if aux_plots    
                hold on; plot(PHOT);
            end
    
            blocks(b).PHOT = [PHOT; PHOT; PHOT];
    
            % find peaks
            %[PKS_PHOT1,LOCS_PHOT1] = findpeaksbase(PHOT, 'MinPeakHeight' , .1 , 'MinPeakDistance' , 1/2*ISI*resampleFreq );
            %[PKS_PHOT2,LOCS_PHOT2] = findpeaksbase(PHOT , 'MinPeakHeight' , peakThreshold , 'MinPeakDistance' , 1/2*ISI*resampleFreq ); 
            [PKS_PHOT1,LOCS_PHOT1] = findpeaksbase(PHOT, 'MinPeakHeight' , peakThresholdOne , 'MinPeakDistance' , 1/2*ISI*resampleFreq );
            [PKS_PHOT2,LOCS_PHOT2] = findpeaksbase(PHOT , 'MinPeakHeight' , peakThreshold , 'MinPeakDistance' , 1/2*ISI*resampleFreq ); 
    
            [LOCS_PHOT1, ind_locs_phot1] = setdiff(LOCS_PHOT1, LOCS_PHOT2);
            PKS_PHOT1 = PKS_PHOT1(ind_locs_phot1);
    
            blocks(b).PKS_PHOT1 = PKS_PHOT1;
            blocks(b).PKS_PHOT2 = PKS_PHOT2;
            
        % two photodiodes
        elseif blocks(b).PHOTType == 2
            %blocks(b)
           
            PHOT1 = blocks(b).PHOT(1,:)/max(blocks(b).PHOT(1,:));
            if blocks(b).dataIsMulti == 0
            PHOT2 = -blocks(b).PHOT(2,:)/max(blocks(b).PHOT(2,:)); %make -blocks if Melvyn rig
            else
            PHOT2 = blocks(b).PHOT(2,:)/max(blocks(b).PHOT(2,:)); %No inversion for multichannel data (Currently)
            end
    
            %QA
            if ~any(~isnan(PHOT1)) || ~any(~isnan(PHOT2))
                ['## No PHOT1 or PHOT2 data found! ##']
                crash = yes
            end
    
            if aux_plots
                photty1 = figure;
                plot(PHOT1)
                title('PHOT1')
                photty2 = figure;
                plot(PHOT2)
                title('PHOT2')
            end
    
            %PHOT1 = movmax(PHOT1,[40 40]);%was [40 40]
            %PHOT2 = movmax(PHOT2,[40 40]);%was [40 40]
            PHOT1 = movmax(PHOT1,photMovMaxWindow);%was [40 40]
            PHOT2 = movmax(PHOT2,photMovMaxWindow);%was [40 40]
            
            blocks(b).PHOT(1,:) = PHOT1;
            blocks(b).PHOT(2,:) = PHOT2;
    
            %find beginning of stimuli
            LOCS_PHOT1 = find(diff(PHOT1 > peakThreshold) > 0) + 1;
            LOCS_PHOT2 = find(diff(PHOT2 > peakThreshold) > 0) + 1;
    
            %first round of double peak detection
            badLOCS_PHOT1 = LOCS_PHOT1([false diff(LOCS_PHOT1) < (0.8*ISI*resampleFreq)]); %Was 0.8
            badLOCS_PHOT2 = LOCS_PHOT2([false diff(LOCS_PHOT2) < (0.8*ISI*resampleFreq)]);
    
            if aux_plots
                %figure; hold on; plot(PHOT1); scatter(LOCS_PHOT1,zeros(size(LOCS_PHOT1)),'b','filled'); scatter(badLOCS_PHOT1,zeros(size(badLOCS_PHOT1)),'m','filled'); title(['Phot 1 data (Max:',num2str(max(blocks(b).PHOT(1,:))),')'])
                %figure; hold on; plot(PHOT2); scatter(LOCS_PHOT2,zeros(size(LOCS_PHOT2)),'r','filled'); scatter(badLOCS_PHOT2,zeros(size(badLOCS_PHOT2)),'m','filled'); title(['Phot 2 data (Max:',num2str(max(blocks(b).PHOT(2,:))),')'])
                figure(photty1); hold on; plot(PHOT1); scatter(LOCS_PHOT1,zeros(size(LOCS_PHOT1)),'b','filled'); scatter(badLOCS_PHOT1,zeros(size(badLOCS_PHOT1)),'m','filled'); title(['Phot 1 data (Max:',num2str(max(blocks(b).PHOT(1,:))),')'])
                figure(photty2); hold on; plot(PHOT2); scatter(LOCS_PHOT2,zeros(size(LOCS_PHOT2)),'r','filled'); scatter(badLOCS_PHOT2,zeros(size(badLOCS_PHOT2)),'m','filled'); title(['Phot 2 data (Max:',num2str(max(blocks(b).PHOT(2,:))),')'])
            end
            
        end

        % fuse locations of PHOT1 and PHOT2 (I figured this was quicker than concatenating and sorting)
        LOCS = zeros(1,length(blocks(b).PHOT)); 
        LOCS(LOCS_PHOT1) = LOCS_PHOT1; LOCS(LOCS_PHOT2) = LOCS_PHOT2;
        LOCS = LOCS(logical(LOCS));
        
        blocks(b).LOCS = LOCS; %This is used in sortSEs to creat ERP windows
        blocks(b).LOCS_PHOT1 = LOCS_PHOT1; %This is used in inferRandomSequence to infer the sequence
        blocks(b).LOCS_PHOT2 = LOCS_PHOT2;
        blocks(b).transProbDesign = 0; %New variable to simplify tracking of these experiments

    elseif n_back <= 0 %Transition probabilities
        disp(['Processing phot for transition probabilities of ',num2str(abs(n_back)),' stimuli'])

        %(Borrowed mostly from above)
        PHOT1 = blocks(b).PHOT(1,:)/max(blocks(b).PHOT(1,:));
        if blocks(b).dataIsMulti == 0
        PHOT2 = -blocks(b).PHOT(2,:)/max(blocks(b).PHOT(2,:)); %make -blocks if Melvyn rig
        else
        PHOT2 = blocks(b).PHOT(2,:)/max(blocks(b).PHOT(2,:)); %No inversion for multichannel data (Currently)
        end
        PHOT3 = blocks(b).PHOT(3,:)/max(blocks(b).PHOT(3,:)); %If this crashes, data was not acquired on the multichannel setup, or was not properly specified as multichannel?

        %QA
        if ~any(~isnan(PHOT1)) || ~any(~isnan(PHOT2)) || ~any(~isnan(PHOT3))
            ['## No PHOT1, PHOT2, or PHOT3 data found! ##']
            crash = yes
        end

        if aux_plots
            photty1 = figure;
            plot(PHOT1)
            title('PHOT1')
            photty2 = figure;
            plot(PHOT2)
            title('PHOT2')
            photty3 = figure;
            plot(PHOT3)
            title('PHOT3')
        end

        PHOT1 = movmax(PHOT1,photMovMaxWindow);%was [40 40]
        PHOT2 = movmax(PHOT2,photMovMaxWindow);%was [40 40]
        PHOT3 = movmax(PHOT3,photMovMaxWindow);%was [40 40]
        
        blocks(b).PHOT(1,:) = PHOT1;
        blocks(b).PHOT(2,:) = PHOT2;
        blocks(b).PHOT(3,:) = PHOT3;

        %find beginning of stimuli
        LOCS_PHOT1 = find(diff(PHOT1 > peakThreshold) > 0) + 1;
        LOCS_PHOT2 = find(diff(PHOT2 > peakThreshold) > 0) + 1;
        LOCS_PHOT3 = find(diff(PHOT3 > peakThreshold) > 0) + 1;

        %first round of double peak detection
            %Note: This is 'double' peaks insofar as double of one channel, not 'double' as in multiple simultaneous phot signals
        badLOCS_PHOT1 = LOCS_PHOT1([false diff(LOCS_PHOT1) < (0.8*ISI*resampleFreq)]); %Was 0.8
        badLOCS_PHOT2 = LOCS_PHOT2([false diff(LOCS_PHOT2) < (0.8*ISI*resampleFreq)]);
        badLOCS_PHOT3 = LOCS_PHOT3([false diff(LOCS_PHOT3) < (0.8*ISI*resampleFreq)]);

        if aux_plots
            %figure; hold on; plot(PHOT1); scatter(LOCS_PHOT1,zeros(size(LOCS_PHOT1)),'b','filled'); scatter(badLOCS_PHOT1,zeros(size(badLOCS_PHOT1)),'m','filled'); title(['Phot 1 data (Max:',num2str(max(blocks(b).PHOT(1,:))),')'])
            %figure; hold on; plot(PHOT2); scatter(LOCS_PHOT2,zeros(size(LOCS_PHOT2)),'r','filled'); scatter(badLOCS_PHOT2,zeros(size(badLOCS_PHOT2)),'m','filled'); title(['Phot 2 data (Max:',num2str(max(blocks(b).PHOT(2,:))),')'])
            figure(photty1); hold on; plot(PHOT1); scatter(LOCS_PHOT1,zeros(size(LOCS_PHOT1)),'b','filled'); scatter(badLOCS_PHOT1,zeros(size(badLOCS_PHOT1)),'m','filled'); title(['Phot 1 data (Max:',num2str(max(blocks(b).PHOT(1,:))),')'])
            figure(photty2); hold on; plot(PHOT2); scatter(LOCS_PHOT2,zeros(size(LOCS_PHOT2)),'r','filled'); scatter(badLOCS_PHOT2,zeros(size(badLOCS_PHOT2)),'m','filled'); title(['Phot 2 data (Max:',num2str(max(blocks(b).PHOT(2,:))),')'])
            figure(photty3); hold on; plot(PHOT3); scatter(LOCS_PHOT3,zeros(size(LOCS_PHOT3)),'r','filled'); scatter(badLOCS_PHOT3,zeros(size(badLOCS_PHOT3)),'m','filled'); title(['Phot 3 data (Max:',num2str(max(blocks(b).PHOT(3,:))),')'])

            %Plot grand merged phot data
            figure
            hold on
            plot(PHOT1,'Color','r')
            plot(PHOT2,'Color','g')
            plot(PHOT3,'Color','b')
            title('All PHOT data')
            legend({'PHOT1','PHOT2','PHOT3'},'AutoUpdate','off')
            scatter(LOCS_PHOT1,zeros(size(LOCS_PHOT1)),'r','filled')
            scatter(LOCS_PHOT2,zeros(size(LOCS_PHOT2)),'g','filled')
            scatter(LOCS_PHOT3,zeros(size(LOCS_PHOT3)),'b','filled')
        end

        %-------------------------------------

        %Bespoke calculations
        %Multiple stimulus system
        multiLOCS = zeros(3,length(blocks(b).PHOT)); %Note: Unlike LOCS, this has n_back rows, not one row; Also, it relates to raw photodiode events, not eventual locs
        allLOCS = nan( abs(n_back) , nanmax([size(LOCS_PHOT1,2),size(LOCS_PHOT2,2),size(LOCS_PHOT3,2)]) );
        allLOCS(1,1:size(LOCS_PHOT1,2)) = LOCS_PHOT1;
        allLOCS(2,1:size(LOCS_PHOT2,2)) = LOCS_PHOT2;
        allLOCS(3,1:size(LOCS_PHOT3,2)) = LOCS_PHOT3;
        %Thus, allLOCS is a full list of all detected LOCS
            %Of critical note though is the fact that there is no guaranteed order synchrony (e.g. 1st apparent PHOT2 event may be after true 7th PHOT1 event etc)
        %Assign all LOCS to array
        for photi = 1:size(allLOCS,1) %Channel
            thisLOCS = allLOCS(photi,:); %This is because zeroes can creep back in
            thisLOCS( thisLOCS == 0 ) = NaN;
            for phote = 1:find( ~isnan(thisLOCS), 1, 'last')%size(thisLOCS,2) %Event
                multiLOCS( photi, thisLOCS(phote) ) = multiLOCS( photi, thisLOCS(phote) ) + 1; %Unambiguously denote this channel as having occurred here
            end
        end
        %Reiterate and find simultaneous events/etc
        windowActual = floor(transMergeWindow*ISI*resampleFreq); %Keeps consistency with multiple uses
        allINDS = find( multiLOCS == 1 ); %Note: Linear indicised for multidimensional matrix
        mergeLOCS = zeros( abs(n_back), length(blocks(b).PHOT) ); 
        diffMag = nan(1,size(allINDS,1));
        oblitInds = nan( size(multiLOCS,1) , size(allINDS,1)); %This is ostensibly nPhotChannels x nEvents, but its values will be derived from allINDs reference frame
        engramFailureCount = 0;
        for phote = 1:size(allINDS,1)
            %Check if already processed
            if ismember( allINDS(phote) , oblitInds )
                continue
            end

            %Reminder: allINDS reference frame is *linearised* full data length
            [r,c] = ind2sub( [size(multiLOCS)], allINDS(phote) ); %R - PHOT channel, C - Index
                %Reminder: c reference frame is full data length
            coords = [ c-windowActual : c+windowActual ];
            coords( coords < 1 ) = []; 
            coords( coords > length(blocks(b).PHOT) ) = [];
            temp = nansum( multiLOCS( :, coords ) , 'all'); %If 1, presumably lone event, if 2, presumably simultaneous event, if 3, presumably triple simultaneous event (Not actually coded in PTB yet), if more, overdetection

            %{{ %Activate this for debugging to suppress everything except diffMag calcs
            %QA
            if temp > 3 
                ['-# Alert: Extreme probability of transition probability merge window being too large #-']
                crash = yes; %Probably overkill
            end
    
            thisEngram = nansum( multiLOCS( :, coords ) , 2);
            %Secondary QA
            if any(thisEngram == 2) %Only possible if more than one apparent event of a particular channel found in a window
                disp( ['-# Phot merge engram failure #-'] )
                %crash = yes
                engramFailureCount = engramFailureCount + 1;
                continue
            end
    
            [yen,ind] = ismember( thisEngram', photTransIndex, 'rows' );
            %Tertiary QA
                %Note: This is basically the only QA standing between 'fake' engram finding and data being correct
            if yen == 0 || numel(ind) > 1
                disp( ['-# Alert: Engram under/overfind in phot possibility index #-'] )
                %crash = yes %Likely if index badly specified or data case unknown
                %Note: May actually be possible with minor anomalies etc?
                engramFailureCount = engramFailureCount + 1;
                continue
            end
            mergeLOCS( ind, c ) = c; %Assign to the photodiode case, rather than channel number
                %And also, as above, assign the index to the index position, because it will be logical'd in a moment
    
            %}
            %Store information about gaps, for posterity
            [subr,subc] = ind2sub( [3,numel(coords)], find(multiLOCS( :, coords )) ); %Find the column distance between phot events across channels (Probably 0 if engram only contains one event?)
                %This line finds 1s in an array of size 3x<numel(coords)>, with a reference frame relative to the first instance of coords
                %Thus, subc reference frame is relative to c-floor(transMergeWindow*ISI*resampleFreq)?
                    %Hence, at least one element should correspond to c
            diffMag(phote) = nanmean( abs(diff(subc)) ); 

            %Obliterate events collected as part of this engram
            subrCorrected = subr; %No actual correction needed, but maintain synchrony
            subcCorrected = subc + c - windowActual - 1; %Will be vector if more than one event in window
                %Also -1 because reasons
            subcCorrected( subr == r ) = []; %Shield original event from being obliterated
            subrCorrected( subr == r ) = []; %Ditto
            [oblind] = sub2ind( [size(multiLOCS)], subrCorrected, subcCorrected ); %Back-calculate allINDS index
            oblitInds(subrCorrected, phote) = oblind;
    
        end

        %Report
        temp = oblitInds;
        temp( isnan(temp) ) = 0;
        temp = logical(temp);
        disp([num2str(nansum(temp,'all')),' events (of ',num2str(size(allINDS,1)),' total detections) were obliterated because of simultaneity reasons'])
            %In theory this should match the number of instances of stimuli denoted by 1 or more simultaneous phot HIGHs
        if engramFailureCount ~= 0
            disp(['Engram failures occurred ',num2str(engramFailureCount),' time/s'])
        end

        %Logicalise
        mergeLOCS = logical(mergeLOCS);
        mergeLOCS = int16(mergeLOCS); %If we do not do this later operations won't work
    
        %QA for potential overfinding
            %(Critically important)
        if any( nansum(mergeLOCS,1) > 1 ) %Indicates whether one point has two or more stimuli apparently having happened there
            ['-# Alert: Potential overfind of stimulus #-']
            crash = yes
        end
        %Note: This QA is the only reason that the sum to follow is valid
    
        %Report
        disp(['Calculated stimulus counts:'])
        disp( num2str(nansum(mergeLOCS,2)) )
    
        %Hist to show apparent accuracy of merging
        if aux_plots
            figure
            hist(diffMag, floor([transMergeWindow*ISI*resampleFreq]) )
            xlim([0,floor([transMergeWindow*ISI*resampleFreq])])
            title(['Simultaneous phot event merge distance (Window: ',num2str(floor([transMergeWindow*ISI*resampleFreq])),'x2 of ',num2str(ISI*resampleFreq),' samples / ',num2str(transMergeWindow*2*100),'%)'])
        end
    
        %Assemble into split form
        %splitLOCS = nan( abs(n_back), nanmax(nansum(mergeLOCS,2)) );

        %transLOCS = nan(1, length(blocks(b).PHOT));
        transLOCS = sum(mergeLOCS,1);
        transLOCS = find(transLOCS == 1); %This mirrors the format/architecture of LOCS under normal circumstances
            %However, it does not contain the stimulus identity, and unlike above, that is not (directly) derivable from LOCS_PHOT1,2 etc

        %Flatten mergeLOCS to one-dimensional
        for stimi = 1:size(mergeLOCS,1)
            mergeLOCS( stimi, find(mergeLOCS(stimi,:)) ) = stimi; %Only works if mergeLOCS back to being int or double
        end
        mergeLOCS = sum(mergeLOCS,1);
        mergeLOCS = mergeLOCS( logical(mergeLOCS) ); %Borrow trick from inferRandomSequence to reduce, but don't change values
            %At the end of this chain, mergeLOCS values will indicate stim identity, with position indicating index in randomSequence
                %(True indices in data are contained in transLOCS, which should be the same length)

        %QA for synchrony between transLOCS and mergeLOCS
        if numel(transLOCS) ~= numel(mergeLOCS)
            ['## Alert: Critical synchrony failure between stimulus indices and identities ##']
            crash = yes; %Not sure how this could happen, but just in case
        end
    
    
        if isempty(mergeLOCS)
            ['## ALERT: NO PHOT LOCS FOUND ##']
            crash = yes %May actually be normal for some conditions?
        end
        
        blocks(b).LOCS = transLOCS; %Alter form of LOCS, given same name for simplicity
        blocks(b).mergeLOCS = mergeLOCS; %New variable that denotes identity (akin to the later-calculated randomSequence)
        blocks(b).LOCS_PHOT1 = LOCS_PHOT1; %Unlike above, these do not unambiguously indicate stimulus identity
        blocks(b).LOCS_PHOT2 = LOCS_PHOT2;
        blocks(b).LOCS_PHOT3 = LOCS_PHOT3;
        blocks(b).transProbDesign = 1;
        blocks(b).transProbAncillary.nStimuli = abs(n_back);


    end %nBack end
    
    

end

end

