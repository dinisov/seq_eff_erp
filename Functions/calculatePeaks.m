function blocks = calculatePeaks(blocks, aux_plots)
%calculatePeaks Calculate photodiode peaks
%   Detailed explanation goes here

for b = 1:length(blocks)
    
    peakThreshold = blocks(b).peakThreshold;
    ISI = blocks(b).ISI;
    resampleFreq = blocks(b).resampleFreq;
    
    % one photodiode (different size peaks)
    if blocks(b).PHOTType == 1

        PHOT = -blocks(b).PHOT(3,:)/max(blocks(b).PHOT(3,:));

        if aux_plots
            figure; plot(PHOT);
        end
        PHOT = movmax(PHOT,[20 20]);
        if aux_plots    
            hold on; plot(PHOT);
        end

        blocks(b).PHOT = [PHOT; PHOT; PHOT];

        % find peaks
        [PKS_PHOT1,LOCS_PHOT1] = findpeaksbase(PHOT, 'MinPeakHeight' , .1 , 'MinPeakDistance' , 1/2*ISI*resampleFreq );
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

        if aux_plots
            photty1 = figure;
            plot(PHOT1)
            title('PHOT1')
            photty2 = figure;
            plot(PHOT2)
            title('PHOT2')
        end

        PHOT1 = movmax(PHOT1,[40 40]);%was [40 40]
        PHOT2 = movmax(PHOT2,[40 40]);%was [40 40]
        
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
    
    blocks(b).LOCS = LOCS;
    blocks(b).LOCS_PHOT1 = LOCS_PHOT1;
    blocks(b).LOCS_PHOT2 = LOCS_PHOT2;

end

end

