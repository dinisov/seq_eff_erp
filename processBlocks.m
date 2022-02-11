% this function processes a set of blocks, usually corresponding to a
% single fly and condition; it concatenates the data for individual blocks
% and conditions (LIT/DARK)
%JITTER VERSION (1 channel)
function R = processBlocks(blocks, aux_plots)
    
    n_blocks = length(blocks);

    for b = 1:n_blocks
        
        PHOT = blocks(b).PHOT(3,:)/max(blocks(b).PHOT(3,:)); 
        resampleFreq = blocks(b).resampleFreq;
        ISI = blocks(b).ISI;

        PHOT = movmax(PHOT,[20 20]);
        blocks(b).PHOT(3,:) = PHOT;

        % find peaks
        [PKS_PHOT1,LOCS_PHOT1] = findpeaksbase(PHOT, 'MinPeakHeight' , .1 , 'MinPeakDistance' , 1/2*ISI*resampleFreq ); %just one stimulus
        [PKS_PHOT2,LOCS_PHOT2] = findpeaksbase(PHOT , 'MinPeakHeight' , .6 , 'MinPeakDistance' , 1/2*ISI*resampleFreq ); %all stimuli

        [LOCS_PHOT1, ind_locs_phot1] = setdiff(LOCS_PHOT1, LOCS_PHOT2);
        PKS_PHOT1 = PKS_PHOT1(ind_locs_phot1);
        
        % remove peak outliers
%         peakSD = 3;
%         LOCS_PHOT1 = LOCS_PHOT1(abs(normalize(PKS_PHOT1)) < peakSD);
%         LOCS_PHOT2 = LOCS_PHOT2(abs(normalize(PKS_PHOT2)) < peakSD);
%         PKS_PHOT1 = PKS_PHOT1(abs(normalize(PKS_PHOT1)) < peakSD);
%         PKS_PHOT2 = PKS_PHOT2(abs(normalize(PKS_PHOT2)) < peakSD);
        
        figure; histogram(PKS_PHOT1);
        figure; histogram(PKS_PHOT2);
        
        % fuse locations of PHOT1 and PHOT2 (I figured this was quicker than concatenating and sorting)
        LOCS = zeros(size(PHOT(1,:)));
        LOCS(LOCS_PHOT1) = LOCS_PHOT1; LOCS(LOCS_PHOT2) = LOCS_PHOT2;
        LOCS = LOCS(logical(LOCS));
        
        % 
%         LOCS = LOCS(2:end-1);
        
        %we must get rid of trials where we could not get a peak and the
        %subsequent four trials
        badLOCS = find((diff(LOCS) > (1.2*ISI*resampleFreq)) | (diff(LOCS) < (0.8*ISI*resampleFreq))) + 1; % index of trials where gap was too long or too short
        badLOCS = LOCS(badLOCS);
        badTrials = zeros(size(PHOT(1,:)));
        badTrials(badLOCS) = 1;

        % infer random sequence (0 - left; 1 - right)
        randomSequence = zeros(size(PHOT(1,:)));
        randomSequence(LOCS_PHOT1) = 2; randomSequence(LOCS_PHOT2) = 1;
        badTrials = badTrials(logical(randomSequence)); % careful order is important here
        randomSequence = randomSequence(logical(randomSequence)) - 1;
        
%         randomSequence = randomSequence(2:end-1);

        %remove 4 trials after a bad one
        badTrialsIndex = find(badTrials);
        badTrials([badTrialsIndex+1 badTrialsIndex+2 badTrialsIndex+3 badTrialsIndex+4]) = 1;
        
        percentDataLost = nnz(badTrials)/length(badTrials);
        disp(['Data lost: ' num2str(percentDataLost*100) '%']);
        
        % histogram of interval between peaks (should have one tight peak)
        % this is a critical check so it is always plotted
        figure;
        histogram(diff(LOCS(~badTrials)));
    
        % add processed data to original blocks structure
        blocks(b).badTrials = badTrials;
        blocks(b).LOCS = LOCS;
        blocks(b).randomSequence = randomSequence;
        
        % peak detection figure
        if aux_plots
            figure
            hold on
            plot(PHOT);% PHOT3
            scatter(LOCS_PHOT1,PKS_PHOT1); % 
            
%             plot(normalize(PHOT(2,:)));% PHOT2
            scatter(LOCS_PHOT2,PKS_PHOT2);

            % sanity check of where peaks were detected and which stimlus
            % (left or right)
%             scatter(LOCS(logical(randomSequence)), 0,40,'r','filled');
%             scatter(LOCS(~logical(randomSequence)), 0,40,'b','filled');

%             scatter([stimulusOnset1 stimulusOnset2], 0,40,'r','filled');
%             scatter([stimulusEnd1 stimulusEnd2], 0,40,'b','filled');
            scatter(badLOCS, 0,40,'m','filled');
        end
            
    end

R = analyseSequentialEffects(blocks, aux_plots);