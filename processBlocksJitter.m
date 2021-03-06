% this function processes a set of blocks, usually corresponding to a
% single fly and condition; it concatenates the data for individual blocks
% and conditions (LIT/DARK)
%JITTER VERSION (1 channel)
function R = processBlocksJitter(blocks, aux_plots)
    
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
        
        % fuse locations of PHOT1 and PHOT2 (I figured this was quicker than concatenating and sorting)
        LOCS = zeros(size(PHOT));
        LOCS(LOCS_PHOT1) = LOCS_PHOT1; LOCS(LOCS_PHOT2) = LOCS_PHOT2;
        LOCS = LOCS(logical(LOCS));
        
        %we must get rid of trials where we could not get a peak and the
        %subsequent four trials
        badLOCS = LOCS([false diff(LOCS) > (1.2*ISI*resampleFreq)] | [false diff(LOCS) < (0.8*ISI*resampleFreq)]); % index of trials where gap was too long or too short

        % infer random sequence (0 - left; 1 - right)
        randomSequence = zeros(size(PHOT(1,:)));
        randomSequence(LOCS_PHOT1) = 2; randomSequence(LOCS_PHOT2) = 1;

        %create logical vector of which trials are bad
        badTrials = zeros(size(PHOT(1,:)));
        badTrials(badLOCS) = 1;

        badTrials = badTrials(logical(randomSequence));
        randomSequence = randomSequence(logical(randomSequence)) - 1;

        %add four trials subsequent to the bad trials vector
        indBadTrials = find(badTrials);
        badTrials([indBadTrials+1 indBadTrials+2 indBadTrials+3 indBadTrials+4]) = 1;
        
        percentDataLost = nnz(badTrials)/length(badTrials);
        disp(['Data lost due to bad peak detection: ' num2str(percentDataLost*100) '%']);
        
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
            scatter(badLOCS, zeros(size(badLOCS)),40,'m','filled');
        end
            
    end

R = analyseSequentialEffectsJitter(blocks, aux_plots);