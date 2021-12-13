% this function processes a set of blocks, usually corresponding to a
% single fly and condition; it concatenates the data for individual blocks
% and con
function R = processBlocks(blocks, aux_plots)
    
    n_blocks = length(blocks);

    for b = 1:n_blocks
        
        PHOT = blocks(b).PHOT; 
        resampleFreq = blocks(b).resampleFreq;
        ISI = blocks(b).ISI;
        peakThreshold = blocks(b).peakThreshold;

        % find peaks (sd of 3 seems to work very well after percentile based cull)
        [PKS_PHOT1,LOCS_PHOT1] = findpeaksbase(normalize(PHOT(1,:)) , 'MinPeakHeight' , peakThreshold , 'MinPeakDistance' , 1/2*ISI*resampleFreq );
        [PKS_PHOT2,LOCS_PHOT2] = findpeaksbase(normalize(PHOT(2,:)) , 'MinPeakHeight' , peakThreshold , 'MinPeakDistance' , 1/2*ISI*resampleFreq );

        % plot peak detection for each photodiode as a sanity check
        if aux_plots
            figure
            plot(normalize(PHOT(1,:)));
            hold on
            scatter(LOCS_PHOT1,PKS_PHOT1)
            xlabel('Index')
            title('Detected peaks in phot1')

            figure
            plot(normalize(PHOT(2,:)));
            hold on
            scatter(LOCS_PHOT2,PKS_PHOT2)
            xlabel('Index')
            title('Detected peaks in phot2')
        end
            
        % fuse locations of PHOT1 and PHOT2 (I figured this was quicker than concatenating and sorting)
        LOCS = zeros(size(PHOT(1,:)));
        LOCS(LOCS_PHOT1) = LOCS_PHOT1; LOCS(LOCS_PHOT2) = LOCS_PHOT2;
        LOCS = LOCS(logical(LOCS));

        %we must get rid of trials where we could not get a peak and the
        %subsequent four trials
        badLOCS = find(diff(LOCS) > 1.5*ISI*resampleFreq) + 1; % index of trials where gap was longer than ISI
        badLOCS = LOCS(badLOCS);
        badTrials = zeros(size(PHOT(1,:)));
        badTrials(badLOCS) = 1;

        % infer random sequence (0 - left; 1 - right)
        randomSequence = zeros(size(PHOT(1,:)));
        randomSequence(LOCS_PHOT1) = 2; randomSequence(LOCS_PHOT2) = 1;
        badTrials = badTrials(logical(randomSequence)); % careful order is important here
        randomSequence = randomSequence(logical(randomSequence)) - 1;

        %remove 4 trials after a bad one
        badTrialsIndex = find(badTrials);
        badTrials([badTrialsIndex+1 badTrialsIndex+2 badTrialsIndex+3 badTrialsIndex+4]) = 1;

        % histogram of interval between peaks (should have one tight peak)
%         if aux_plots
            figure;
            histogram(diff(LOCS));
%         end
    
        % add processed data to original blocks structure
        blocks(b).badTrials = badTrials;
        blocks(b).LOCS = LOCS;
        blocks(b).randomSequence = randomSequence;
            
    end

R = analyseSequentialEffects(blocks, aux_plots);