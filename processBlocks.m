% this function processes a set of blocks, usually corresponding to a
% single fly and condition; it concatenates the data for individual blocks
% and conditions (LIT/DARK)
function R = processBlocks(blocks, aux_plots)
    
    n_blocks = length(blocks);

    for b = 1:n_blocks
        
        PHOT1 = blocks(b).PHOT(1,:)/max(blocks(b).PHOT(1,:));
        PHOT2 = blocks(b).PHOT(2,:)/max(blocks(b).PHOT(2,:)); 
        resampleFreq = blocks(b).resampleFreq;
        ISI = blocks(b).ISI;
%         disp(resampleFreq);
        PHOT1 = movmax(PHOT1,[20 20]);
        PHOT2 = movmax(PHOT2,[20 20]);
        blocks(b).PHOT(1,:) = PHOT1;
        blocks(b).PHOT(2,:) = PHOT2;

%         figure; plot(PHOT2); hold on; plot(PHOT1);

        %find beginning of stimuli
        LOCS_PHOT1 = find(diff(PHOT1 > .4) > 0) + 1;
        LOCS_PHOT2 = find(diff(PHOT2 > .4) > 0) + 1;

%         figure; hold on; plot(PHOT1); scatter(LOCS_PHOT1,zeros(size(LOCS_PHOT1)),'filled');
%         figure; hold on; plot(PHOT2); scatter(LOCS_PHOT2,zeros(size(LOCS_PHOT2)),'filled'); 
        
        % fuse locations of PHOT1 and PHOT2 (I figured this was quicker than concatenating and sorting)
        LOCS = zeros(size(PHOT1));
        LOCS(LOCS_PHOT1) = LOCS_PHOT1; LOCS(LOCS_PHOT2) = LOCS_PHOT2;
        LOCS = LOCS(logical(LOCS));
        
        %we must get rid of trials where we could not get a peak and the
        %subsequent four trials
        badLOCS = LOCS([false diff(LOCS) > (1.2*ISI*resampleFreq)] | [false diff(LOCS) < (0.8*ISI*resampleFreq)]); % index of trials where gap was too long or too short

        % infer random sequence (0 - left; 1 - right)
        randomSequence = zeros(size(PHOT1(1,:)));
        randomSequence(LOCS_PHOT1) = 2; randomSequence(LOCS_PHOT2) = 1;

        %create logical vector of which trials are bad
        badTrials = zeros(size(PHOT1(1,:)));
        badTrials(badLOCS) = 1;

        badTrials = badTrials(logical(randomSequence));
        randomSequence = randomSequence(logical(randomSequence)) - 1;

        %add four trials subsequent to the bad trials vector
%         indBadTrials = find(badTrials);
%         badTrials([indBadTrials+1 indBadTrials+2 indBadTrials+3 indBadTrials+4]) = 1;
        
        percentDataLost = nnz(badTrials)/length(badTrials);
        disp(['Data lost due to bad peak detection: ' num2str(percentDataLost*100) '%']);
        
        % histogram of interval between peaks (should have one tight peak)
        % this is a critical check so it is always plotted
        figure;
        histogram(diff(LOCS));
    
        % add processed data to original blocks structure
        blocks(b).badTrials = badTrials;
        blocks(b).LOCS = LOCS;
        blocks(b).randomSequence = randomSequence;
        
        % peak detection figure
        if aux_plots
            figure
            hold on
            plot(PHOT1); plot(PHOT2);
            scatter(LOCS_PHOT1,zeros(size(LOCS_PHOT1)),'b','filled');
            scatter(LOCS_PHOT2,zeros(size(LOCS_PHOT2)),'r','filled');

            scatter(badLOCS, zeros(size(badLOCS)),40,'m','filled');
        end
            
    end

R = analyseSequentialEffects(blocks, aux_plots);