function blocks = calculateBadTrials(blocks, aux_plots)
%calculateBadTrials Summary of this function goes here
%   Detailed explanation goes here

for b = 1:length(blocks)
    
        ISI = blocks(b).ISI;
        resampleFreq = blocks(b).resampleFreq;
        LOCS = blocks(b).LOCS;
   
        %we must get rid of trials where we could not get a peak and the
        %subsequent four trials
        badLOCS = LOCS([false diff(LOCS) > (1.5*ISI*resampleFreq)] | [false diff(LOCS) < (0.5*ISI*resampleFreq)]); % index of trials where gap was too long or too short

        %create logical vector of which trials are bad
        badTrials = zeros(1,length(blocks(b).PHOT));
        badTrials(badLOCS) = 1;

        aux = zeros(1,length(blocks(b).PHOT)); aux(LOCS) = 1;
        badTrials = badTrials(logical(aux));

        %add four trials subsequent to the bad trials vector
        indBadTrials = find(badTrials);
        badTrials([indBadTrials+1 indBadTrials+2 indBadTrials+3 indBadTrials+4]) = 1;

        percentDataLost = nnz(badTrials)/length(badTrials);
        disp(['Data lost due to bad peak detection: ' num2str(percentDataLost*100) '%']);

        % add processed data to original blocks structure
        blocks(b).badLOCS = badLOCS;
        blocks(b).badTrials = badTrials;
        
        % histogram of interval between peaks (should have one tight peak)
        if aux_plots
            figure;
            histogram(diff(LOCS(~badTrials)));title('Good peaks histogram');
        end
    
end

end

