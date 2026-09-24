function blocks = inferRandomSequence(blocks)
%inferRandomSequence Summary of this function goes here
%   Detailed explanation goes here

for b = 1:length(blocks)
   
    if blocks(b).transProbDesign ~= 1
        % infer random sequence (0 - left; 1 - right)
        randomSequence = zeros(1,length(blocks(b).PHOT)); %randomSequence starts as length n(Samples)
        randomSequence(blocks(b).LOCS_PHOT1) = 2; randomSequence(blocks(b).LOCS_PHOT2) = 1;
        
        randomSequence = randomSequence(logical(randomSequence)) - 1; %Here randomSequence becomes length n(Stim events) and value [0,1]
            %The above, rather obtuse line works by using the fact that logical(1) and logical(2) are both 1, and therefore randomSequence is
            %reduced to only instances of '1' and '2', which are then instantaneously converted to 0 and 1, respectively
        
        blocks(b).randomSequence = randomSequence;
    else
        disp(['Inferring sequence for ',num2str(blocks(b).transProbAncillary.nStimuli),' stimuli'])
        randomSequence = blocks(b).mergeLOCS; %Pull directly from calculatePeaks product

        blocks(b).randomSequence = randomSequence;
    end

end

end

