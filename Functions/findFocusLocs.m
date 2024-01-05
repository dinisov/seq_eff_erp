function blocks = findFocusLocs(blocks)
%findFocusLocs Summary of this function goes here
%   Detailed explanation goes here

for b = 1:length(blocks)
    
    if ~isnan(blocks(b).blockLength)
        resampleFreq = blocks(b).resampleFreq;
        LOCS = blocks(b).LOCS;
        interBlockPeriod = blocks(b).interBlockPeriod;
        blockLength = blocks(b).blockLength;
        focusPeak = blocks(b).focusPeak;

        % focus on specific peak in train (starting at 5 up to blockLength)
        focusLocs = LOCS(find(diff(LOCS) > 0.8*interBlockPeriod*resampleFreq)-(blockLength-focusPeak));

        %create logical vector of sequence elements corresponding to last
        %peaks
        focusPeaks = zeros(1,length(blocks(b).PHOT));
        focusPeaks(focusLocs) = 1;
        aux = zeros(1,length(blocks(b).PHOT)); aux(LOCS) = 1;
        focusPeaks = focusPeaks(logical(aux));

        blocks(b).focusLocs = focusLocs;
        blocks(b).focusPeaks = focusPeaks;
    end
        
end
        
end

