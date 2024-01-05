function blocks = inferRandomSequence(blocks)
%inferRandomSequence Summary of this function goes here
%   Detailed explanation goes here

for b = 1:length(blocks)
   
    % infer random sequence (0 - left; 1 - right)
    randomSequence = zeros(1,length(blocks(b).PHOT));
    randomSequence(blocks(b).LOCS_PHOT1) = 2; randomSequence(blocks(b).LOCS_PHOT2) = 1;
    
    randomSequence = randomSequence(logical(randomSequence)) - 1;
    
    blocks(b).randomSequence = randomSequence;
end

end

