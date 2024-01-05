function [allERPs, allPHOTs, goodTrials, focusPeaks] = groupBlocks(blocks,window,n_back)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    n_seq = 2^n_back;

    n_blocks = length(blocks);

    total_length = 0;

    for b = 1:n_blocks
        total_length = total_length + size(blocks(b).ERPS,3);
    end 
    
    % concatenate ERPs from different "experiments" (blocks)
    allERPs = zeros(length(window(1):window(2)), n_seq, total_length);
    allPHOTs = zeros(length(window(1):window(2)), n_seq, total_length);
    
    start_index = 0; goodTrials = [];

    % these are already separated by sequence so in order to group by block
    % it is only necessary to stack along third dimension
    for b = 1:n_blocks

        allERPs(:,:,start_index + 1:start_index + size(blocks(b).ERPS,3)) = blocks(b).ERPS;
        allPHOTs(:,:,start_index+ 1:start_index + size(blocks(b).ERPS,3)) = blocks(b).seqPHOT;

        start_index = start_index + size(blocks(b).ERPS,3);

        goodTrials = [goodTrials 1-blocks(b).badTrials]; %#ok<AGROW> 
        
    end
    
    if isfield(blocks,'focusPeaks')
        focusPeaks = [];
        for b = 1:n_blocks
            focusPeaks = [focusPeaks blocks(b).focusPeaks]; %#ok<AGROW> 
        end
    else
        focusPeaks = [];
    end

end

