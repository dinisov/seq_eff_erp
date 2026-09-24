function [allERPs, allPHOTs, goodTrials, focusPeaks, allTIMEs] = groupBlocks(blocks,n_back)
%groupBlocks Summary of this function goes here
%   Detailed explanation goes here

    %n_seq = 2^n_back; %Moved lower, for transition probability reasons

    n_blocks = length(blocks);

    total_length = 0;

    transProbDesign = [];
    for b = 1:n_blocks
        total_length = total_length + size(blocks(b).ERPS,3);
        transProbDesign = [transProbDesign,blocks(b).transProbDesign];
    end 
    %QA for variation in block type
    if numel(unique(transProbDesign)) > 1
        ['## Alert: Potential mix of original and transition probability blocks ##']
        crash = yes %It is currently not supported to swap between these in one analysis run
    end
    transProbDesign = unique(transProbDesign);

    if ~transProbDesign
        n_seq = 2^n_back;
    else
        nStimuli = [];
        nBackActual = [];
        for b = 1:n_blocks
            nStimuli = [nStimuli,blocks(b).transProbAncillary.nStimuli];
            nBackActual = [nBackActual,blocks(b).transProbAncillary.nBackActual];
        end
        %QA
        if numel(unique(nStimuli)) > 1 || numel(unique(nBackActual)) > 1
            ['## Cannot analyse transition probability experiments with varying numbers of stimuli AND/OR nBack ##']
            crash = yes
        end
        nStimuli = unique(nStimuli);
        nBackActual = unique(nBackActual);
        n_seq = nStimuli^nBackActual;
        %Note: If there is variability between 
    end
    
    window = floor(blocks(1).window * blocks(1).resampleFreq); %Assumes first block is representative
    
    % concatenate ERPs from different "experiments" (blocks)
    allERPs = zeros(length(window(1):window(2)), n_seq, total_length);
    allPHOTs = zeros(length(window(1):window(2)), n_seq, total_length);
    allTIMEs = zeros(length(window(1):window(2)), n_seq, total_length);
    
    start_index = 0; goodTrials = [];

    % these are already separated by sequence so in order to group by block
    % it is only necessary to stack along third dimension
    for b = 1:n_blocks

        allERPs(:,:,start_index + 1:start_index + size(blocks(b).ERPS,3)) = blocks(b).ERPS;
        allPHOTs(:,:,start_index+ 1:start_index + size(blocks(b).ERPS,3)) = blocks(b).seqPHOT;
        allTIMEs(:,:,start_index+ 1:start_index + size(blocks(b).ERPS,3)) = blocks(b).seqTIME;

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

