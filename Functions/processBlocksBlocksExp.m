% this function processes a set of blocks, usually corresponding to a
% single fly and condition; it concatenates the data for individual blocks
% and conditions (LIT/DARK)
function R = processBlocksBlocksExp(blocks, aux_plots)
    
    blocks = calculatePeaks(blocks, aux_plots);

    blocks = inferRandomSequence(blocks);

    blocks = calculateBadTrials(blocks, aux_plots);
    
    blocks = findFocusLocs(blocks);

    for b = 1:length(blocks)
        
        % peak detection figure
        if aux_plots
            figure
            hold on
            plot(blocks(b).PHOT(1,:));% PHOT3
            scatter(blocks(b).LOCS_PHOT1,blocks(b).PKS_PHOT1); % 

            scatter(blocks(b).LOCS_PHOT2,blocks(b).PKS_PHOT2);

              scatter(blocks(b).focusLocs, zeros(size(blocks(b).focusLocs)),40,'m','filled');
        end
            
    end

R = analyseSequentialEffectsBlocksExp(blocks, aux_plots);