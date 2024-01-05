% this function processes a set of blocks, usually corresponding to a
% single fly and condition; it concatenates the data for individual blocks
% and conditions (e.g. LIT/DARK)
function R = processBlocks(blocks, aux_plots, mode)

blocks = calculatePeaks(blocks, aux_plots);

blocks = inferRandomSequence(blocks);

blocks = calculateBadTrials(blocks, aux_plots);

blocks = findFocusLocs(blocks);

% peak detection figure
if aux_plots

    for b = 1:length(blocks)
        figure;
        hold on

        if blocks(b).PHOTType == 1
            plot(blocks(b).PHOT(1,:));% PHOT3
            scatter(blocks(b).LOCS_PHOT1,blocks(b).PKS_PHOT1);
            scatter(blocks(b).LOCS_PHOT2,blocks(b).PKS_PHOT2);
        elseif blocks(b).PHOTType == 2
            plot(blocks(b).PHOT(1,:)); plot(blocks(b).PHOT(2,:));
            scatter(blocks(b).LOCS_PHOT1,zeros(size(blocks(b).LOCS_PHOT1)),'b','filled');
            scatter(blocks(b).LOCS_PHOT2,zeros(size(blocks(b).LOCS_PHOT2)),'r','filled');
        end

        if isfield(blocks,'focusPeaks')
            scatter(blocks(b).focusLocs, zeros(size(blocks(b).focusLocs)),40,'m','filled');
        else
            scatter(blocks(b).badLOCS, zeros(size(blocks(b).badLOCS)),40,'m','filled');
        end
    end

end

switch mode
    case 'time'
        R = analyseSequentialEffects(blocks, aux_plots);
    case 'frequency'
        R = analyseSequentialEffectsFrequency(blocks, aux_plots);
    case 'timefrequency'
        R = analyseSequentialEffectsTimeFrequency(blocks, aux_plots);
end
