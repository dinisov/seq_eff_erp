function R = processBlocks(blocks, aux_plots)
%processBlocks Summary of this function goes here
% this function processes a set of blocks, usually corresponding to a
% single fly and condition; it concatenates the data for individual blocks
% and conditions (e.g. LIT/DARK)

blocks = calculatePeaks(blocks, aux_plots);

blocks = inferRandomSequence(blocks);

blocks = calculateBadTrials(blocks, aux_plots);

% only used for block experiments (automatic)
blocks = findFocusLocs(blocks);

% peak detection figure (different depending on photodiode method)
if aux_plots
    peakDetectionFigure(blocks);
end

R = analyseSequentialEffects(blocks, aux_plots);

end