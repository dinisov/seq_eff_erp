function R = processBlocks(blocks, aux_plots, plotSelector, reOrder, n_back, options)
%processBlocks Summary of this function goes here
% this function processes a set of blocks, usually corresponding to a
% single fly and condition; it concatenates the data for individual blocks
% and conditions (e.g. LIT/DARK)
arguments
    blocks struct
    aux_plots double
    plotSelector double
    reOrder double
    n_back double
    options.suppressANOVA double = 0
    options.plotIndividualFlies double = 1
    options.scramLevel double = 0
end

%Extra
suppressANOVA = options.suppressANOVA;
plotIndividualFlies = options.plotIndividualFlies;
scramLevel = options.scramLevel;

blocks = calculatePeaks(blocks, aux_plots);

blocks = inferRandomSequence(blocks);

blocks = calculateBadTrials(blocks, aux_plots);

% only used for block experiments (automatic)
blocks = findFocusLocs(blocks);

% peak detection figure (different depending on photodiode method)
if aux_plots
    peakDetectionFigure(blocks);
end

if scramLevel == 1
    blocks.randomSequence = blocks.randomSequence( randperm( numel(blocks.randomSequence) ) );
    disp(['Data scrambled at raw sequence level (Scram 1)'])
end

%R = analyseSequentialEffects(blocks, aux_plots);
%R = analyseSequentialEffects(blocks, aux_plots, plotSelector, reOrder, n_back);
R = analyseSequentialEffects(blocks, aux_plots, plotSelector, reOrder, n_back, ...
    'suppressANOVA',suppressANOVA, 'plotIndividualFlies',plotIndividualFlies,...
    'scramLevel',scramLevel);

end