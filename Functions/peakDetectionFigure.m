function peakDetectionFigure(blocks)
%peakDetectionFigure Summary of this function goes here
%   Detailed explanation goes here
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