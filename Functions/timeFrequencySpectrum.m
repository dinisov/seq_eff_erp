function R = timeFrequencySpectrum(R, blocks)
%timeFrequencyAnalysis Summary of this function goes here
%   Detailed explanation goes here
    cwtERPs = [];
    
    %300 Hz seems like a sensible upper limit
    for s = 1:16
       [wt,f] = cwt(R.meanERPs(:,s),blocks(1).resampleFreq,'FrequencyLimits',[0 300]);
       cwtERPs(:,s,:) = wt; %#ok<AGROW>
    end
    
    magnitudeSEs = abs(cwtERPs);
    phaseSEs = angle(cwtERPs);
     
    R.magnitudeSEs = magnitudeSEs;
    R.phaseSEs = phaseSEs;
    R.f = f;
end