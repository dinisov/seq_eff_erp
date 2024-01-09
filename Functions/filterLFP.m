function LFP = filterLFP(LFP, resampleFreq)
%filterLFP Summary of this function goes here
%   Detailed explanation goes here

    % butterworth filter LFP
    [b_f,a_f] = butter(6,100/resampleFreq*2);
    LFP = filter(b_f,a_f,LFP.').';

    %notch filter at 50Hz
    wo = 50/(resampleFreq/2);  
    bw = wo/10;
    [b_f,a_f] = iirnotch(wo,bw);

    LFP = filter(b_f,a_f,LFP.').';
    
    % savitsky-golay filter
%     LFP = smoothdata(LFP,'sgolay');

end