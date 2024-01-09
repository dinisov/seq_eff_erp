function LFP = trimLFP(LFP,n_sd, aux_plots)
%trimLFP Summary of this function goes here
%   Detailed explanation goes here
    if ~isnan(n_sd)
    
        n_sd_lfp = n_sd;
        sd_lfp = std(LFP); mean_lfp = mean(LFP);
        
        % do a nice plot
        if aux_plots
            figure; plot(LFP); hold on; %#ok<*UNRCH>
            x_lim = xlim;
            plot([x_lim(1) x_lim(2)],[mean_lfp-n_sd_lfp*sd_lfp mean_lfp-n_sd_lfp*sd_lfp],'r');
            plot([x_lim(1) x_lim(2)],[mean_lfp+n_sd_lfp*sd_lfp mean_lfp+n_sd_lfp*sd_lfp],'r');
        end
            
        % nuke everything beyond +/- the set number of sd
        LFP(LFP > mean_lfp+n_sd_lfp*sd_lfp | LFP < mean_lfp-n_sd_lfp*sd_lfp) = NaN;

    end
end