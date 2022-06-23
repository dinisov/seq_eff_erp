% analyse in which frequency bands SLRP and LRPR live

clear; close all;

load slrp_lrpr.mat

load data_flies_frequency.mat

frequency_bounds = [0 60];

r_squared_slrp = [];
r_squared_lrpr = [];
r_squared_overall = [];

options = optimset('Algorithm','interior-point','FinDiffType','central');

resampleFreq = 3000;%CHECK

fly_index = 1;

for fly = 1:length(FLIES)
    
    if ~isempty(FLIES(fly).LIT)
   
        L = size(FLIES(fly).LIT.magnitudeSEs,1);
        
        for coeff = 1:ceil(L/2)
            
                freq = resampleFreq*(coeff-1)/L;
                
                if freq > frequency_bounds(1) && freq < frequency_bounds(2)
                    
                    seq_eff_pattern = FLIES(fly).LIT.magnitudeSEs(coeff,:).';
                    
                    %fit only to slrp or lrpr
                   [x_slrp,sse_slrp] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[-inf  0 -inf -inf],[inf 0 inf inf],[],options);
                   [x_lrpr,sse_lrpr] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 1 1 0],[],[],[],[],[0  -inf -inf -inf],[0 inf inf inf],[],options);
                    
                   %fit overall
                   [x_overall,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 1 0 0],[],[],[],[],[-inf  -inf 0 -inf],[inf inf 0 inf],[],options);

                   
                   sse_total = sum((seq_eff_pattern-mean(seq_eff_pattern)).^2);
                   
                   r_squared_slrp(fly_index,coeff) = 1-(sse_slrp/sse_total); %#ok<SAGROW>
                   r_squared_lrpr(fly_index,coeff) = 1-(sse_lrpr/sse_total); %#ok<SAGROW>
                   r_squared_overall(fly_index,coeff) = 1-(sse_overall/sse_total); %#ok<SAGROW>
                   
                end
                
        end
        
        fly_index = fly_index+1;
        
    end
    
end

imagesc(r_squared_slrp); colorbar;
figure; imagesc(r_squared_lrpr); colorbar;
figure; imagesc(r_squared_overall); colorbar;

%waterfall plots
figure; waterfall(r_squared_slrp);
figure; waterfall(r_squared_lrpr);

%% Grouped analysis

