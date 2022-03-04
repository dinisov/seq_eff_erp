function [exp_filter] = exp_filter(b,a,alpha)

all_seq =  [1 1 1 1 1
            1 0 0 0 0
            1 1 0 0 0
            1 0 1 1 1
            1 1 1 0 0
            1 0 0 1 1
            1 1 0 1 1
            1 0 1 0 0
            1 1 1 1 0
            1 0 0 0 1
            1 1 0 0 1
            1 0 1 1 0
            1 1 1 0 1
            1 0 0 1 0
            1 1 0 1 0
            1 0 1 0 1];

%% frequency filter

exponents = repmat([4 3 2 1],16,1);

exp_filt = all_seq(:,1:4).*(alpha.^exponents);

exp_filt = sum(exp_filt,2)/sum(alpha.^(1:4));

exp_filt(all_seq(:,5) == 0) = 1-exp_filt(all_seq(:,5) == 0);

%% alternation/repetition filter

all_seq = 1-abs(diff(all_seq,1,2));

exponents = repmat([3 2 1],16,1);

exp_filt_ar = all_seq(:,1:3).*(alpha.^exponents);

exp_filt_ar = sum(exp_filt_ar,2)/sum(alpha.^(1:3));

exp_filt_ar(all_seq(:,4) == 0) = 1-exp_filt_ar(all_seq(:,4) == 0);

%figure; create_seq_eff_plot(1-(exp_filt_ar+exp_filt));

%CAREFUL ONLY ONE FILTER HERE AT THE MOMENT
exp_filter = b+a*(exp_filt);

end