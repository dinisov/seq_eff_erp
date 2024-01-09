% function [exp_filter] = exp_filter(b,a,alpha)

all_seq =  [1 1 1 1 1
            0 1 1 1 1
            0 0 1 1 1
            1 0 1 1 1
            0 0 0 1 1
            1 0 0 1 1
            1 1 0 1 1
            0 1 0 1 1
            0 0 0 0 1
            1 0 0 0 1
            1 1 0 0 1
            0 1 0 0 1
            1 1 1 0 1
            0 1 1 0 1
            0 0 1 0 1
            1 0 1 0 1];
        
% all_seq = ~all_seq;

%% filter

alpha = 0.7;

exp_filt = zeros(16,1);

for s = 1:16
    x_0 = 0.5;
    for i = 1:5
       x_0 = alpha*x_0 + (1-alpha)*all_seq(s,i); 
    end
    exp_filt(s) = x_0;
end

create_seq_eff_plot(-exp_filt,[]);

% exponents = repmat([4 3 2 1],16,1);
% 
% exp_filt = all_seq(:,1:4).*(alpha.^exponents);
% 
% exp_filt = sum(exp_filt,2)/sum(alpha.^(1:4));
% 
% exp_filt(all_seq(:,5) == 0) = 1-exp_filt(all_seq(:,5) == 0);
% 
% %% alternation/repetition filter
% 
% all_seq = 1-abs(diff(all_seq,1,2));
% 
% exponents = repmat([3 2 1],16,1);
% 
% exp_filt_ar = all_seq(:,1:3).*(alpha.^exponents);
% 
% exp_filt_ar = sum(exp_filt_ar,2)/sum(alpha.^(1:3));
% 
% exp_filt_ar(all_seq(:,4) == 0) = 1-exp_filt_ar(all_seq(:,4) == 0);

%figure; create_seq_eff_plot(1-(exp_filt_ar+exp_filt));

%CAREFUL ONLY ONE FILTER HERE AT THE MOMENT
% exp_filter = b+a*(exp_filt);
