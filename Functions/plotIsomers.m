function plotIsomers(allERPs, isomer, window, n_back,resampleFreq)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    isomer1 = isomer; isomer2 = 1-isomer;
    
    index1 = zeros(1,(2^n_back)); index2 = zeros(1,(2^n_back));
    
    %convert to index
    for n = 5:length(isomer)
        
       idx1 = bin2dec(num2str(isomer1(n-n_back+1:n))) + 1;
       idx2 = bin2dec(num2str(isomer2(n-n_back+1:n))) + 1;
       
       index1(idx1) = idx1;
       index2(idx2) = idx2;
        
    end
    
    index1 = index1+fliplr(index1); index1(17:end) = []; index1 = index1(seq_eff_order(n_back));
    index2 = index2+fliplr(index2); index2(17:end) = []; index2 = index2(seq_eff_order(n_back));

    %uncomment to check if your isomer is sound (should be 1:32)
% 	disp(sort(union(index1,index2)));
    
    %grab all the ERPs for each isomer, already in SE standard order
    allERPs1 = allERPs(:,index1,:); 
    allERPs2 = allERPs(:,index2,:);
    
    R1 = calculateSEs(allERPs1,[],0,window, resampleFreq);
    R2 = calculateSEs(allERPs2,[],0,window, resampleFreq);
    
    figure; create_seq_eff_plot(R1.amplitudeSEs.',R2.amplitudeSEs.');
    h = legend({'First Isomer','Second Isomer'});
    set(h,'FontSize',6);

end

