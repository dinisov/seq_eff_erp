function plotIsomers(allERPs, isomer, n_back)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    isomer1 = isomer; isomer2 = 1-isomer;
    
    index1 = zeros(1,(2^n_back)/2); index2 = zeros(1,(2^n_back)/2);
    
    %convert to index
    for n = 5:length(isomer)
       
       index1(n-n_back+1) = bin2dec(num2str(isomer1(n-n_back+1:n))) + 1;
       index2(n-n_back+1) = bin2dec(num2str(isomer2(n-n_back+1:n))) + 1;
        
    end
    
    % uncomment to check if isomers are sound
%     disp(sort(union(index1,index2)));

    allERPs1 = allERPs(:,index1,:); allERPs2 = allERPs(:,index2,:);


end

