%this function finds the predictive probabilities for a given data set depending on a set of parameters
%this is the function that is minimised when fitting the model to the data
%it returns a sum of least square errors between the data and the model
function [least2] = least_squares_exp_filters(b,a,alpha,data)

rt_prediction = exp_filter(b,a,alpha);

least2 = sum((rt_prediction-data).^2);
    
%disp(least2);

end