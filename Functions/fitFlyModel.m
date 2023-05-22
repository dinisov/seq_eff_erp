function FIT = fitFlyModel(data)
%fitFlyModel fit a model to a fly (or flies') data
%   right now this fits a linear combination of lrpr, slrp and weird

    % fit options (choices here are historical but keeping them for now)
    options = optimset('Algorithm','interior-point','FinDiffType','central');

    load('slrp_lrpr','lrpr','slrp','weird');

    % these are the standard ALT/REP components from the literature (Jentzsch 2002)
    lrpr = normalize(-lrpr); slrp = normalize(slrp); weird = normalize(weird);

    % fit to a linear combination of SE components 
    [x,~] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,data.amplitudeSEs.'),[1 1 1 0],[],[],[],[],[-inf  -inf -inf -inf],[inf inf inf inf],[],options);

    FIT.model_fit_amplitude = x(1)*slrp + x(2)*lrpr + x(3)*weird + x(4);
    FIT.fit_params = x;
end

