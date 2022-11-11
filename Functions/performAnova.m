function performAnova(allERPs, siz)

% full data for maxima and minima across all stimuli
aux = squeeze(reshape(allERPs,[1,size(allERPs,1)*size(allERPs,2),size(allERPs,3)])).';

% data tables for performing an ANOVA to check for the effect of
% sequence
dataMaxima = aux(:,sub2ind(siz,ind_max_erp,1:16));
dataMinima = aux(:,sub2ind(siz,ind_min_erp,1:16));
dataAmplitude = dataMaxima-dataMinima;

anova1(dataAmplitude);