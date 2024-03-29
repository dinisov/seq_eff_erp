function ERROR = calculateError(FLIES, chosenFlies)

amplitudeSEs = zeros(16,length(chosenFlies));
positiveAmplitudeSEs = zeros(16,length(chosenFlies));
negativeAmplitudeSEs = zeros(16,length(chosenFlies));

%collect profiles in matrix
for fly = 1:length(chosenFlies)        
    %sequential effects results
    amplitudeSEs(:,fly) = FLIES(chosenFlies(fly)).PROFILE.amplitude.';
    positiveAmplitudeSEs(:,fly) = FLIES(chosenFlies(fly)).PROFILE.positiveAmplitude.';
    negativeAmplitudeSEs(:,fly) = FLIES(chosenFlies(fly)).PROFILE.negativeAmplitude.';
end

%subtract mean
amplitudeSEs = amplitudeSEs-mean(amplitudeSEs);
positiveAmplitudeSEs = positiveAmplitudeSEs-mean(positiveAmplitudeSEs);
negativeAmplitudeSEs = negativeAmplitudeSEs-mean(negativeAmplitudeSEs);

%onion plot of profiles
figure; create_seq_eff_plot(amplitudeSEs,[]);

%do ANOVA on subtracted mean data (for amplitude)
anova1(amplitudeSEs.');

semAmplSEs = std(amplitudeSEs,[],2)/sqrt(length(chosenFlies));
semPosAmplSEs = std(positiveAmplitudeSEs,[],2)/sqrt(length(chosenFlies));
semNegAmplSEs = std(negativeAmplitudeSEs,[],2)/sqrt(length(chosenFlies));

figure; create_seq_eff_plot(mean(amplitudeSEs,2),[],'errors',semAmplSEs);

%propagated errors
ERROR.amplitude = semAmplSEs.';
ERROR.positiveAmplitude = semPosAmplSEs.';
ERROR.negativeAmplitude = semNegAmplSEs.';
ERROR.latencyToPeak = [];
ERROR.latencyToTrough = [];