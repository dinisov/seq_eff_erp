function ALLFLIES = groupFlies(FLIES, chosenFlies,errorMethod)

ALLFLIES = struct;

%matrices to hold SE profiles for each fly
amplitudeSEs = zeros(16,length(chosenFlies));
negativeAmplitudeSEs = zeros(16,length(chosenFlies));
positiveAmplitudeSEs = zeros(16,length(chosenFlies));
latencyToPeakSEs = zeros(16,length(chosenFlies));
latencyToTroughSEs = zeros(16,length(chosenFlies));
nERPsFly = zeros(1,length(chosenFlies));

%collect profiles and errors in matrix
for fly = 1:length(chosenFlies)
        %number of ERPs for each fly (summing over sequences)
        nERPsFly(fly) = sum(FLIES(chosenFlies(fly)).nERPs);
        
        %sequential effects results
        amplitudeSEs(:,fly) = (FLIES(chosenFlies(fly)).PROFILE.amplitude.')*nERPsFly(fly);
        positiveAmplitudeSEs(:,fly) = (FLIES(chosenFlies(fly)).PROFILE.positiveAmplitude.')*nERPsFly(fly);
        negativeAmplitudeSEs(:,fly) = (FLIES(chosenFlies(fly)).PROFILE.negativeAmplitude.')*nERPsFly(fly);
        latencyToPeakSEs(:,fly) = (FLIES(chosenFlies(fly)).PROFILE.latencyToPeak.')*nERPsFly(fly);
        latencyToTroughSEs(:,fly) = (FLIES(chosenFlies(fly)).PROFILE.latencyToTrough.')*nERPsFly(fly);
end

%preserve matrices with all profiles
allProfilesAmplitude = amplitudeSEs;
allProfilesPositiveAmplitude = positiveAmplitudeSEs;
allProfilesNegativeAmplitude = negativeAmplitudeSEs;
allProfilesLatencyToPeak= latencyToPeakSEs;
allProfilesLatencyToTrough = latencyToTroughSEs;

%divide SE profiles by total number of ERPs to finish weighted average
amplitudeSEs = sum(amplitudeSEs,2)/sum(nERPsFly);
positiveAmplitudeSEs = sum(positiveAmplitudeSEs,2)/sum(nERPsFly);
negativeAmplitudeSEs = sum(negativeAmplitudeSEs,2)/sum(nERPsFly);
latencyToPeakSEs = sum(latencyToPeakSEs,2)/sum(nERPsFly);
latencyToTroughSEs = sum(latencyToTroughSEs,2)/sum(nERPsFly);

%propagate errors
if errorMethod
    ALLFLIES.ERROR = propagateError(FLIES, chosenFlies);
else
    ALLFLIES.ERROR = calculateError(FLIES, chosenFlies);
end

% profiles for all flies in group (ignores LIT vs DARK)
ALLFLIES.PROFILE.amplitude = amplitudeSEs.';
ALLFLIES.PROFILE.positiveAmplitude = positiveAmplitudeSEs.';           
ALLFLIES.PROFILE.negativeAmplitude = negativeAmplitudeSEs.';
ALLFLIES.PROFILE.latencyToPeak = latencyToPeakSEs.';
ALLFLIES.PROFILE.latencyToTrough = latencyToTroughSEs.';

ALLFLIES.nERPsFly = nERPsFly.';

ALLFLIES.allProfilesAmplitude = allProfilesAmplitude;
ALLFLIES.allProfilesPositiveAmplitude = allProfilesPositiveAmplitude;
ALLFLIES.allProfilesNegativeAmplitude = allProfilesNegativeAmplitude;
ALLFLIES.allProfilesLatencyToPeak = allProfilesLatencyToPeak;
ALLFLIES.allProfilesLatencyToTrough = allProfilesLatencyToTrough;