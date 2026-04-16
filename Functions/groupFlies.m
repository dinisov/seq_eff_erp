function ALLFLIES = groupFlies(FLIES, chosenFlies,errorMethod, reOrder, n_back,options)

arguments
    FLIES struct
    chosenFlies double
    errorMethod double
    reOrder double
    n_back double
    options.suppressANOVA double = 0
end

%Extra
suppressANOVA = options.suppressANOVA;

ALLFLIES = struct;

%matrices to hold SE profiles for each fly
%amplitudeSEs = zeros(16,length(chosenFlies));
%negativeAmplitudeSEs = zeros(16,length(chosenFlies));
%positiveAmplitudeSEs = zeros(16,length(chosenFlies));
%latencyToPeakSEs = zeros(16,length(chosenFlies));
%latencyToTroughSEs = zeros(16,length(chosenFlies));
%transectSEs = zeros(16,length(chosenFlies)); %Borrow above structure
%avTransectWindowSEs = zeros(16,length(chosenFlies)); %Borrow above structure
amplitudeSEs = zeros(0.5*2^n_back,length(chosenFlies));
negativeAmplitudeSEs = zeros(0.5*2^n_back,length(chosenFlies));
positiveAmplitudeSEs = zeros(0.5*2^n_back,length(chosenFlies));
latencyToPeakSEs = zeros(0.5*2^n_back,length(chosenFlies));
latencyToTroughSEs = zeros(0.5*2^n_back,length(chosenFlies));
transectSEs = zeros(0.5*2^n_back,length(chosenFlies)); %Borrow above structure
avTransectWindowSEs = zeros(0.5*2^n_back,length(chosenFlies)); %Borrow above structure
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
        transectSEs(:,fly) = (FLIES(chosenFlies(fly)).PROFILE.transect.')*nERPsFly(fly);
        avTransectWindowSEs(:,fly) = (FLIES(chosenFlies(fly)).PROFILE.avTransectWindow.')*nERPsFly(fly);
end

%preserve matrices with all profiles
allProfilesAmplitude = amplitudeSEs;
allProfilesPositiveAmplitude = positiveAmplitudeSEs;
allProfilesNegativeAmplitude = negativeAmplitudeSEs;
allProfilesLatencyToPeak= latencyToPeakSEs;
allProfilesLatencyToTrough = latencyToTroughSEs;
allProfilesTransect = transectSEs;
allProfilesAvTransectWindow = avTransectWindowSEs;

%divide SE profiles by total number of ERPs to finish weighted average
amplitudeSEs = sum(amplitudeSEs,2)/sum(nERPsFly);
positiveAmplitudeSEs = sum(positiveAmplitudeSEs,2)/sum(nERPsFly);
negativeAmplitudeSEs = sum(negativeAmplitudeSEs,2)/sum(nERPsFly);
latencyToPeakSEs = sum(latencyToPeakSEs,2)/sum(nERPsFly);
latencyToTroughSEs = sum(latencyToTroughSEs,2)/sum(nERPsFly);
transectSEs = sum(transectSEs,2)/sum(nERPsFly);
avTransectWindowSEs = sum(avTransectWindowSEs,2)/sum(nERPsFly);

%propagate errors
if errorMethod
    ALLFLIES.ERROR = propagateError(FLIES, chosenFlies);
else
    %ALLFLIES.ERROR = calculateError(FLIES, chosenFlies);
    ALLFLIES.ERROR = calculateError(FLIES, chosenFlies, reOrder, n_back,'suppressANOVA',suppressANOVA);
end

% profiles for all flies in group (ignores LIT vs DARK)
ALLFLIES.PROFILE.amplitude = amplitudeSEs.';
ALLFLIES.PROFILE.positiveAmplitude = positiveAmplitudeSEs.';           
ALLFLIES.PROFILE.negativeAmplitude = negativeAmplitudeSEs.';
ALLFLIES.PROFILE.latencyToPeak = latencyToPeakSEs.';
ALLFLIES.PROFILE.latencyToTrough = latencyToTroughSEs.';
ALLFLIES.PROFILE.transect = transectSEs.';
ALLFLIES.PROFILE.avTransectWindow = avTransectWindowSEs.';

ALLFLIES.nERPsFly = nERPsFly.';

ALLFLIES.allProfilesAmplitude = allProfilesAmplitude;
ALLFLIES.allProfilesPositiveAmplitude = allProfilesPositiveAmplitude;
ALLFLIES.allProfilesNegativeAmplitude = allProfilesNegativeAmplitude;
ALLFLIES.allProfilesLatencyToPeak = allProfilesLatencyToPeak;
ALLFLIES.allProfilesLatencyToTrough = allProfilesLatencyToTrough;
ALLFLIES.allProfilesTransect = allProfilesTransect;
ALLFLIES.allProfilesAvTransectWindow = allProfilesAvTransectWindow;