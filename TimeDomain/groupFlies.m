function ALLFLIES = groupFlies(FLIES, chosenFlies)

ALLFLIES = struct;

amplitudeSEs = zeros(16,length(chosenFlies));
negativeAmplitudeSEs = zeros(16,length(chosenFlies));
positiveAmplitudeSEs = zeros(16,length(chosenFlies));
latencyToPeakSEs = zeros(16,length(chosenFlies));
latencyToTroughSEs = zeros(16,length(chosenFlies));
nERPsFly = zeros(1,length(chosenFlies));

semAmplSEs = zeros(16,length(chosenFlies));
semPosAmplSEs = zeros(16,length(chosenFlies));
semNegAmplSEs = zeros(16,length(chosenFlies));

for fly = 1:length(chosenFlies)

        %sequential effects results
        amplitudeSEs(:,fly) = (FLIES(chosenFlies(fly)).amplitudeSEs.');
        negativeAmplitudeSEs(:,fly) = (FLIES(chosenFlies(fly)).negativeAmplitudeSEs.');
        positiveAmplitudeSEs(:,fly) = (FLIES(chosenFlies(fly)).positiveAmplitudeSEs.');
        latencyToPeakSEs(:,fly) = (FLIES(chosenFlies(fly)).latencyToPeakSEs.');
        latencyToTroughSEs(:,fly) = (FLIES(chosenFlies(fly)).latencyToTroughSEs.');

        %standard errors for each fly
        semAmplSEs(:,fly) = FLIES(chosenFlies(fly)).semAmplSEs;
        semPosAmplSEs(:,fly) = FLIES(chosenFlies(fly)).semPosAmplSEs;
        semNegAmplSEs(:,fly) = FLIES(chosenFlies(fly)).semNegAmplSEs; 

        nERPsFly(fly) = sum(FLIES(chosenFlies(fly)).nERPs);

        % SEs weighted by how many ERPs were in each fly
        amplitudeSEs(:,fly) = amplitudeSEs(:,fly)*nERPsFly(fly);
        positiveAmplitudeSEs(:,fly) = positiveAmplitudeSEs(:,fly)*nERPsFly(fly);
        negativeAmplitudeSEs(:,fly) = negativeAmplitudeSEs(:,fly)*nERPsFly(fly);
        latencyToPeakSEs(:,fly) = latencyToPeakSEs(:,fly)*nERPsFly(fly);
        latencyToTroughSEs(:,fly) = latencyToTroughSEs(:,fly)*nERPsFly(fly);

        % part of the error propagation calculation (n_A^2*sem_A^2)
        semAmplSEs(:,fly) = semAmplSEs(:,fly).^2 * nERPsFly(fly).^2;
        semPosAmplSEs(:,fly) = semPosAmplSEs(:,fly).^2 * nERPsFly(fly).^2;
        semNegAmplSEs(:,fly) = semNegAmplSEs(:,fly).^2 * nERPsFly(fly).^2;

        %finish calculating error propagation
        % sem_{(n_A*A + n_B*B)/(n_A+n_B)^2} = sqrt(n_A^2/(n_A+n_B)^2 sem_A^2 + n_B^2/(n_A+n_B)^2 sem_B^2)
        semAmplSEs = sqrt(sum(semAmplSEs/(sum(nERPsFly)^2),2));
        semPosAmplSEs = sqrt(sum(semPosAmplSEs/(sum(nERPsFly)^2),2));
        semNegAmplSEs = sqrt(sum(semNegAmplSEs/(sum(nERPsFly)^2),2));
            
end

%divide SE profiles by total number of ERPs
amplitudeSEs = sum(amplitudeSEs,2)/sum(nERPsFly);
positiveAmplitudeSEs = sum(positiveAmplitudeSEs,2)/sum(nERPsFly);
negativeAmplitudeSEs = sum(negativeAmplitudeSEs,2)/sum(nERPsFly);
latencyToPeakSEs = sum(latencyToPeakSEs,2)/sum(nERPsFly);
latencyToTroughSEs = sum(latencyToTroughSEs,2)/sum(nERPsFly);

% add matrices with all profiles to ALLFLIES structure (ignores LIT vs DARK)
ALLFLIES.amplitudeSEs = amplitudeSEs.';
ALLFLIES.positiveAmplitudeSEs = positiveAmplitudeSEs.';           
ALLFLIES.negativeAmplitudeSEs = negativeAmplitudeSEs.';
ALLFLIES.latencyToPeakSEs = latencyToPeakSEs.';
ALLFLIES.latencyToTroughSEs = latencyToTroughSEs.';
ALLFLIES.nERPsFly = nERPsFly.';

%propagated errors
ALLFLIES.semAmplSEs = semAmplSEs.';
ALLFLIES.semPosAmplSEs = semPosAmplSEs.';
ALLFLIES.semNegAmplSEs = semNegAmplSEs.';