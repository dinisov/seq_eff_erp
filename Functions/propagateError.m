function ERROR = propagateError(FLIES, chosenFlies)

semAmplSEs = zeros(16,length(chosenFlies));
semPosAmplSEs = zeros(16,length(chosenFlies));
semNegAmplSEs = zeros(16,length(chosenFlies));

nERPsFly = zeros(1,length(chosenFlies));

%collect errors in matrix
for fly = 1:length(chosenFlies)        

        %number of ERPs for each fly (summing over sequences)
        nERPsFly(fly) = sum(FLIES(chosenFlies(fly)).nERPs);
        
        %standard errors for each fly
        semAmplSEs(:,fly) = FLIES(chosenFlies(fly)).ERROR.amplitude;
        semPosAmplSEs(:,fly) = FLIES(chosenFlies(fly)).ERROR.positiveAmplitude;
        semNegAmplSEs(:,fly) = FLIES(chosenFlies(fly)).ERROR.negativeAmplitude; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%error propagated according to
%sem(a1*mean_1 + .. + an*mean_n) = sqrt(a1^2*sem_n + ... + an^2*sem_n)
%for a weigthed mean the a,...,an are nERP1/(nERP1 + ... + nERPn) so
%sqrt((nERP1/(nERP1 + ... + nERPn))^2*sem_1 + ... + (nERP1/(nERP1 + ... + nERPn))^2*sem_n)
%or equally
%%sqrt((1/(nERP1 + ... + nERPn))^2*(nERP1^2*sem_1 + ... + nERP1^2*sem_n))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate errors by propagating errors from individual flies
for fly = 1:length(chosenFlies)
        % part of the error propagation calculation (nERP1^2*sem_1 + ... + nERP1^2*sem_n)
        semAmplSEs(:,fly) = semAmplSEs(:,fly).^2 * nERPsFly(fly).^2;
        semPosAmplSEs(:,fly) = semPosAmplSEs(:,fly).^2 * nERPsFly(fly).^2;
        semNegAmplSEs(:,fly) = semNegAmplSEs(:,fly).^2 * nERPsFly(fly).^2;
end

%finish calculating error propagation
%sqrt((1/(nERP1 + ... + nERPn))^2*(nERP1^2*sem_1 + ... + nERP1^2*sem_n))
semAmplSEs = sqrt(sum(semAmplSEs/(sum(nERPsFly)^2),2));
semPosAmplSEs = sqrt(sum(semPosAmplSEs/(sum(nERPsFly)^2),2));
semNegAmplSEs = sqrt(sum(semNegAmplSEs/(sum(nERPsFly)^2),2));

%propagated errors
ERROR.amplitude = semAmplSEs.';
ERROR.positiveAmplitude = semPosAmplSEs.';
ERROR.negativeAmplitude = semNegAmplSEs.';
ERROR.latencyToPeak = [];
ERROR.latencyToTrough = [];