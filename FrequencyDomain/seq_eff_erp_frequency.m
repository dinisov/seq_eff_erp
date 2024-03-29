% analyse sequential dependencies in fly ERPs
close all;
clear;

%% load auxiliary functions
addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');
addpath('./Functions');

%% load data
homeDirectory = 'D:\group_swinderen\Dinis';

resultsDirectory = [homeDirectory '\Results_Frequency'];

fly_record = readtable('fly_record');

%% restrict to some frequency
fly_record = fly_record(fly_record.Frequency == 1.25,:);

%% restrict to LIT or DARK 
fly_record = fly_record(convertCharsToStrings(fly_record.Condition) == 'LIT',:);

% remove flies to be excluded (usually because data is unsound for some obvious reason)
fly_record = fly_record(~logical(fly_record.Exclude),:);

%remove jittering flies
fly_record = fly_record(~contains(fly_record.Comments,'jittering','IgnoreCase',true),:);

%remove red light flies
fly_record = fly_record(~contains(fly_record.Comments,'red','IgnoreCase',true),:);

%% frequency bounds to calculate profiles for

frequency_bounds = [0 60];

%% choose flies and experiments
whichFly =      fly_record.Fly.';
flySet = unique(whichFly);

% choose which flies to run here
chosenFlies = [6];
% chosenFlies = flySet; % choose all flies

% choose which blocks to run
%NOTE: while unlikely as a request, this does not handle the case where two
%flies have a block with the same number but we would like to look at both
%flies but not one of the blocks with the same number
chosenBlocks = [8 10];
% chosenBlocks = unique(fly_record.Block.');% do not choose specific blocks

chosenOnes = ismember(fly_record.Block.', chosenBlocks) & ismember(fly_record.Fly.', chosenFlies);

%% experimental parameters and window calculation
ISI = fly_record.ISI + fly_record.SDT;

% stimulus duration time
SDT = fly_record.SDT;

%this is the window to look at around each peak
time_before_peak = zeros(size(ISI));
time_after_peak = fly_record.SDT;

%peak detection threshold
peak_threshold = fly_record.Threshold;

%light_on_dark = 1 means a bright bar over a dark background was used
light_on_dark = strcmp(fly_record.Condition,'LIT').';
% chosenOnes = 1-light_on_dark; % choose all lit flies

%% whether to plot auxiliary plots
aux_plots = 1;

%% add data to structure according to block
% the structure may be largely empty if analysing only one fly
% the index here is that of the fly_record table, *not* the original block number

BLOCKS = struct;

% loop though all blocks, irrespective of experiment
for b = find(chosenOnes)
    
    date = datestr(fly_record.Date(b),'ddmmyy');
    block = num2str(fly_record.Block(b));

    % load this block's data
    load([homeDirectory '/Output/' date '/LFP/Analyzed_TagTrials_block' block '/' date '_chunk_0']);
    
    % photodiode and lfp data
    LFP = EEG.LFP1.data(1,:);
    PHOT = EEG.PHOT.data;
    rawPHOT = EEG.PHOT.data; %preserve raw unfiltered photodiode data
    resampleFreq = EEG.srate;
    
    % correct for the fact that the left photodiode is inverted
    % such that peaks are always upward for peak detection
    if light_on_dark(b)
        PHOT(2,:) = -PHOT(2,:);
        rawPHOT(2,:) = -rawPHOT(2,:);
    else
        PHOT(1,:) = -PHOT(1,:); 
        rawPHOT(1,:) = -rawPHOT(1,:);
    end

    PHOT = trim_phot_outliers(PHOT, 5);

    BLOCKS(b).LFP = LFP;
    BLOCKS(b).PHOT = PHOT;
    BLOCKS(b).rawPHOT = rawPHOT;
    
    BLOCKS(b).times = EEG.times;
    BLOCKS(b).resampleFreq = resampleFreq;
    
    BLOCKS(b).time_before_peak = time_before_peak(b);
    BLOCKS(b).time_after_peak = time_after_peak(b);
    BLOCKS(b).ISI = ISI(b);
    BLOCKS(b).SDT = SDT(b);
    BLOCKS(b).peakThreshold = peak_threshold(b);

end

%% analyse data per fly and whether experiment is lit or unlit

FLIES = struct;

for fly = chosenFlies
    
    lit_dark = {'DARK','LIT'};
   
   for lit = [0 1]
       
       % index of the blocks belonging to this fly if they are to be
       % included in the analysis
       thisFlyBlocks = BLOCKS(whichFly == fly & light_on_dark == lit & chosenOnes);
       
       % in case there are no LIT/DARK blocks for this fly, or they were
       % not selected
       if ~isempty(thisFlyBlocks)
       
           R = processBlocks(thisFlyBlocks, aux_plots, 'frequency');

           % sequential effects results
           FLIES(fly).(lit_dark{lit+1}).magnitudeSEs = R.magnitudeSEs;
           FLIES(fly).(lit_dark{lit+1}).phaseSEs = R.phaseSEs;
           
           %SEM for the different cases (except latency)
%            FLIES(fly).(lit_dark{lit+1}).semAmplSEs = R.semAmplSEs;
%            FLIES(fly).(lit_dark{lit+1}).semPosAmplSEs = R.semPosAmplSEs;
%            FLIES(fly).(lit_dark{lit+1}).semNegAmplSEs = R.semNegAmplSEs;
           
           % mean ERPs, PHOT and number of ERPs for each sequence
           FLIES(fly).(lit_dark{lit+1}).nERPs = R.nERPs;
           FLIES(fly).(lit_dark{lit+1}).meanERPs = R.meanERPs;
           FLIES(fly).(lit_dark{lit+1}).meanPHOTs = R.meanPHOTs;
           
       end
       
   end
 
end

save('data_flies_frequency','FLIES');

%% plot sequential dependencies per fly

% lit_dark = {'DARK','LIT'};
% 
% for lit = [0 1]
% 
%     seq_eff_amplitude = zeros(16,length(chosenFlies));
%     seq_eff_negative_amplitude = zeros(16,length(chosenFlies));
%     n = zeros(1,length(chosenFlies));
%     
%     for fly = chosenFlies
%         
%         resultsDirectoryFly = [resultsDirectory '\Fly' num2str(fly)];
%         
%         if ~exist(resultsDirectoryFly, 'dir')
%             mkdir([resultsDirectoryFly '\Magnitude']);
%             mkdir([resultsDirectoryFly '\Phase']);
%         end
%         
%         delete([resultsDirectoryFly '\Magnitude\*']);
%         delete([resultsDirectoryFly '\Phase\*']);
%     
%         % if DARK/LIT field is not empty for this fly
%         if isfield(FLIES(fly),lit_dark{lit+1}) && ~isempty(FLIES(fly).(lit_dark{lit+1}))
%             
%             L = size(FLIES(fly).(lit_dark{lit+1}).magnitudeSEs,1);
%             
%             for coeff = 1:ceil(L/2)
%             
%                 freq = resampleFreq*(coeff-1)/L;
%                 
%                 if freq > frequency_bounds(1) && freq < frequency_bounds(2)
%                 
%                     %magnitude sequential effects
%                     figure('Name',['Magnitude_fly_' num2str(fly) '_' lit_dark{lit+1}],'NumberTitle','off');
%                     create_seq_eff_plot(FLIES(fly).(lit_dark{lit+1}).magnitudeSEs(coeff,:).',[]);%FLIES(fly).(lit_dark{lit+1}).amplitudeSEs
%                     saveas(gcf,[resultsDirectoryFly '\Magnitude\fly_' num2str(fly) '_' lit_dark{lit+1} '_magnitude_' num2str(freq) 'Hz.png']);
% 
%                     %phase sequential effects
%                     figure('Name',['Phase_fly_' num2str(fly) '_' lit_dark{lit+1}],'NumberTitle','off');
%                     create_seq_eff_plot(FLIES(fly).(lit_dark{lit+1}).phaseSEs(coeff,:).',[]);
%                     saveas(gcf,[resultsDirectoryFly '\Phase\fly_' num2str(fly) '_' lit_dark{lit+1} '_phase_' num2str(freq) 'Hz.png']);
% 
%                     close all;
%                 
%                 end
% 
%             end
%                 
%         end
%         
%     end
%     
% end

%% calculate SEs for all flies by averaging SE profiles

% if length(chosenFlies) > 1
% 
%     lit_dark = {'DARK','LIT'};
% 
%     for lit = [0 1]
% 
%         amplitudeSEs = zeros(16,length(chosenFlies));
%         negativeAmplitudeSEs = zeros(16,length(chosenFlies));
%         positiveAmplitudeSEs = zeros(16,length(chosenFlies));
%         latencySEs = zeros(16,length(chosenFlies));
%         nERPsFly = zeros(1,length(chosenFlies));
% 
%         semAmplSEs = zeros(16,length(chosenFlies));
%         semPosAmplSEs = zeros(16,length(chosenFlies));
%         semNegAmplSEs = zeros(16,length(chosenFlies));
% 
%         for fly = 1:length(chosenFlies)
% 
%             % if DARK/LIT field is not empty for this fly
%             if isfield(FLIES(chosenFlies(fly)),lit_dark{lit+1}) && ~isempty(FLIES(chosenFlies(fly)).(lit_dark{lit+1}))
% 
%                 %sequential effects results
%                 amplitudeSEs(:,fly) = FLIES(chosenFlies(fly)).(lit_dark{lit+1}).amplitudeSEs.';
%                 negativeAmplitudeSEs(:,fly) = FLIES(chosenFlies(fly)).(lit_dark{lit+1}).negativeAmplitudeSEs.';
%                 positiveAmplitudeSEs(:,fly) = FLIES(chosenFlies(fly)).(lit_dark{lit+1}).positiveAmplitudeSEs.';
%                 latencySEs(:,fly) = FLIES(chosenFlies(fly)).(lit_dark{lit+1}).latencySEs.';
% 
%                 %standard errors
%                 semAmplSEs(:,fly) = FLIES(chosenFlies(fly)).(lit_dark{lit+1}).semAmplSEs;
%                 semPosAmplSEs(:,fly) = FLIES(chosenFlies(fly)).(lit_dark{lit+1}).semPosAmplSEs;
%                 semNegAmplSEs(:,fly) = FLIES(chosenFlies(fly)).(lit_dark{lit+1}).semNegAmplSEs; 
% 
%                 nERPsFly(fly) = sum(FLIES(chosenFlies(fly)).(lit_dark{lit+1}).nERPs);
% 
%                 amplitudeSEs(:,fly) = amplitudeSEs(:,fly).*nERPsFly(:,fly);
%                 positiveAmplitudeSEs(:,fly) = positiveAmplitudeSEs(:,fly).*nERPsFly(:,fly);
%                 negativeAmplitudeSEs(:,fly) = negativeAmplitudeSEs(:,fly).*nERPsFly(:,fly);
%                 latencySEs(:,fly) = latencySEs(:,fly).*nERPsFly(:,fly);
% 
%                 % part of the error propagation calculation (n_A^2*sem_A^2)
%                 semAmplSEs(:,fly) = semAmplSEs(:,fly).^2 * nERPsFly(fly).^2;
%                 semPosAmplSEs(:,fly) = semPosAmplSEs(:,fly).^2 * nERPsFly(fly).^2;
%                 semNegAmplSEs(:,fly) = semNegAmplSEs(:,fly).^2 * nERPsFly(fly).^2;
% 
%             end
% 
%         end
% 
%         %finish calculating error propagation
%         % sem_{(n_A*A + n_B*B)/(n_A+n_B)^2} = sqrt(n_A^2/(n_A+n_B)^2 sem_A^2 + n_B^2/(n_A+n_B)^2 sem_B^2)
%         semAmplSEs = sqrt(sum(semAmplSEs/(sum(nERPsFly)^2),2));
%         semPosAmplSEs = sqrt(sum(semPosAmplSEs/(sum(nERPsFly)^2),2));
%         semNegAmplSEs = sqrt(sum(semNegAmplSEs/(sum(nERPsFly)^2),2));
% 
%         figure('Name',['Amplitude_all_flies_method_2' lit_dark{lit+1}],'NumberTitle','off');
%         create_seq_eff_plot(sum(amplitudeSEs,2)./sum(nERPsFly,2),[],'errors',semAmplSEs);
%         saveas(gcf,[resultsDirectory '/All flies 2/all_flies_' lit_dark{lit+1} '_amplitude.png']);
% 
%         figure('Name',['Positive_amplitude_all_flies_method_2' lit_dark{lit+1}],'NumberTitle','off');
%         create_seq_eff_plot(sum(positiveAmplitudeSEs,2)./sum(nERPsFly,2),[],'errors',semPosAmplSEs);
%         saveas(gcf,[resultsDirectory '/All flies 2/all_flies_' lit_dark{lit+1} '_positive_amplitude.png']);
% 
%         figure('Name',['Negative_amplitude_all_flies_method_2' lit_dark{lit+1}],'NumberTitle','off');
%         create_seq_eff_plot(sum(negativeAmplitudeSEs,2)./sum(nERPsFly,2),[],'errors',semNegAmplSEs);
%         saveas(gcf,[resultsDirectory '/All flies 2/all_flies_' lit_dark{lit+1} 'negative_amplitude.png']);
% 
%         figure('Name',['Latency_all_flies_method_2' lit_dark{lit+1}],'NumberTitle','off');
%         create_seq_eff_plot(sum(latencySEs,2)./sum(nERPsFly,2),[]);
%         saveas(gcf,[resultsDirectory '/All flies 2/all_flies_' lit_dark{lit+1} '_latency.png']);
% 
%     end
% 
% end
        
%% calculate SEs for all flies by stacking all ERPs

% if length(chosenFlies) > 1
% 
%     for lit = [0 1]
% 
%         lit_dark = {'DARK','LIT'};
% 
%         whichBlocks = BLOCKS(light_on_dark == lit & chosenOnes);
% 
%         if ~isempty(whichBlocks)
% 
%             R = processBlocks(whichBlocks, aux_plots);
% 
%             figure('Name',['Amplitude_all_flies_' lit_dark{lit+1}],'NumberTitle','off');
%             create_seq_eff_plot(R.amplitudeSEs.',[],'errors',R.semAmplSEs.');
%             saveas(gcf,[resultsDirectory '/All flies/all_flies_' lit_dark{lit+1} '_amplitude.png']);
% 
%             figure('Name',['Positive_amplitude_all_flies_' lit_dark{lit+1}],'NumberTitle','off');
%             create_seq_eff_plot(R.positiveAmplitudeSEs.',[],'errors',R.semPosAmplSEs.');
%             saveas(gcf,[resultsDirectory '/All flies/all_flies_' lit_dark{lit+1} '_positive_amplitude.png']);
% 
%         end
% 
%     end
% 
% end

%% plot mean ERP for lit vs unlit for each fly

% lit_dark = {'DARK','LIT'};
% 
% for fly = 1:length(chosenFlies)
%     
%     for lit = [0 1]
%     
%         figure('Name',['meanERP_fly_' num2str(chosenFlies(fly)) '_' lit_dark{lit+1} '.png'],'NumberTitle','off');
%         hold on
%         plot(normalize(sum(FLIES(chosenFlies(fly)).(lit_dark{lit+1}).meanERPS,2)/16));
%         plot(normalize(sum(FLIES(chosenFlies(fly)).(lit_dark{lit+1}).meanPHOT,2)/16));
%         legend(['ERP ' (lit_dark{lit+1})],['PHOT ' (lit_dark{lit+1})]);
%         saveas(gcf,['meanERP_fly_' num2str(chosenFlies(fly)) '_' lit_dark{lit+1} '.png']);
%     
%     end
%     
% end