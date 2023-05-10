% analyse sequential dependencies in fly ERPs
close all;
clear;

%% load auxiliary functions
addpath('../');
addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');
addpath('../Functions');

load slrp_lrpr

load results

% these are the standard ALT/REP components from the literature (Jentzsch 2002)
lrpr = normalize(-lrpr); slrp = normalize(slrp); weird = normalize(weird);

%% type of analysis

% 1 - regular SEs; 2 - block experiments
analysisType = 1;

% if analysisType is 2, choose which stimulus to look in the train (starting at 5 up to train length)
focusPeak = 5;

%% toggles

% fit model
fitModel = true;

%which SEs to plot (in order: amplitude, pos amplitude, neg amplitude, latency to peak, latency to trough)
plotSelector = [1 0 0 0 0];

% whether to plot auxiliary plots (some are always plotted)
aux_plots = 1;

%% load data
homeDirectory = '../..';

resultsDirectory = [homeDirectory '\Results\12dot5Hz'];

struct_name = 'freq12dot5hz';

fly_record = readtable([homeDirectory '\Fly record\fly_record']);

%% restrict to some frequency
% fly_record = fly_record(fly_record.Frequency == 6.25,:);

%% restrictions on data of interest; always check this before running

% remove DARK flies
fly_record = fly_record(convertCharsToStrings(fly_record.Condition) == 'LIT',:);

% remove flies to be excluded (usually because data is unsound for some obvious reason)
fly_record = fly_record(~logical(fly_record.Exclude),:);

%remove jittering flies
fly_record = fly_record(~contains(fly_record.Comments,'jitter','IgnoreCase',true),:);

%remove red light flies
fly_record = fly_record(~contains(fly_record.Comments,'red','IgnoreCase',true),:);

%remove block paradigm flies
fly_record = fly_record(~contains(fly_record.Comments,'block','IgnoreCase',true),:);

%% 
% fly_record = fly_record(contains(fly_record.Comments,'adjfront','IgnoreCase',true),:);

%% choose flies and experiments
whichFly =      fly_record.Fly.';
flySet = unique(whichFly);

% choose which flies to run here
chosenFlies = [64];
% chosenFlies = setdiff(flySet, [24 25]);
% chosenFlies = flySet; % choose all flies
% chosenFlies = setdiff(chosenFlies, 24:29);

% choose which blocks to run
%NOTE: while unlikely as a request, this does not handle the case where two
%flies have a block with the same number but we would like to look at both
%flies but not one of the blocks with the same number
chosenBlocks = [4];
% chosenBlocks = unique(fly_record.Block.');% do not choose specific blocks

chosenOnes = ismember(fly_record.Block.', chosenBlocks) & ismember(fly_record.Fly.', chosenFlies);

%% experimental parameters in fly_record

% inter-stimulus interval (ISI in fly_record is the stimulus-off period, 
%i.e. the time between the end of the stimulus and the beginning of the
%next one
ISI = fly_record.ISI + fly_record.SDT; 
SDT = fly_record.SDT; % stimulus duration time

% relevant LFP channel
LFPChannel = fly_record.LFPChannel;

% whether to do a 1 or 2 photodiode analysis
% 2 is better currently but 1 is more universal/convenient
PHOTType = fly_record.PHOTType;

% for block experiments
interBlockPeriod = fly_record.InterBlockPeriod;
blockLength = fly_record.BlockLength;

%this is the window to look at around each peak
time_before_peak = fly_record.SDT.*fly_record.Window1;
time_after_peak = fly_record.ISI.*fly_record.Window2;

%peak detection threshold
peak_threshold = fly_record.Threshold;

%light_on_dark = 1 means a bright bar over a dark background was used
light_on_dark = strcmp(fly_record.Condition,'LIT').';
% chosenOnes = 1-light_on_dark; % choose all lit flies

%number of sd to trim from LFP
LFPsd = fly_record.LFPsd;

%% add data to structure according to block
% the structure may be largely empty if analysing only one fly
% the index here is that of the fly_record table, *not* the original block number

BLOCKS = struct;

% loop though all blocks, irrespective of experiment
for b = find(chosenOnes)
    
    date = datestr(fly_record.Date(b),'ddmmyy');
    block = num2str(fly_record.Block(b));

    % load this block's data
    load([homeDirectory '/SEoutput/' date '/LFP/Analyzed_TagTrials_block' block '/' date '_chunk_0']);
    
    % photodiode and lfp data
    LFP = EEG.LFP1.data(LFPChannel(b),:);
    PHOT = EEG.PHOT.data;
    rawPHOT = EEG.PHOT.data; %preserve raw unfiltered photodiode data
    resampleFreq = EEG.srate;
    
    % remove first/last n_sec seconds of recording
    n_sec = 1;
    chunk_a_time = n_sec*resampleFreq;
    LFP = LFP(:,chunk_a_time:end-chunk_a_time);
    PHOT = PHOT(:,chunk_a_time:end-chunk_a_time);
    rawPHOT = rawPHOT(:,chunk_a_time:end-chunk_a_time);
    
    %%
    
    % correct for the fact that the left photodiode is inverted
    % such that peaks are always upward for peak detection
    if light_on_dark(b)
        PHOT(2,:) = -PHOT(2,:);
        rawPHOT(2,:) = -rawPHOT(2,:);
    else
        PHOT(1,:) = -PHOT(1,:); 
        rawPHOT(1,:) = -rawPHOT(1,:);
    end

%     PHOT(3,:) = -PHOT(3,:);
    
    % butterworth filter LFP
    [b_f,a_f] = butter(6,100/resampleFreq*2);
    LFP = filter(b_f,a_f,LFP.').';

    wo = 50/(resampleFreq/2);  
    bw = wo/10;
    [b_f,a_f] = iirnotch(wo,bw);

    LFP = filter(b_f,a_f,LFP.').';
        
    LFP = smoothdata(LFP,'sgolay');

    PHOT = trim_phot_outliers(PHOT, 12);
    
    %% remove anything from LFP beyond some sd
    
    figure; plot(LFP); hold on;
    
    if ~isnan(LFPsd(b))
    
        % do a nice plot
        n_sd_lfp = LFPsd(b);
        sd_lfp = std(LFP); mean_lfp = mean(LFP);
        x_lim = xlim;
        plot([x_lim(1) x_lim(2)],[mean_lfp-n_sd_lfp*sd_lfp mean_lfp-n_sd_lfp*sd_lfp],'r');
        plot([x_lim(1) x_lim(2)],[mean_lfp+n_sd_lfp*sd_lfp mean_lfp+n_sd_lfp*sd_lfp],'r');

        % nuke everything beyond +/- the set number of sd
        LFP(LFP > mean_lfp+n_sd_lfp*sd_lfp | LFP < mean_lfp-n_sd_lfp*sd_lfp) = NaN;

    end
    
    %% 
    
%     figure; plot(rawPHOT(2,:)); hold on; plot(PHOT(2,:));

%     figure; plot(PHOT(1,:)); hold on; plot(xlim, [photSD photSD]); plot(xlim, [-photSD -photSD]);
%     figure; plot(PHOT(2,:)); hold on; plot(xlim, [photSD photSD]); plot(xlim, [-photSD -photSD]);

    BLOCKS(b).LFP = LFP;
    BLOCKS(b).PHOT = PHOT;
    BLOCKS(b).rawPHOT = rawPHOT;
    
    BLOCKS(b).times = EEG.times;
    BLOCKS(b).resampleFreq = resampleFreq;
    
    % only for block experiments
    BLOCKS(b).interBlockPeriod = interBlockPeriod(b);
    BLOCKS(b).focusPeak = focusPeak;
    BLOCKS(b).blockLength = blockLength(b);
    
    BLOCKS(b).time_before_peak = time_before_peak(b);
    BLOCKS(b).time_after_peak = time_after_peak(b);
    BLOCKS(b).ISI = ISI(b);
    BLOCKS(b).SDT = SDT(b);
    BLOCKS(b).peakThreshold = peak_threshold(b);
    
    %type of photodiode analysis
    BLOCKS(b).PHOTType = PHOTType(b);

end

%% if some block requires specific parameters, set them here
% careful index here is that of the fly_record table *not* the block number

%% analyse data per fly and whether experiment is lit or unlit

FLIES = struct;

% fit options (choices here are historical but keeping them for now)
options = optimset('Algorithm','interior-point','FinDiffType','central');

n_fly = 1;

for fly = chosenFlies
    
    lit_dark = {'DARK','LIT'};
   
   for lit = [0 1]
       
       % index of the blocks belonging to this fly if they are to be
       % included in the analysis
       thisFlyBlocks = BLOCKS(whichFly == fly & light_on_dark == lit & chosenOnes);
       
       % in case there are no LIT/DARK blocks for this fly, or they were
       % not selected
       if ~isempty(thisFlyBlocks)

           disp(['Processing fly #' num2str(fly)]);
           
           switch analysisType
               case 1

                   if thisFlyBlocks(1).PHOTType == 1
                    R = processBlocksOneChannel(thisFlyBlocks, aux_plots,'time');
                   elseif thisFlyBlocks(1).PHOTType == 2
                    R = processBlocks(thisFlyBlocks, aux_plots,'time');  
                   end
                   
               case 2
                    R = processBlocksBlocksExp(thisFlyBlocks, aux_plots);
           end
                    
           % sequential effects results
           FLIES(fly).(lit_dark{lit+1}).amplitudeSEs = (R.amplitudeSEs);
           FLIES(fly).(lit_dark{lit+1}).negativeAmplitudeSEs = R.negativeAmplitudeSEs;
           FLIES(fly).(lit_dark{lit+1}).latencyToPeakSEs = R.latencyToPeakSEs;
           FLIES(fly).(lit_dark{lit+1}).latencyToTroughSEs = R.latencyToTroughSEs;
           FLIES(fly).(lit_dark{lit+1}).positiveAmplitudeSEs = R.positiveAmplitudeSEs;
           
           %SEM for the different cases (except latency)
           FLIES(fly).(lit_dark{lit+1}).semAmplSEs = R.semAmplSEs;
           FLIES(fly).(lit_dark{lit+1}).semPosAmplSEs = R.semPosAmplSEs;
           FLIES(fly).(lit_dark{lit+1}).semNegAmplSEs = R.semNegAmplSEs;
           
           % mean ERPs, PHOT and number of ERPs for each sequence
           FLIES(fly).(lit_dark{lit+1}).nERPs = R.nERPs;
           FLIES(fly).(lit_dark{lit+1}).meanERPs = R.meanERPs;
           FLIES(fly).(lit_dark{lit+1}).meanPHOTs = R.meanPHOTs;

           % perform fits to a simple exponential filter
           if fitModel
               [x,~] = fmincon(@(x) least_squares_exp_filters(x(1),x(2),x(3),R.amplitudeSEs.'),[0 1 0.5],[],[],[],[],[-inf  -inf 0],[inf inf 1],[],options); %#ok<UNRCH>
               [x,~] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,(R.amplitudeSEs).'),[1 1 0 0],[],[],[],[],[-inf  -inf 0 -inf],[inf inf 0 inf],[],options);

               FLIES(fly).(lit_dark{lit+1}).FITS.model_fit_amplitude = x(4) + x(1)*slrp + x(2)*lrpr + x(3)*weird;
               FLIES(fly).(lit_dark{lit+1}).FITS.fit_params = x;

               results.(struct_name)(n_fly).model_fit = FLIES(fly).(lit_dark{lit+1}).FITS.model_fit_amplitude;
               results.(struct_name)(n_fly).fit_params = x;
           end
           
           % dump results in a structure to be saved in the end (may need to improve this)
           results.(struct_name)(n_fly).seq_eff_profile = (R.amplitudeSEs).';

       end
       
   end
 
   n_fly = n_fly + 1;
   
end

% save('results','results');

%% plot sequential dependencies per fly

lit_dark = {'DARK','LIT'};

for lit = [0 1]

    seq_eff_amplitude = zeros(16,length(chosenFlies));
    seq_eff_negative_amplitude = zeros(16,length(chosenFlies));
    seq_eff_positive_amplitude = zeros(16,length(chosenFlies));
    seq_eff_latency = zeros(16,length(chosenFlies));
    n = zeros(1,length(chosenFlies));
    
    for fly = chosenFlies
    
        % if DARK/LIT field is not empty for this fly
        if isfield(FLIES(fly),lit_dark{lit+1}) && ~isempty(FLIES(fly).(lit_dark{lit+1}))
            
            if plotSelector(1)
                %amplitude sequential effects
                figure('Name',['Amplitude_fly_' num2str(fly) '_' lit_dark{lit+1}],'NumberTitle','off');
                if fitModel
                    create_seq_eff_plot(FLIES(fly).(lit_dark{lit+1}).amplitudeSEs.',FLIES(fly).(lit_dark{lit+1}).FITS.model_fit_amplitude,'errors',FLIES(fly).(lit_dark{lit+1}).semAmplSEs.','scores',FLIES(fly).(lit_dark{lit+1}).FITS.fit_params(1:3)); %#ok<UNRCH>
                else
                      create_seq_eff_plot(FLIES(fly).(lit_dark{lit+1}).amplitudeSEs.',[],'errors',FLIES(fly).(lit_dark{lit+1}).semAmplSEs.');               %#ok<UNRCH>
                end
                    saveas(gcf,[resultsDirectory '/Amplitude/fly_' num2str(fly) '_' lit_dark{lit+1} '.png']);
    %             close(gcf);
            end
            
            if plotSelector(2)
                %positive amplitude SEs
                figure('Name',['Positive_amplitude_fly_' num2str(fly) '_' lit_dark{lit+1}],'NumberTitle','off');
                create_seq_eff_plot(FLIES(fly).(lit_dark{lit+1}).positiveAmplitudeSEs.',[],'errors',FLIES(fly).(lit_dark{lit+1}).semPosAmplSEs.');
                saveas(gcf,[resultsDirectory '/Positive amplitude/fly_' num2str(fly) '_' lit_dark{lit+1} '_positive_amplitude.png']);
    %             close(gcf);
            end
            
            if plotSelector(3)
                %negative amplitude SEs
                figure('Name',['Negative_amplitude_fly_' num2str(fly) '_' lit_dark{lit+1}],'NumberTitle','off')
                create_seq_eff_plot(FLIES(fly).(lit_dark{lit+1}).negativeAmplitudeSEs.',[],'errors',FLIES(fly).(lit_dark{lit+1}).semNegAmplSEs.');
                saveas(gcf,[resultsDirectory '/Negative amplitude/fly_' num2str(fly) '_' lit_dark{lit+1} '_negative_amplitude.png']);
    %             close(gcf);
            end
    
            if plotSelector(4)
                %latency to peak sequential effects
                figure('Name',['Latency_to_peak_fly_' num2str(fly) '_' lit_dark{lit+1}],'NumberTitle','off');
                create_seq_eff_plot(FLIES(fly).(lit_dark{lit+1}).latencyToPeakSEs.',[]);
                saveas(gcf,[resultsDirectory '/Latency/fly_' num2str(fly) '_' lit_dark{lit+1} '_latency_to_peak.png']);
    %             close(gcf);
            end
            
            if plotSelector(5)
                %latency to trough sequential effects
                figure('Name',['Latency_to_trough_fly_' num2str(fly) '_' lit_dark{lit+1}],'NumberTitle','off');
                create_seq_eff_plot(FLIES(fly).(lit_dark{lit+1}).latencyToTroughSEs.',[]);
                saveas(gcf,[resultsDirectory '/Latency/fly_' num2str(fly) '_' lit_dark{lit+1} '_latency_to_trough.png']);
    %             close(gcf);
            end

%             close all

        end
        
    end
    
end

%% calculate SEs for all flies by averaging SE profiles

if length(chosenFlies) > 1

    lit_dark = {'DARK','LIT'};

    for lit = [0 1]

        % if DARK/LIT field is not empty
        if isfield(FLIES,lit_dark{lit+1})

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

                if ~isempty(FLIES(chosenFlies(fly)).(lit_dark{lit+1})) %&& isfield(FLIES(chosenFlies(fly)),lit_dark{lit+1})

                    %sequential effects results
                    amplitudeSEs(:,fly) = (FLIES(chosenFlies(fly)).(lit_dark{lit+1}).amplitudeSEs.');
                    negativeAmplitudeSEs(:,fly) = (FLIES(chosenFlies(fly)).(lit_dark{lit+1}).negativeAmplitudeSEs.');
                    positiveAmplitudeSEs(:,fly) = (FLIES(chosenFlies(fly)).(lit_dark{lit+1}).positiveAmplitudeSEs.');
                    latencyToPeakSEs(:,fly) = (FLIES(chosenFlies(fly)).(lit_dark{lit+1}).latencyToPeakSEs.');
                    latencyToTroughSEs(:,fly) = (FLIES(chosenFlies(fly)).(lit_dark{lit+1}).latencyToTroughSEs.');
    
                    %standard errors for each fly
                    semAmplSEs(:,fly) = FLIES(chosenFlies(fly)).(lit_dark{lit+1}).semAmplSEs;
                    semPosAmplSEs(:,fly) = FLIES(chosenFlies(fly)).(lit_dark{lit+1}).semPosAmplSEs;
                    semNegAmplSEs(:,fly) = FLIES(chosenFlies(fly)).(lit_dark{lit+1}).semNegAmplSEs; 
    
                    nERPsFly(fly) = sum(FLIES(chosenFlies(fly)).(lit_dark{lit+1}).nERPs);
    
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

            end
            
            [x,~] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,(sum(amplitudeSEs,2)/sum(nERPsFly))),[1 1 1 0],[],[],[],[],[-inf  -inf -inf -inf],[inf inf inf inf],[],options);
%             [x,~] = fmincon(@(x) least_squares_slrp_lrpr(x(1),x(2),x(3),slrp,lrpr,sum(amplitudeSEs,2)/sum(nERPsFly).'),[1 1 0],[],[],[],[],[-inf  -inf -inf],[inf inf inf],[],options);
           
            if plotSelector(1)
                figure('Name',['Amplitude_all_flies_method_2' lit_dark{lit+1}],'NumberTitle','off');
                best_fit = x(4) + x(1)*slrp + x(2)*lrpr + x(3)*weird;
                create_seq_eff_plot((sum(amplitudeSEs,2)/sum(nERPsFly)),best_fit,'errors',semAmplSEs,'scores',x(1:3));
                saveas(gcf,[resultsDirectory '/All flies 2/all_flies_' lit_dark{lit+1} '_amplitude.png']);
            end
            
            if plotSelector(2)
                figure('Name',['Positive_amplitude_all_flies_method_2' lit_dark{lit+1}],'NumberTitle','off');
                create_seq_eff_plot(sum(positiveAmplitudeSEs,2)/sum(nERPsFly),[],'errors',semPosAmplSEs);
                saveas(gcf,[resultsDirectory '/All flies 2/all_flies_' lit_dark{lit+1} '_positive_amplitude.png']);
            end
            
            if plotSelector(3)
                figure('Name',['Negative_amplitude_all_flies_method_2' lit_dark{lit+1}],'NumberTitle','off');
                create_seq_eff_plot(sum(negativeAmplitudeSEs,2)/sum(nERPsFly),[],'errors',semNegAmplSEs);
                saveas(gcf,[resultsDirectory '/All flies 2/all_flies_' lit_dark{lit+1} 'negative_amplitude.png']);
            end
            
            if plotSelector(4)
                figure('Name',['Latency_all_flies_method_2' lit_dark{lit+1}],'NumberTitle','off');
                create_seq_eff_plot(sum(latencyToPeakSEs,2)/sum(nERPsFly),[]);
                saveas(gcf,[resultsDirectory '/All flies 2/all_flies_' lit_dark{lit+1} '_latency.png']);
            end
            
            if plotSelector(5)
                figure('Name',['Latency_all_flies_method_2' lit_dark{lit+1}],'NumberTitle','off');
                create_seq_eff_plot(sum(latencyToTroughSEs,2)/sum(nERPsFly),[]);
                saveas(gcf,[resultsDirectory '/All flies 2/all_flies_' lit_dark{lit+1} '_latency.png']);
            end

        end

    end

end
        
%% calculate SEs for all flies by stacking all ERPs (this is probably not sound)

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
%             figure('Name',['Negative_amplitude_all_flies_' lit_dark{lit+1}],'NumberTitle','off');
%             create_seq_eff_plot(R.negativeAmplitudeSEs.',[],'errors',R.semNegAmplSEs.');
%             saveas(gcf,[resultsDirectory '/All flies/all_flies_' lit_dark{lit+1} '_negative_amplitude.png']);
% 
%             figure('Name',['Latency_all_flies_' lit_dark{lit+1}],'NumberTitle','off');
%             create_seq_eff_plot(R.latencySEs.',[]);
%             saveas(gcf,[resultsDirectory '/All flies/all_flies_' lit_dark{lit+1} '_latency.png']);
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