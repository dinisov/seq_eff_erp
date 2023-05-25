% analyse sequential dependencies in fly ERPs
close all;
clear;

%% load auxiliary functions
addpath('../');
addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');
addpath('../Functions');

% load slrp_lrpr

% load results

% these are the standard ALT/REP components from the literature (Jentzsch 2002)
% lrpr = normalize(-lrpr); slrp = normalize(slrp); weird = normalize(weird);

%% type of analysis

% 1 - regular SEs; 2 - block experiments
analysisType = 1;

% if analysisType is 2, choose which stimulus to look in the train (starting at 5 up to train length)
focusPeak = 5;

%% toggles

%error method (0 - a la Bruno and Matt; 1- propagation )
errorMethod = 0;

% fit model
fitModel = 1;

%plot individual flies?
plotIndividualFlies = 1;

%which SEs to plot (in order: amplitude, pos amplitude, neg amplitude, latency to peak, latency to trough)
plotSelector = [1 0 0 0 0];

% whether to plot auxiliary plots (some are always plotted)
aux_plots = 0;

%% load data
homeDirectory = '../..';

resultsDirectory = [homeDirectory '/Results/12dot5Hz/'];

struct_name = 'freq12dot5hz';

fly_record = readtable([homeDirectory '/Fly record/fly_record']);

%% restrict to some frequency
fly_record = fly_record(fly_record.Frequency == 13.3,:);

%% restrictions on data of interest; always check this before running

% remove DARK flies
fly_record = fly_record(convertCharsToStrings(fly_record.Condition) == 'LIT',:);

% remove flies to be excluded (usually because data is unsound for some obvious reason)
fly_record = fly_record(~logical(fly_record.Exclude),:);

%remove jittering flies
fly_record = fly_record(~contains(fly_record.Comments,'jitter','IgnoreCase',true),:);

%remove red light flies
% fly_record = fly_record(~contains(fly_record.Comments,'red','IgnoreCase',true),:);

%restrict to red1 flies
fly_record = fly_record(contains(fly_record.Comments,'red1','IgnoreCase',true),:);

%remove block paradigm flies
fly_record = fly_record(~contains(fly_record.Comments,'block','IgnoreCase',true),:);

%% 
% fly_record = fly_record(contains(fly_record.Comments,'adjfront','IgnoreCase',true),:);

%% choose flies and experiments
whichFly =      fly_record.Fly.';
flySet = unique(whichFly);

% choose which flies to run here
% chosenFlies = [43];
% chosenFlies = setdiff(flySet, [24 25]);
chosenFlies = flySet; % choose all flies
% chosenFlies = setdiff(chosenFlies, 24:29);

%for testing
% chosenFlies = chosenFlies(1:2);

% choose which blocks to run
%NOTE: while unlikely as a request, this does not handle the case where two
%flies have a block with the same number but we would like to look at both
%flies but not one of the blocks with the same number
% chosenBlocks = [88];
chosenBlocks = unique(fly_record.Block.');% do not choose specific blocks

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
        
    if ~isnan(LFPsd(b))
    
        n_sd_lfp = LFPsd(b);
        sd_lfp = std(LFP); mean_lfp = mean(LFP);
        
        % do a nice plot
        if aux_plots
            figure; plot(LFP); hold on; %#ok<*UNRCH>
            x_lim = xlim;
            plot([x_lim(1) x_lim(2)],[mean_lfp-n_sd_lfp*sd_lfp mean_lfp-n_sd_lfp*sd_lfp],'r');
            plot([x_lim(1) x_lim(2)],[mean_lfp+n_sd_lfp*sd_lfp mean_lfp+n_sd_lfp*sd_lfp],'r');
        end
            
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

%% analyse data per fly and whether experiment is lit or unlit

for fly = chosenFlies
       
   % index of the blocks belonging to this fly if they are to be
   % included in the analysis
   thisFlyBlocks = BLOCKS(whichFly == fly & chosenOnes);

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
   
   FLIES(fly) = R; %#ok<SAGROW>
   
end

%% fit individual flies (can easily be expanded to include fits not just to amplitude)

if fitModel
    for fly = chosenFlies
        FLIES(fly).FIT = fitFlyModel(FLIES(fly));
    end
end

%% plot sequential dependencies per fly (if model results exist they are plotted)

if plotIndividualFlies
    plotFlies(FLIES, chosenFlies, plotSelector, resultsDirectory);
end

%% calculate SEs for all flies by averaging SE profiles

if length(chosenFlies) > 1
    ALLFLIES = groupFlies(FLIES, chosenFlies, errorMethod);
end

%% fit model to grouped data

if fitModel
    ALLFLIES.FIT = fitFlyModel(ALLFLIES);
end

%% plot grouped SE profiles

plotFlies(ALLFLIES,chosenFlies,plotSelector,[resultsDirectory 'All flies/']);
