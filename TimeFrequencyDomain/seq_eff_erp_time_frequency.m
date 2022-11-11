% analyse sequential dependencies in fly ERPs
close all;
clear;

%% load auxiliary functions
addpath('../');
addpath('D:\group_swinderen\Dinis\Scripts\Global functions\');
addpath('D:\group_swinderen\Dinis\Scripts\Indexes and legends\');
addpath('../Functions');

%% load data
homeDirectory = '../..';

resultsDirectory = [homeDirectory '\Results_Time_Frequency'];

fly_record = readtable([homeDirectory '\Fly record\fly_record']);

%% restrict to some frequency

frequency = '6dot25';%change this for the right file name

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
chosenBlocks = [8];
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
%for a time-frequency analysis this will always be the same
time_before_peak = fly_record.SDT;
time_after_peak = fly_record.ISI;

%peak detection threshold
peak_threshold = fly_record.Threshold;

%light_on_dark = 1 means a bright bar over a dark background was used
light_on_dark = strcmp(fly_record.Condition,'LIT').';
% chosenOnes = 1-light_on_dark; % choose all lit flies

%% whether to plot auxiliary plots
aux_plots = 0;

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

    %type of photodiode analysis
    BLOCKS(b).PHOTType = PHOTType(b);
    
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
       
            if thisFlyBlocks(1).PHOTType == 1
                R = processBlocksOneChannel(thisFlyBlocks, aux_plots,'timefrequency');
            elseif thisFlyBlocks(1).PHOTType == 2
                R = processBlocks(thisFlyBlocks, aux_plots,'timefrequency');  
            end

           FLIES(fly).(lit_dark{lit+1}) = R;
           
       end
       
   end
 
end

save(['data_flies_time_frequency_wavelet_' frequency],'FLIES');