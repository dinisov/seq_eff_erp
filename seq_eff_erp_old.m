% analyse sequential dependencies in fly ERPs
close all;
clear;

%% load auxiliary functions
addpath('../');
addpath('..\Scripts\Global functions\');
addpath('..\Scripts\Indexes and legends\');
addpath('Functions');
%addpath('C:\Users\uqmvan13\2p\2P SEs\Functions');


% load peak finder (Matt's code)
addpath('..\Scripts\Toolboxes\basefindpeaks\');

%% type of analysis

% if block experiment, choose which stimulus to look in the train (starting at 5 up to BlockLength)
focusPeak = 10;

% whether to perform time-frequency analysis (takes a long time)
timeFrequency = 1;

%% toggles

%error method (0 - a la Bruno and Matt; 1- propagation )
errorMethod = 0;
                                     
% fit model
fitModel = 1;

% plot individual flies?
plotIndividualFlies = 1;
plotComponents = 0; %Currently only applies to time/freq analysis, enables/disables SLRP, etc graphs

% which SEs to plot (in order: amplitude, pos amplitude, neg amplitude, latency to peak, latency to trough, transect, average transect in window)
plotSelector = [1 0 0 0 0 1 0];
%if plotSelector(6)
%%transectTime = 122; %frames; Always calculate
%end
%Note: This data is collated in groupFlies, collected initially in [processBlocks->analyseSequentialEffects->]calculateSEs
additionalTransectPlots = 1; %Whether to do additional transect calcs (all-timepoint transect sig, transect window, etc)

% whether to plot auxiliary plots (some are always plotted)
aux_plots = 1;

%% load data
homeDirectory = '../../Bruno';

resultsDirectory = [homeDirectory '/Results/12dot5Hz/'];

struct_name = 'freq12dot5hz';

fly_record = readtable([homeDirectory '/Fly record/fly_record']);

%% remove flies to be excluded (usually because data is unsound for some obvious reason)

fly_record = fly_record(~logical(fly_record.Exclude),:);

%% restrict to some frequency if desired

frequency = [12.5];

if ~isempty(frequency)
    fly_record = fly_record(fly_record.Frequency == frequency,:);
end

%% filter in keywords

%filter out example: {'wiChr','baseline','255'};
%filterIn = {'sri'}; %Multiple flies
filterIn = {}; %Singular
if ~isempty(filterIn)
    for i = 1:length(filterIn)
        fly_record = fly_record(contains(fly_record.Comments,filterIn{i},'IgnoreCase',true),:);
    end
end

%% filter out keywords

% cell array of keywords
filterOut = {};

fly_record = fly_record(~contains(fly_record.Comments,filterOut,'IgnoreCase',true),:);

%% choose flies and experiments
whichFly = fly_record.Fly.';
flySet = unique(whichFly);

% choose which flies to run here
chosenFlies = [24]; %Singular
%chosenFlies = setdiff(flySet, [24]); 
%chosenFlies = flySet; % choose all flies; Multiple
% chosenFlies = setdiff(chosenFlies, 24:29); 
                                                                                          
%for testing                                                                                                                                                                                                                    
% chosenFlies = chosenFlies(1:2);

% choose which blocks to run
%NOTE: while unlikely as a request, this does not handle the case where two
%flies have a block with the same number but we would like to look at both
%flies but not one of the blocks with the same number
chosenBlocks = [5];
%chosenBlocks = unique(fly_record.Block.');% do not choose specific blocks; Multiple

%chosenChannels = 1; %Only applies to multichannel data; Which channel (Singular currently) to analyse
    %Note: Not iterable across different flies/blocks (For now)
    %Currently redundant with "LFPChannel" in fly record, so unused

chosenOnes = ismember(fly_record.Block.', chosenBlocks) & ismember(fly_record.Fly.', chosenFlies);
%QA
if nansum(chosenOnes) == 0
    ['## Alert: No fly data set to be analysed ##']
    crash = yes
end
chosenFind = find(chosenOnes);
disp(['List of fly/s to be analysed:'])
for chosin = chosenFind
    disp(['Fly ',num2str(fly_record.Fly(chosin)),' Block ',num2str(fly_record.Block(chosin))])    
end
%justice

%% add data to structure according to block
% the structure may be largely empty if analysing only one fly
% the index here is that of the fly_record table, *not* the original block number

BLOCKS = collateEphysData(fly_record,chosenOnes,focusPeak,timeFrequency,homeDirectory,aux_plots);

%% analyse data for each fly (can include multiple blocks)
% IMPORTANT: make sure window is the same for all blocks belonging to the
% same fly
for fly = chosenFlies
       
   % index of the blocks belonging to this fly if they are to be
   % included in the analysis
   thisFlyBlocks = BLOCKS(whichFly == fly & chosenOnes);
   %QA
   if isempty(thisFlyBlocks)
       ['## Alert: No valid blocks found for fly #',num2str(fly),' ##']
       crash = yes
   end

   disp(['Processing fly #' num2str(fly)]);
   %disp(['(Block/s ',num2str(thisFlyBlocks),')'])

   %Quick splice in of transect params if necessary
   %if plotSelector(6)
   %thisFlyBlocks.transectTime = transectTime; %Deprecated since flyRecord encoded
   %end

   R = processBlocks(thisFlyBlocks, aux_plots, plotSelector);
   
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
    plotFlies(FLIES, chosenFlies, plotSelector, resultsDirectory); %#ok<*UNRCH>
end

%% calculate SEs for all flies by averaging SE profiles
if length(chosenFlies) > 1
    ALLFLIES = groupFlies(FLIES, chosenFlies, errorMethod);
    
    % fit model to grouped data
    if fitModel
        ALLFLIES.FIT = fitFlyModel(ALLFLIES);
    end
    
    % plot grouped SE profiles
    plotFlies(ALLFLIES,chosenFlies,plotSelector,[resultsDirectory 'All flies/']);
end

%% time-frequency analysis
if timeFrequency
    %timeFrequencyAnalysis(FLIES, '..', plotIndividualFlies);
    timeFrequencyAnalysis(FLIES, chosenFlies, '..', plotIndividualFlies, plotComponents, [8,9 ; 1,9]);
end

%% Additional transect stuff
if additionalTransectPlots
    extraTransectAnalyses(FLIES, chosenFlies, resultsDirectory,'correctForNTimepoints',1,'showPValPlots',1,...
        'patchMethod','lowestP')
end