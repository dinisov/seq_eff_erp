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
timeFrequency = 0;

%% toggles

%n-back
n_back = 5;

%error method (0 - a la Bruno and Matt; 1- propagation )
errorMethod = 0;

%Scramble level
scramLevel = 0; %What stage to scramble sequences at (0 - None, 1 - Raw sequence calculation, 2 - Derived sequence)
    %1 - Most rigorous; Scrambles shortly after sequence calculated from phot, per fly 
    %2 - Less rigorous; Scrambles the arrangement of the 16 (or etc) 
        %Use mode 1 for rigorous scrambling, mode 2 to preserve within-sequence consistency
        %Note: Mode 2 as written currently will not preserve the isomer balance (i.e. 75% of left-ending may end up right-labelled etc)
        %Secondary note: This is always applied at a per-fly level, so even with mode 2, the scrambled ordering of seqs will differ across flies

%Ordering for plots
%reOrder = [9,10,11,12,13,14,15,16,1,2,3,4,5,6,7,8]; %Reorders plots (and labels ofc) 
if n_back == 5
    reOrder = [9,10,11,12,13,14,15,16,1,2,3,4,5,6,7,8];
elseif n_back == 6
    reOrder = [17:32,1:16]; %Reorders plots (and labels ofc) 

elseif n_back == 7
    reOrder = [33:64,1:32];

elseif n_back == 4
    reOrder = [5,6,7,8,1,2,3,4];

elseif n_back == 3
    reOrder = [3,4,1,2];

else
    reOrder = [];
end
if isempty(reOrder)
    reOrder = [1:0.5*2^n_back];
end

%Whether to artificially 'shift' all ERP data above zero, to assist with certain metrics
zeroShiftMode = 1; %1 - Add the abs of the literal minimum to all elements

%Override channel
overrideChannel = []; %If non-empty will override the fly_record LFPChannel for collateEphysData and beyond
                                     
% fit model
fitModel = 0;

%Suppress ANOVA calcs
suppressANOVA = 1;

% plot individual flies?
plotIndividualFlies = 0;
plotComponents = 0; %Currently only applies to time/freq analysis, enables/disables SLRP, etc graphs

% which SEs to plot (in order: amplitude, pos amplitude, neg amplitude, latency to peak, latency to trough, transect, average transect in window)
plotSelector = [0 1 0 0 0 1 0];
%if plotSelector(6)
%%transectTime = 122; %frames; Always calculate
%end
%Note: This data is collated in groupFlies, collected initially in [processBlocks->analyseSequentialEffects->]calculateSEs
additionalTransectPlots = 0; %Whether to do additional transect calcs (all-timepoint transect sig, transect window, etc)
additionalIsomerPlots = 1; %Whether to calculate isomer correlations across time, similar to above extra transect analyses

% whether to plot auxiliary plots (some are always plotted)
aux_plots = 1;
rawDataPlot = 0;

%%

%QAs
if numel(reOrder) ~= 0.5*(2^n_back)
    ['-# Alert: Mismatch between requested order and n_back #-']
    crash = yes
end

%% load data
homeDirectory = '../../Bruno';

resultsDirectory = [homeDirectory '/Results/12dot5Hz/'];

struct_name = 'freq12dot5hz';

fly_record = readtable([homeDirectory '/Fly record/fly_record.xlsx']);

%% remove flies to be excluded (usually because data is unsound for some obvious reason)

fly_record = fly_record(~logical(fly_record.Exclude),:);

%%

selectionMode = 'keywords'; %keywords or manual; Modify this

%%

switch selectionMode

    case 'keywords' %Do not modify

        %----------------------
        %Specify keywords to select flies/blocks
        %----------------------

        % restrict to some frequency if desired
        %frequency = [12];
        frequency = [12.5];
        
        if ~isempty(frequency)
            fly_record = fly_record(fly_record.Frequency == frequency,:);
        end
        
        % filter in keywords
        filterIn = {'tsh/wichr','255','no atr','baseline'};
        %filterIn = {'Baseline'};

        if ~isempty(filterIn)
            for i = 1:length(filterIn)
                fly_record = fly_record(contains(fly_record.Comments,filterIn{i},'IgnoreCase',true),:);
            end
        end
        
        % filter out keywords
        %filter out example: {'wiChr','baseline','255'};
        % cell array of keywords
        filterOut = {};
        
        fly_record = fly_record(~contains(fly_record.Comments,filterOut,'IgnoreCase',true),:);
        
        % choose flies and experiments
        whichFly = fly_record.Fly.';
        %%flySet = unique(whichFly); %Disabled, because unique causes sorting
        
        % choose which flies to run here
        chosenFlies = whichFly; % choose all flies; Multiple
                                                                                                          
        % choose which blocks to run
        %NOTE: while unlikely as a request, this does not handle the case where two flies have a block with the same number but we would like to look at both
        %chosenBlocks = unique(fly_record.Block.');% do not choose specific blocks; Multiple; %Disabled, because do not trust unique-induced sorting
        chosenBlocks = fly_record.Block.';%Matt modification of above; No enforced uniqueness
        
        %chosenChannels = 1; %Only applies to multichannel data; Which channel (Singular currently) to analyse
            %Note: Not iterable across different flies/blocks (For now)
            %Currently redundant with "LFPChannel" in fly record, so unused
        
        chosenOnes = ismember(fly_record.Block.', chosenBlocks) & ismember(fly_record.Fly.', chosenFlies); %This implementation relies on too many assumptions (imo)
            %Unless mistaken, with uniques disabled, this line is redundant, because by definition, chosenBlocks comes from fly_record and so too for chosenFlies

    case 'manual' %Do not modify

        %----------------------
        %Specify flies/blocks manually
        %----------------------

        chosenFlies = [258]; %Singular
        %chosenBlocks = [];
        %chosenBlocks = {[26,28],[3,4,6]}; %If non-empty, must specify a block for each element of chosenFlies in the format {[<fly 1 block/s>],[<fly 2 blocks/s>], [etc]}, where multiple blocks can be selected for each fly if requested
        chosenBlocks = {[13]}; %Specify one block per fly (e.g. {[13],[17]}
            %...theoretically all aspects of this system support multiple blocks per fly (e.g. {[13,18],[1,3,5]}), but Dinis' analysis does not
                % ^ Mildly incorrect; groupBlocks (via analyseSequentialEffects) seems to support multiple blocks

        %QA
        if ~isempty(chosenBlocks) && size(chosenFlies,2) ~= size(chosenBlocks,2)
            ['## Alert: Disparity between specified chosenFlies (N=',num2str(size(chosenFlies,2)),') and chosenBlocks (N=',num2str(size(chosenBlocks,2)),') ##']
            crash = yes
        end

        %Subselect fly_record similar to other case
        chosenOnes = zeros(size(fly_record,1),1);
        chosenFliesActual = []; %This will become important in a moment
        for flyIn = 1:size(chosenFlies,2)
            thisFly = chosenFlies(flyIn);
            if ~isempty(chosenBlocks)
                for blockIn = 1:size(chosenBlocks{flyIn},2)
                    thisBlock = chosenBlocks{flyIn}(blockIn);
                    targetRow = intersect( find( fly_record.Fly == thisFly ),find( fly_record.Block == thisBlock ) ); %Find only row that matches fly/block combination
                    if isempty(targetRow) || size(targetRow,1) > 1
                        ['## Alert: Critical underfind (Non-existent fly/block?) or overfind (Ambiguous fly/block?) in fly_record for fly #',num2str(thisFly),' Block',num2str(thisBlock),' ##']
                        crash = yes
                    end  
                    chosenOnes(targetRow) = 1; %Assign this row to be collected
                    chosenFliesActual = [chosenFliesActual,thisFly];
                end
            else
                chosenOnes( find( fly_record.Fly == thisFly ) ) = 1; %Assign all rows associated with this fly (regardless of block) to be collected
                chosenFliesActual = [chosenFliesActual,repmat(thisFly,1,[nansum(fly_record.Fly == thisFly)])];
            end
        end
        %Report
        disp([num2str(nansum(chosenOnes)),' fly/block combinations will be analysed'])
        fly_record( ~chosenOnes,: ) = []; %Trim to only applicable fly/blocks
        chosenOnes( ~chosenOnes ) = []; %Trim self
        chosenOnes = chosenOnes'; %Because reasons

        %Adjust structure of chosenFlies to match how it is for keywords (namely, that multiple blocks for one fly result in multiplication of fly number)
        chosenFlies = chosenFliesActual;
        whichFly = chosenFlies; %Necessary for synchrony with keywords system


end

%%

%QA
if nansum(chosenOnes) == 0
    ['## Alert: No fly data set to be analysed ##']
    crash = yes
end
if numel(unique(chosenFlies)) ~= numel(chosenFlies)
    ['## Alert: Multiple blocks requested for one or more flies, which is not supported ##']
    crash = yes
end
chosenFind = find(chosenOnes);
disp(['List of fly/s to be analysed:'])
for chosin = chosenFind
    disp(['Fly ',num2str(fly_record.Fly(chosin)),' Block ',num2str(fly_record.Block(chosin)),' (',datestr(fly_record.Date(chosin)),')'])    
end
%justice

%% add data to structure according to block
% the structure may be largely empty if analysing only one fly
% the index here is that of the fly_record table, *not* the original block number

%BLOCKS = collateEphysData(fly_record,chosenOnes,focusPeak,timeFrequency,homeDirectory,aux_plots);
BLOCKS = collateEphysData(fly_record,chosenOnes,focusPeak,timeFrequency,homeDirectory,aux_plots, zeroShiftMode,'+',overrideChannel,rawDataPlot); %'-' for multichannel, '+' for single

%% analyse data for each fly (can include multiple blocks)
% IMPORTANT: make sure window is the same for all blocks belonging to the
% same fly
for fly = chosenFlies
       
   % index of the blocks belonging to this fly if they are to be
   % included in the analysis
   thisFlyBlocks = BLOCKS(whichFly == fly & chosenOnes); %Note: chosenOnes here mostly redundant, since pretty much just fullsize series of positive logicals

   %QA
   if isempty(thisFlyBlocks)
       ['## Alert: No valid blocks found for fly #',num2str(fly),' ##']
       crash = yes
   end

   disp(['Processing fly #' num2str(fly)]);
   %disp(['(Block/s ',num2str(thisFlyBlocks),')'])
   for bInd = 1:size(thisFlyBlocks,2)
    disp(['Block/s: ',thisFlyBlocks(bInd).block])
   end

   %Quick splice in of transect params if necessary
   %if plotSelector(6)
   %thisFlyBlocks.transectTime = transectTime; %Deprecated since flyRecord encoded
   %end

   R = processBlocks(thisFlyBlocks, aux_plots, plotSelector, reOrder, n_back, ...
       'suppressANOVA',suppressANOVA, 'plotIndividualFlies',plotIndividualFlies,...
       'scramLevel',scramLevel);
   
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
    plotFlies(FLIES, chosenFlies, plotSelector, resultsDirectory, reOrder, n_back, scramLevel); %#ok<*UNRCH>
end

%% calculate SEs for all flies by averaging SE profiles
if length(chosenFlies) > 1
    %ALLFLIES = groupFlies(FLIES, chosenFlies, errorMethod);
    ALLFLIES = groupFlies(FLIES, chosenFlies, errorMethod, reOrder, n_back,'suppressANOVA',suppressANOVA);
    
    % fit model to grouped data
    if fitModel
        ALLFLIES.FIT = fitFlyModel(ALLFLIES);
    end
    
    % plot grouped SE profiles
    plotFlies(ALLFLIES,chosenFlies,plotSelector,[resultsDirectory 'All flies/'], reOrder, n_back, scramLevel);
end

%% time-frequency analysis
if timeFrequency
    %timeFrequencyAnalysis(FLIES, '..', plotIndividualFlies);
    timeFrequencyAnalysis(FLIES, chosenFlies, '..', plotIndividualFlies, plotComponents, [9,1 ; 8,1 ; 9,16 ; 8,16],n_back);
end

%% Additional transect stuff
if additionalTransectPlots
    extraTransectAnalyses(FLIES, chosenFlies, resultsDirectory,'correctForNTimepoints',1,'showPValPlots',1,...
        'patchMethod','lowestP', 'n_back',n_back,'reOrder',reOrder, 'plotIndividualFlies',plotIndividualFlies)
end

%% Additional isomer stuff
if additionalIsomerPlots
    [CROSSISOMER] = extraIsomerAnalyses(FLIES, chosenFlies, resultsDirectory,'correctForNTimepoints',1,'showPValPlots',1,...
        'patchMethod','lowestP', 'n_back',n_back,'reOrder',reOrder,'plotSelector',plotSelector,...
        'plotIndividualFlies',plotIndividualFlies, ...
        'timeStep',20, 'doAnimatedPlot',0, 'saveFigVid',0, 'extraFigSubplots',2);
end

%% Report about potential scrambling

if scramLevel ~= 0
    ['-# Caution: Scrambling was enabled for this analysis #-']
end

