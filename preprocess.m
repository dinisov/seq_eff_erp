%The script formerly known as preprocess_01_converttdt_Calibration_AllFlies_NoFilter1000Hz
%Mk 1 - Core functionality
%Mk 2 - Folder generalisation (for multiple experiments on one day)
%Mk 3 - Trimmed down use of cd
%    .25 - Network mode, Fixed chunk/block naming issues
%    .5 - Synapse support, MATLAB 2014b compatibility
%    .75 - Empty block armouring improvements
%    .85 - Synapse/Bruno improvements (12/5/21)
%Mk 4 - Integration of calib functionality
%    .25 - Generalisation of data import from Synapse (But not non-Synapse data)
%    .5 - Added subsampling alternative for resampling (04/06/21)

%To do: Add subsampling for inputData

%%clc 
clear
close all

%% Load all paths here.. 
%-------------------------------------------------
%Params (Mostly)
networkMode = 0; %Whether to load data from network rather than local
if networkMode == 0
    %dataPath = 'D:\group_swinderen\Bruno\Data\'
    %outPath = 'D:\group_swinderen\Bruno\Processed\'
    dataPath = 'D:\group_swinderen\Dinis\Data\';
    outPath = 'D:\group_swinderen\Dinis\Output\';
    %dataPath = 'D:\group_swinderen\Matthew\TDTs\Data\'
    %outPath = 'D:\group_swinderen\Matthew\TDTs\Processed\'
else
    dataPath = 'I:\PHDMVDP002-Q1471\LFP\Analysis\Data\';
    outPath = 'I:\PHDMVDP002-Q1471\LFP\Analysis\Processed\';
end
%-----------------------------------------------
%Flags (Mostly)
expName = 'TagTrials'; %This is the name that will be looked for in all detected folders that exist within dataPath (probably)
    %Note: With Mk 4.25 wildcards no longer necessary here (And will in fact cause issues for tank detection)
calibMode = 0; %Whether to analyse calib data
if calibMode == 1
    findPolRevers = 1; %Whether to try find the polarity reversal for calb data
end
dataIsFromSynapse = 1; %Whether data was collected from Synapse
if dataIsFromSynapse == 1
    specialSubsampleFields = []; %These fields will be subsampled with a coordinate technique, rather than resampled
        %This is intended for TTL and TTL-like fields that react...poorly to traditional resampling methods
end
blockConvention = 'block'; %This is the string that will be searched for in folders (e.g. if your convention is 'C:\TDT\Tanks\290421\block1' etc, then the blockConvention is 'block' and the value immediately after is the block number)
reSampleFreq = NaN; % Desired Sampling Frequency (Set as NaN to use source framerate)
skipExistingFiles = 1; %Whether to skip analysis of already processed files
%-------------------------------------------------

%cd(dataPath)
%paths
%addpath(genpath(['D:\group_swinderen\Matthew\Scripts\toolboxes\TDTMatlabSDK'])) %Required
%addpath(genpath(['D:\group_swinderen\Matthew\Scripts\toolboxes\basefindpeaks']))
addpath(genpath(['D:\group_swinderen\Dinis\Scripts\Toolboxes\TDTMatlabSDK'])) %Required
addpath(genpath(['D:\group_swinderen\Dinis\Scripts\Toolboxes\basefindpeaks']))
%cd(raw_path);

%Begin actual processing
%fly_list = dir('*290520*'); % changed from Analyzed*
fly_list = dir([dataPath, filesep, '*']); % changed from Analyzed*
%expName = '*290421*';

fly_list = fly_list(convertCharsToStrings({fly_list.name}) == '200122',:);

if length(fly_list) == 0
    ['## ERROR: NO DATA FOUND ##']
    crash = yes
end
%Remove operating system 'folders' from fly_list
osExcludeList = [{'.'},{'..'}];
for osInd = 1:size(osExcludeList,2)
    for flyInd = size(fly_list,1):-1:1
        %if ( contains( fly_list(flyInd).name , osExcludeList{osInd} ) == 1 && size( fly_list(flyInd).name,2 ) == size(osExcludeList{osInd},2) )
        if ( isempty( strfind( fly_list(flyInd).name , osExcludeList{osInd} ) ) ~= 1 && size( fly_list(flyInd).name,2 ) == size(osExcludeList{osInd},2) )
            fly_list(flyInd) = [];
        end
    end
end
%Remove non-directory hits
for flyInd = size(fly_list,1):-1:1
    if fly_list(flyInd).isdir ~= 1
        fly_list(flyInd) = [];
    end
end

for fly_number = 1:length(fly_list)
    disp([char(10),'#######################################################################################################'])
    % locpathappend = Path1;

    %% Dataset path here..
    datasetname = [fly_list(fly_number).name];

    datasetname

    %datasetPath = [dataPath,filesep, datasetname];
    datasetPath = [dataPath, datasetname]; %dataPath contains a filesep under normal conditions anyway
        %Folder structure reminder:
        %"...\Data\<Date of experiment (i.e. 031220)>\<Block name (i.e. 031220_Overnight)>\Block <1 : etc>"

    %{
    %logfile_path = ['G:\Rhiannon\Ephys Data\Sleep\LFP\' datasetname filesep];

    %  if ~isdir(logfile_path)
    %             mkdir(logfile_path);
    %  end

    % cd(logfile_path);
    % diary('tdt_calibration_diary.txt');
    % fprintf('---The log started at %s ----',datestr(now));%Displaying this for the log file..
    % 

    % try
    % data_path = [data_root filesep datasetname filesep 'lfp' filesep];
     %cd(data_path)
     %catch
     %    disp('Trying other directory structure...')
     %   data_path = [data_root filesep datasetname filesep];
     %cd(data_path) 
     %end

    %%cd(datasetname)

    %all_dir = dir('*oddball');
    %all_dir = dir(expName);
    
    %{
    try
    allsetname = all_dir(1).name;
    catch
        warning('No all folder detected. Skipping file.')
        continue
    end
    %}
    %}
    %all_dir = dir([datasetPath, filesep, expName]);
    all_dir = dir([datasetPath, filesep, '*']);
    validAllSets = [];
    a = 1;
    for i = 1:size(all_dir,1)
        if calibMode == 0 %"Not running in calib preprocess mode"
            if isempty(strfind(all_dir(i).name, 'Calib')) == 1 && isempty(strfind(all_dir(i).name, '.')) == 1 && isempty(strfind(all_dir(i).name, expName)) ~= 1
                validAllSets{a} = all_dir(i).name;
                a = a + 1;
            end
        else %"Yes calib preprocess"
            if isempty(strfind(all_dir(i).name, 'Calib')) ~= 1 && isempty(strfind(all_dir(i).name, '.')) == 1
                validAllSets{a} = all_dir(i).name;
                a = a + 1;
            end
        end
    end
    disp(['-- ', num2str(a-1), ' valid datasets detected for ',datasetname, ' --'])
    
    if size(validAllSets,2) == 0
        ['## ALERT: NO VALID "',expName,'" DATA DETECTED FOR ', datasetname, ' ##']
        continue
    end

    %Synapse compatibility
    %{
    if synapseMode == 1
        subExpTimes = struct;
        if nansum( contains( [ validAllSets ] , 'block' ) ) > 0.5 * size(validAllSets,2) %"At least half of all folders have friendly block naming convention"
            disp(['-- Calculating block numbers from folder names... --'])
            for subExpNum = 1:size(validAllSets,2)
                unPos = strfind( validAllSets{subExpNum} , 'block' );
                %QA
                if isempty(unPos) == 1
                    ['## Invalid block name! ##']
                    crash = yes
                end
                blockNum = str2num( validAllSets{subExpNum}(unPos + size('block',2) : end) );
                subExpTimes = crash %Code not finished
            end
        else
            disp(['-- Calculating block numbers from Synapse naming convention... --'])
            %%subExpTimes = struct;
            for subExpNum = 1:size(validAllSets,2)
                thisSubExpName = validAllSets{subExpNum};
                potentialBlockList = dir([datasetPath, filesep, thisSubExpName]);
                for potInd = size(potentialBlockList,1):-1:1
                    if potentialBlockList(potInd).isdir == 1
                        filesInPotentialBlock = dir( [datasetPath, filesep, thisSubExpName, filesep, potentialBlockList(potInd).name, filesep, '*.tev'] ); %This will be attempted for the OS special hits ("." and ".."), but shouldn't match
                            %(Unless there are .tev files in the root block directory for some reason...)
                        if isempty(filesInPotentialBlock) == 1
                            potentialBlockList(potInd) = [];
                        end
                    else
                        potentialBlockList(potInd) = [];
                    end 
                end
                for blockInd = 1:size(potentialBlockList,1)
                    thisBlockName = potentialBlockList(blockInd).name;
                    hypLocs = strfind( thisBlockName , '-'); %Find locations of hyphens
                    if size(hypLocs,2) < 2
                        ['## Alert: Invalid naming format detected for ',validAllSets{subExpNum}, ' ##']
                        crash = yes %Will happen if filename is not something like "290421_f1-210429-121244"
                    end
                    dataTimeStamp = thisBlockName( hypLocs(end-1) + 1 : end ); %Use second-to-last hyphen as start of timestamp, collect from there to end
                    dataTime = datetime(strcat(dataTimeStamp),'Format', 'yyMMdd-HHmmss', 'TimeZone', '+10:00');

                    %thisBlockFieldName = ['Exp_',thisBlockName];
                    %thisBlockFieldName = strrep( thisBlockFieldName , '-' , '_' );
                    thisSubFieldName = ['SubExp_',thisSubExpName];
                    thisSubFieldName = strrep( thisSubFieldName , '-' , '_' );
                    eval(['subExpTimes.',thisSubFieldName,'(blockInd).subExpNum = subExpNum;']);
                    eval(['subExpTimes.',thisSubFieldName,'(blockInd).dirInd = blockInd;']);
                    %eval(['subExpTimes.',thisSubFieldName,'(blockInd).name = thisBlockFieldName;']);
                    eval(['subExpTimes.',thisSubFieldName,'(blockInd).name = thisBlockName;']); %Unmodified, because unnecessary
                    eval(['subExpTimes.',thisSubFieldName,'(blockInd).time = dataTime;']);
                    eval(['subExpTimes.',thisSubFieldName,'(blockInd).posixtime = posixtime(dataTime);']);
                    %subExpTimes.(thisSubFieldName)(blockInd).name = thisBlockFieldName;
                    %subExpTimes.(thisSubFieldName)(blockInd).time = dataTime;
                    %subExpTimes.(thisSubFieldName)(blockInd).posixtime = posixtime(dataTime);
                end

                %Sort by time (Honestly probably redundant but good protocol)
                temp = struct2table(subExpTimes.(thisSubFieldName)); % convert the struct array to a table
                temp = sortrows(temp, 'posixtime'); % sort the table by blocknumber
                    %Note: Unstable use of temp but w/e
                subExpTimes.(thisSubFieldName) = table2struct(temp);

                %Now append correct block numbers
                for subEx = 1:size(subExpTimes.(thisSubFieldName),1)
                    subExpTimes.(thisSubFieldName)(subEx).blockNum = subEx;
                end

            end
        end
        
    end
    %}

    varSaveList = who;
    varSaveList = [varSaveList; {'varSaveList'}; {'subExpNum'}];
    
    for subExpNum = 1:size(validAllSets,2)        
        close all
        clearvars('-except',varSaveList{:})
        thisSubExpName = validAllSets{subExpNum};
        disp(['------------------------------------------------------------------------'])
        disp(['-- Now processing ',thisSubExpName, ' (',num2str(subExpNum), ' of ', num2str(size(validAllSets,2)), ') --'])
        
        allsetname = validAllSets{subExpNum};
        %cd(allsetname)

        %%cd([dataPath, filesep, datasetname, filesep, allsetname]) %Full location referencing because relative referencing fails in loops
        % calib_data_path = [data_path filesep calibsetname];
        % cd(calib_data_path)
        %allsetPath = [dataPath, filesep, datasetname, filesep, allsetname];
        allsetPath = [dataPath, datasetname, filesep, allsetname]; %Removed filesep because terminal filesep already existing in dataPath

        detExpName = strrep(allsetname, datasetname, ''); %Determine experiment title
        while nanmin(strfind(detExpName,'_')) == 1
            detExpName(1) = ''; %Remove proximal underscores
        end

        %block_list = dir('Block*'); % changed from Analyzed*
        %Old
        %{
        if synapseMode ~= 1 %"Not Synapse mode"
            blockMode = 0; %Determines whether blocks are readily detectable (OpenEx) or if they are in a more...arcane format (Synapse)
            block_list = dir([allsetPath,filesep,'*',blockConvention,'*']); % changed from Analyzed*
            if isempty(block_list) == 1 && synapseMode == 1
                disp(['#- Warning: No detected blocks with name "Block"; Attempting with expName ("',expName,'") instead #-'])
                block_list = dir([allsetPath,filesep,expName]);
                blockMode = 1; %Switch to Synapse mode
            elseif isempty(block_list) == 1 && synapseMode == 0
                ['## Alert: No valid block detected ##']
                crash = yes
            end
        else
            thisSubFieldName = ['SubExp_',thisSubExpName];
            thisSubFieldName = strrep( thisSubFieldName , '-' , '_' );
            
            %Assemble block_list
            %block_list = dir([allsetPath,filesep,expName]);
            block_list = subExpTimes.(thisSubFieldName);
            
            %{
            %Find position in subExpTimes, to derive block number
            match = [];
            for subEx = 1:size(subExpTimes.(thisSubFieldName),1)
                if subExpTimes.(thisSubFieldName)(subEx).dirInd == subExpNum
                    match = subEx;
                    continue
                end
            end
            thisBlockNum = subExpTimes(match).blockNum;
            %}            
        end
        %}
        %New
        blockMode = 0; %Determines whether blocks are readily detectable (OpenEx) or if they are in a more...arcane format (Synapse)
        block_list = dir([allsetPath,filesep,'*',blockConvention,'*']); % changed from Analyzed*
        if isempty(block_list) && dataIsFromSynapse == 1
            disp(['#- Warning: No detected blocks with name "Block"; Attempting with expName ("',expName,'") instead #-'])
            block_list = dir([allsetPath,filesep,'*',expName,'*']);
            blockMode = 1; %Switch to Synapse mode
        elseif isempty(block_list) == 1 && dataIsFromSynapse == 0
            ['## Alert: No valid block detected ##']
            crash = yes
        end

        numErrorChunks = 0; %Hopefully will not iterate, but may do if there are empty data blocks/etc
        thisChunkIsError = 0; %Will be overwritten later if case
        
        totalchunks = length(block_list);
        %Old
        %{
        if synapseMode == 0
            %Clean and sort block_list
            for blockInd = 1:size(block_list,1)
                %unPos = strfind(block_list(blockInd).name,'-'); %Old, searched for hyphen
                unPos = strfind(block_list(blockInd).name,blockConvention); %New, searches for word 'block'
                if isempty(unPos) ~= 1
                    block_list(blockInd).blocknumber = str2num( block_list(blockInd).name(unPos+size(blockConvention,2):end) );
                else
                    ['## Alert: Block identity could not be derived from block name ##']
                    crash = yes
                end
            end
        else
            %Figure out block numbers from subExpTimes
            for blockInd = 1:size(block_list,1)
                block_list(blockInd).blocknumber = block_list(blockInd).blockNum; %Mirror blockNum across all block files
            end
        end
        %}
        %New
        for blockInd = 1:size(block_list,1)
            %unPos = strfind(block_list(blockInd).name,'-'); %Old, searched for hyphen
            unPos = strfind(block_list(blockInd).name,blockConvention); %New, searches for word 'block'
            if isempty(unPos) ~= 1
                block_list(blockInd).blocknumber = str2num( block_list(blockInd).name(unPos+size(blockConvention,2):end) );
            else
                ['## Alert: Block identity could not be derived from block name ##']
                crash = yes
            end
        end
        temp = struct2table(block_list); % convert the struct array to a table
        temp = sortrows(temp, 'blocknumber'); % sort the table by blocknumber
            %Note: Unstable use of temp but w/e
        block_list = table2struct(temp);

        for chunkidx = 1:totalchunks
            
            processingStartTime = clock;
            %%Step 1: Extract data in segments of 1 hour each..
            %blockname = ['Block-' char(sprintf("%d",chunkidx))]; %Old
            blockname = [block_list(chunkidx).name]; %New
            %blockpath = [allsetname filesep blockname];
            %%blockPathFull = [dataPath, filesep, datasetname, filesep, blockpath];
            blockPathFull = [allsetPath, filesep, blockname];
            blockID = blockname(strfind(blockname, '-')+1:end);
            if size(blockID,2) < 2
                blockID = ['0',blockID]; %Cover for up to block 99
            end

            disp([char(10),'-- Reading detected block #',num2str(chunkidx),' (',blockname,') of ',num2str(totalchunks),' --'])
            %blockpath = char(horzcat(blockpath(1:length(blockpath))));

            locpathappend = [outPath, datasetname];
            %outputfolder = [locpathappend filesep 'LFP' filesep 'Analyzed_LFP_' blockname];
            outputfolder = [locpathappend filesep 'LFP' filesep 'Analyzed_' detExpName '_' blockname];
            S.output_folder  = outputfolder;
            %S.eeg_filename = [datasetname  '_chunk_' char(sprintf("%0.2d",chunkidx))];
            S.eeg_filename = [datasetname,  '_chunk_', blockID]; %Dynamic detection of block/chunk number
            mat_name = [S.output_folder filesep S.eeg_filename '.mat'];
            
            try
                fileIsExist = isfile(mat_name); %Will fail on MATLAB 2014b, but that's why this is in a try-catch
            catch
                fileIsExist = exist(mat_name);
            end
            if fileIsExist ~= 0 && skipExistingFiles == 1
                warning('Block already processed. Skipping file.');
                continue
            end
            %check the timeduration to start and end..
            tdt_dur_load = 0;
            a = 0; %Iterator to prevent infiniloop
            %tdt_load_count =0;
            tic
            while tdt_dur_load == 0 && a < 10
                try
                 %[timeduration, info] = TDTduration(blockpath);
                    [timeduration, info] = TDTduration(blockPathFull);
                    disp('TDTduration read successful.');
                    tdt_dur_load = 1;
                catch
                    disp(['TDTduration read failed; Retrying (in 10s)...'])
                    pause(10)
                    %{
                    tdt_load_count =tdt_load_count+1;
                    pause(30);
                    warning(['TDTduration read failed. Trying again after 30 seconds. Attempt number ' num2str(tdt_load_count)])
                    if tdt_load_count == 30;          
                    warning(['Error processing file.' blockpath 'Skipping to the next one.']);
                    tdt_dur_load = 1;
                    continue
                    end
                    %}
                    a = a + 1;
                end
            end
            if tdt_dur_load == 0
                ['## Alert: TDT duration could not be successfully read ##']
                crash = yes
            end
            
            disp(['TDTduration performed in ',num2str(toc),'s'])

            %fprintf('Reading detected block %d of %d..\n', chunkidx, totalchunks);
            

            %fprintf('Total time duration of recording is: %0.2f secs..\n', timeduration);

            if ~isdir(locpathappend)
                mkdir(locpathappend);
            end

            %check the timeduration to start and end.

            %fprintf('Reading block %d of %d..\n', chunkidx, totalchunks);
 
            fprintf('Total time duration of recording is: %0.2f secs..\n', timeduration);
            fprintf('Equivalent to: %0.2f minutes..\n', timeduration/60);

            starttime = 0;
            endtime = timeduration;
            
            %#########################
            %QA for corrupted chunk
            if isfield(info,'headerstoptime') ~= 1
                ['-# Warning: headerstoptime is missing from chunk; Approximating value #-']
                info.headerstoptime = info.headerstarttime + endtime; %May not be identical to proper data
            end
            %#########################
                        
            %data = TDTbin2mat(blockpath, 'T1', starttime, 'T2', endtime);
            block_read = 0;
            %block_read_count = 0;
            a = 0;
            while block_read == 0 && a < 10
                try
                    data = TDTbin2mat(blockPathFull);
                    disp('TDTbin2mat read successful.')
                    block_read = 1;
                catch
                    disp(['TDTbin2mat read failed; Retrying (in 10s)...'])
                    pause(10)
                    a = a + 1;
                    %TDTbin2mat
                    %{
                    warning('TDTbin2mat read unsuccessful. Trying again in 30 seconds...')
                    pause(30)
                    block_read = 0;
                    block_read_count = block_read_count +1;
                    if block_read_count == 5
                        warning('TDTbin2mat read failed too many times. Skipping file.');
                        block_read = 1;
                        continue
                    end
                    %}
                end
            end
            if block_read == 0
                ['## Alert: Detected chunk #',num2str(chunkidx),' - ',blockname,' could not be read ##']
                crash = yes
            end

            %Returned data
            try
                if dataIsFromSynapse == 0
                    waveData = data.streams.Wave;
                else
                    %{
                    try
                        waveData = data.streams.Wav1;
                    catch
                        waveData = data.streams.LFP1;
                    end
                    %}
                    dataFiels = fieldnames(data.streams);
                    for fielInd = 1:size(dataFiels,1)
                        waveData.(dataFiels{fielInd}) = data.streams.(dataFiels{fielInd}); %Note: This means that Synapse data has a fundamentally different (Read: Expanded) architecture compared to non-Synapse
                    end
                end
            catch %Most likely in case of real data not existing
                ['## Warning: Error collecting real data for chunk #',num2str(chunkidx),' - ',blockname,' ##']
                thisChunkIsError = 1; %Flag
                waveData = struct;
                waveData.data = []; %Blank these for less effort error handling
                waveData.fs = reSampleFreq;  %This cannot be blank else the rat fails
                chandata_resamp = []; %Generate a blank value for this ahead of time, since it won't be generated by the normal loop (Because chandata size 0)
            end
            
            if dataIsFromSynapse == 0
                chandata = waveData.data; %Old, assumed Single format
            else
                chandata = waveData; %New
            end
            
            %Input data
            try
                if dataIsFromSynapse == 0
                    %stimdata = data.streams.InpP.data(:,:);
                    inpData = data.streams.InpP;
                else
                    if isfield(data.streams,'Phot') == 1
                        inpData = data.streams.Phot; %No generalisation here currently
                    else
                        inpData = [];
                    end
                end
            catch
                ['## Warning: Error collecting stim data for chunk #',num2str(chunkidx),' - ',blockname,' ##']
                thisChunkIsError = 1; %Flag
                inpData = struct;
                inpData.data = []; %Blank these for less effort error handling
                inpData.fs = reSampleFreq; %This cannot be blank else the rat fails
                stimdata_resamp = []; %Generate a blank value for this ahead of time, since it won't be generated by the normal loop (Because chandata size 0)
            end
            try
                stimdata = inpData.data(:,:);
            catch
               stimdata = []; 
            end

            %     % Keep the channel input streams the same as is.
            %     switch 0
            %         case 1
            %     %now just to do some remapping of stimulus resampling..
            %     remapidx = find(abs(stimdata(:,:))<1400000);
            %     stimdata(1,remapidx) = 0;
            %     remapidx = find(stimdata(:,:)<-1400000);
            %     stimdata(1,remapidx) = -500;
            %     remapidx = find(stimdata(:,:)>1400000);
            %     stimdata(1,remapidx) = 500;
            %     end

            %set the starttime in the eeglab struct..
            eegtimestart = info.headerstarttime + starttime;
            eegtimeend = info.headerstarttime + endtime;

            %%Step 2:  Downsample the data..
            
            %Switch between modes depending on whether chandata contains all data (Synapse) or just one channel (non-Synapse)
            if ~isstruct(chandata)
                %Just one channel
                %inputfreq = data.streams.Wave.fs;% Actual Sampling Frequency
                inputfreq = waveData.fs;% Actual Sampling Frequency (Now dynamic between Synapse and OpenEx)
                %resamplefreq = reSampleFreq; 
                if isnan(reSampleFreq) ~= 1
                    resamplefreq = reSampleFreq; % Desired Sampling Frequency (Why is this a Re;Zero reference? Because I can that's why)
                else
                    resamplefreq = inputfreq;
                end
                [N,D] = rat(resamplefreq/inputfreq); % Rational Fraction Approximation

                for idx = 1:size(chandata,1)
                    chandata_resamp(idx,:) = resample(double(chandata(idx,:)), N, D);% Resampled Signal   
                end
            
            else
                %Iterate across all
                chandata_resamp = []; 
                for fielInd = 1:size(dataFiels,1) %Use dataFiels from above
                    
                    thisFiel = dataFiels{fielInd};
                    inputfreq = waveData.(thisFiel).fs;% Actual Sampling Frequency (Now dynamic between Synapse and OpenEx)
                    thisData = chandata.(thisFiel).data;
                    if ~isnan(reSampleFreq)
                        resamplefreq = reSampleFreq; % Desired Sampling Frequency (Why is this a Re;Zero reference? Because I can that's why)
                    else
                        resamplefreq = inputfreq;
                    end
                    [N,D] = rat(resamplefreq/inputfreq); % Rational Fraction Approximation
                    %Aborted alternative method system
                    %totalDuration = size(thisData,2) / inputfreq;
                    %timepoints = linspace(0,totalDuration,size(thisData,2));
                    %blorg = resample( double(thisData) , timepoints, 200 );
                    subSampleProceed = 0; %Will be adjusted if such
                    for subInd = 1:size(specialSubsampleFields,2)
                        if nansum( thisFiel == specialSubsampleFields{subInd} ) == size(thisFiel,2) %"Field name is perfect match for an element of specialSubsampleFields"
                            subSampleProceed = 1;
                            disp(['-- Field ',thisFiel,' will be subsampled as per request --'])
                        end
                    end
                    
                    if subSampleProceed == 1
                        %'New' artefact-free subsampling
                        subCoords = round([1 : inputfreq/reSampleFreq : size(thisData,2)]); 
                            %Note: Will introduce variably sized timing errors depending on much of a not-multiple reSampleFreq is into inputfreq
                        %Quick supersampling QA
                        if reSampleFreq > inputfreq
                            ['-# Caution: reSampleFreq (',num2str(reSampleFreq),'Hz) exceeds inputfreq (',num2str(inputfreq),'Hz); Oversampling has occurred #-']                            
                        end
                        chandata_resamp.(thisFiel) = []; %Clear each time
                        for idx = 1:size(thisData,1)
                            chandata_resamp.(thisFiel)(idx,:) = nan( 1, size( resample(double(thisData(idx,:)), N, D) ,2 ) ); %Use alternative system to derive proper length
                            chandata_resamp.(thisFiel)(idx,1:size(subCoords,2)) = double( thisData(idx,subCoords) );% Subsampled Signal
                                %Note: This method can be up to inputfreq/reSampleFreq points different in length from data sampled the other way
                                    %E.g. Subsampling a 3Khz signal at 200Hz means that every ~15th point is grabbed and thus the resampled data
                                    %can be only ever exactly as long or up to 15 points shorter (based on how MATLAB iterators work)
                        end    
                    else
                        %'Old' system
                        %[N,D] = rat(resamplefreq/inputfreq); % Rational Fraction Approximation
                        chandata_resamp.(thisFiel) = []; %Clear each time
                        for idx = 1:size(thisData,1)
                            %chandata_resamp.(thisFiel)(idx,:) = resample(double(chandata.(thisFiel).data(idx,:)), N, D);% Resampled Signal  
                            chandata_resamp.(thisFiel)(idx,:) = resample(double(thisData(idx,:)), N, D);% Resampled Signal  
                        end
                    end
                    
                end
                
            end

            try
                %inputfreq = data.streams.InpP.fs;% Actual Sampling Frequency
                inputfreq = inpData.fs;% Actual Sampling Frequency
                [N,D] = rat(resamplefreq/inputfreq); % Rational Fraction Approximation

                for idx = 1: size(stimdata,1)
                    stimdata_resamp(idx,:) = resample(double(stimdata(idx,:)), N, D);% Resampled Signal   
                        %Currently no subsampling alternative for inputData
                end
            catch
                disp(['#- InpP data not existing for Chunk #',num2str(chunkidx),' #-'])
                stimdata_resamp = [];
            end
            
            %Calculate size disparity between virtual and real channels, if any
            if isempty(stimdata) ~= 1
                virtualDur =  ( size(inpData.data,2) / inpData.fs );
                if dataIsFromSynapse == 0
                    dataDur = ( size(waveData.data,2) / waveData.fs );
                        %Note: Any slight Fs inaccuracies will result in a potentially false positive here
                    if abs(dataDur - virtualDur) >= 0.001 %Disparity larger than 1ms
                        ['# Warning: Size disparity of approx. ',num2str(dataDur - virtualDur),'s exists between real and virtual channels #']
                    end
                else
                    for fielInd = 1:size(dataFiels,1)
                        dataDur = ( size(waveData.(dataFiels{fielInd}).data,2) / waveData.(dataFiels{fielInd}).fs );
                        %Note: Any slight Fs inaccuracies will result in a potentially false positive here
                        if abs(dataDur - virtualDur) >= 0.001 %Disparity larger than 1ms
                            disp(['# Warning: Size disparity of approx. ',num2str(dataDur - virtualDur),'s exists between ',dataFiels{fielInd},' and virtual channel/s #'])
                        end
                    end
                end
            end
            
            %length(chandata)
            %%clear chandata stimdata %data %data not cleared so as to be used slightly lower
            %LFP = [chandata_resamp; stimdata_resamp];
            LFP = chandata_resamp;
            STIM = stimdata_resamp;
            %%clear chandata_resamp stimdata_resamp


            %%Step 3:  Store it in a EEGlab file.. % note RJ: this is not clock time, but time of recording
            samplerate = resamplefreq; 
            timeres = 1/samplerate;
            Tottime = endtime - starttime;
%             timepoints = 0: timeres : Tottime;
            timepoints = linspace(0,Tottime,length(chandata_resamp.LFP1));%DINIS: used LFP1 assuming the other channels are the same.

            %%Step 4: Create EEGlab file for LfP data..
            EEG = [];
            EEG.setname = [blockname];
            EEG.filename = [datasetname];
            % EEG.filepath = [calib_data_path];
            if ~isstruct(LFP) %Probably non-Synapse
                EEG.nbchan = size(LFP,1);
                EEG.data = double(LFP);
                EEG.pnts = length(LFP);
                EEG.chanlocs(1).labels = 'LfP';
            else
                for fielInd = 1:size(dataFiels,1)
                    thisFiel = dataFiels{fielInd};
                    EEG.(thisFiel).nbchan = size(LFP.(thisFiel),1);
                    EEG.(thisFiel).data = double(LFP.(thisFiel));
                    EEG.(thisFiel).pnts = length(LFP.(thisFiel));
                    EEG.(thisFiel).chanlocs(1).labels = thisFiel;
                end
            end
            
            EEG.stims = double(STIM);
            EEG.times  = timepoints;
            EEG.xmax = max(EEG.times);
            EEG.xmin = min(EEG.times);
            

            EEG.icawinv =[];
            EEG.icaweights =[];
            EEG.icasphere =[];
            EEG.icaact = [];
            EEG.trials = 1;
            EEG.srate = samplerate;
            
            %Matt additions of useful
            try
                %EEG.sourceFramerates.InpP = data.streams.InpP.fs;
                EEG.sourceFramerates.InpP = inpData.fs;
            catch
                EEG.sourceFramerates.InpP = [];
            end
            %EEG.sourceFramerates.Wave = data.streams.Wave.fs;
            if isfield(waveData,'fs') == 1
                EEG.sourceFramerates.Wave = waveData.fs;
            else
                for fielInd = 1:size(dataFiels,1)
                    thisFiel = dataFiels{fielInd};
                    EEG.(thisFiel).sourceFramerates.Wave = waveData.(thisFiel).fs;
                end
            end

            %evalexp = 'eeg_checkset(EEG)'; %Seems to be redundant
            % [T,EEG] = evalc(evalexp);

            tmpval = datestr(datenum([1970, 1, 1, 0, 0, eegtimestart]),'HH:MM:SS');
            d = datetime(tmpval,'TimeZone','UTC');
            d.TimeZone = 'Australia/Brisbane';
            eegtimestart = datestr(d,'HH:MM:SS');
            tmpval = datestr(datenum([1970, 1, 1, 0, 0, eegtimeend]),'HH:MM:SS');
            d = datetime(tmpval,'TimeZone','UTC');
            d.TimeZone = 'Australia/Brisbane';
            eegtimeend = datestr(d,'HH:MM:SS');

            EEG.timestart = eegtimestart;
            EEG.timeend = eegtimeend;

            % Added by RJ 01/05/2019 % the info structure from TDTduration has
            % posixtime information inside, so I will use this for relating the video
            % to the data file.
            EEG.info = info;
            %date_format = 'yyyy:mm:dd';% '2016-07-29 10:05:24'; from posixtime Matlab website

            EEG.epoch_start = EEG.info.headerstarttime;
            EEG.epoch_end = EEG.info.headerstoptime;

            if isstruct(LFP) == 0 %Probably non-Synapse
                EEG.epoch_times = linspace(EEG.info.headerstarttime, EEG.info.headerstoptime, length(LFP));
            else
                for fielInd = 1:size(dataFiels,1)
                    thisFiel = dataFiels{fielInd};
                    EEG.(thisFiel).epoch_times = linspace(EEG.info.headerstarttime, EEG.info.headerstoptime, length(LFP.(thisFiel)));
                end
            end

            % EEG.stimdata_resamp = stimdata_resamp;

            %EEG.data(17,:) = stimdata_resamp;

            %fprintf('Block %d start time: %s...\n',chunkidx,eegtimestart);
            %fprintf('Block %d end time: %s...\n',chunkidx,eegtimeend);
            disp(['Block #',num2str(chunkidx),' start time: ',eegtimestart,'...']);
            disp(['Block #',num2str(chunkidx),' end time: ',eegtimeend,'...']);

            %#######
            %Step 4.5 - Find pol. reversal if requested for calib data (Note: Will probably crash horribly on Synapse data)
            if calibMode == 1 && findPolRevers == 1
                polReverseChan = [];
                %Find stimulus peaks (all)
                stimSepTime = 1; %Assumed time between calibration peaks (s)
                stimData = EEG.stims(1,:); %Only take row 1 of stims
                stimDataMean = nanmean(stimData);
                stimDataSD = nanstd(stimData(1:EEG.srate)); %Use 1st second of data to find std

                %##
                [stimDataHist, stimDataHistCenters] = hist(stimData,256); %Make hist
                stimDataHist(1:floor(size(stimDataHist,2)*0.75)) = NaN; %Remove lower 75% of data
                [ ~ , stimSignalConsistentPeakHeightIdx] = nanmax(stimDataHist); %Find bin position of peak of hist (Coincides with peak of sine wave)
                stimSignalConsistentPeakHeight = stimDataHistCenters(stimSignalConsistentPeakHeightIdx); %Find Y value of bin position
                %##

                %[stimPKS,stimLOCS] = findpeaks(stimData,'MinPeakHeight',stimDataMean+2*stimDataSD, 'MinPeakDistance',stimSepTime*EEG.srate*0.75);
                    %"Find peaks more than mean+2*SD in height, separated by at least 75% of a cycle"
                [stimPKS,stimLOCS] = findpeaksbase(stimData,'MinPeakHeight',stimSignalConsistentPeakHeight-0.25*stimSignalConsistentPeakHeight, 'MinPeakDistance',stimSepTime*EEG.srate*0.25);
                    %Same as above version, except MinPeakDistance is reduced to 25% to catch both start and end spikes (this will factor into following filtering)
                %[stimPKS,stimLOCS] = findpeaks(stimData,'MinPeakProminence',8*stimDataSD, 'MinPeakDistance',stimSepTime*EEG.srate*0.25);
                    %Same again, but swapped MinPeakHeight for MinPeakProminence, to account for non-zero baselines
                %Filter stimPeaks according to whether they are onset or offset
                stimLOCSProc = []; %Filtered copies of parents
                stimPKSProc = [];
                for i = 1:size(stimLOCS,2)
                    if stimData( floor(stimLOCS(i) + stimSepTime*EEG.srate*0.25) ) >= stimDataMean+4*stimDataSD
                        stimLOCSProc = [stimLOCSProc,stimLOCS(i)];
                        stimPKSProc = [stimPKSProc,stimPKS(i)];
                    end
                end
                %QA
                if isempty(stimLOCSProc) == 1 || size(stimLOCSProc,2) > 36*3 %First statement is complete miss of detection (or no stims), Second statement is potential overfind
                    ['## Warning: Critical failure in calibration stimulus peak detection ##']
                    crash = yes
                end
                %Testatory figure
                figure
                plot(stimData, 'b')
                hold on
                scatter(stimLOCS,stimPKS,'k')
                scatter(stimLOCSProc,stimPKSProc,'g')
                %title([datasetname,' - ',allsetname,' - Stim data (all peaks - black, valid peaks - green)'])
                title([datasetname,' - ',allsetname,' - Stim data (all peaks - black, valid peaks - green)'])

                %Collect channel data at all stim peaks into hyperData
                captureWindowSize = floor(stimSepTime*EEG.srate*0.4); %Collect 90% of cycle following stimulus onset

                lfpData = EEG.data;
                lfpDataProc = lfpData - nanmean(lfpData,2); %Adjust to mean (Note dimension specification)
                for chanInd = 1:size(lfpDataProc,1)
                    lfpDataProc(chanInd,:) = smooth(lfpDataProc(chanInd,:),20); %Arbitrarily smooth
                        %Note: Smoothing may affect detections
                end

                stimHyperData = [];
                for i = 1:size(stimLOCSProc,2)
                    stimHyperData(:,:,i) = lfpDataProc(:, stimLOCSProc(i):stimLOCSProc(i)+captureWindowSize );
                end

                meanHyperData = nanmean(stimHyperData(:,:,:),3);

                %Calculate local maxima within capture window (for each channel)
                tfData = [];
                firstTFs = [];
                for chanInd = 1:size(stimHyperData,1)
                    %%tfData(chanInd,:) = islocalmax(meanHyperData(chanInd,:),'MinProminence',nanmax(meanHyperData(chanInd,:))/2); %Use half max as threshold
                    tfData(chanInd,:) = islocalmax( abs(meanHyperData(chanInd,:)),'MinProminence',nanmax(meanHyperData(chanInd,:))/2); %Use half max as threshold
                    tfData(chanInd,1:10) = NaN; %Artefact detection removal
                    if nansum(tfData(chanInd,:)) >= 1
                        firstTFs(chanInd) = find(tfData(chanInd,:) == 1,1);
                    else
                        firstTFs(chanInd) = NaN; %No TF found
                    end
                end

                %Find channel values at median first TF location
                medianFirstTF = floor(nanmedian(firstTFs));
                channelVals = meanHyperData(:,medianFirstTF);
                channelValsDeflection = abs(channelVals - nanmean(meanHyperData,2)); %These values represent which channels were most deviated from their own mean during the first TF
                [~,prosPolReverseChanMin] = nanmin(channelValsDeflection);

                %Testatory to show vals
                figure
                plot(channelVals, '-or')
                hold on
                plot(channelValsDeflection, '-ob')
                title([datasetname,' - ',allsetname,' - Values of channels at first TF'])

                %QA
                if prosPolReverseChanMin >= 15 || prosPolReverseChanMin <= 8
                    ['#- Warning: Detected polarity reversal channel index (',num2str(prosPolReverseChanMin),') outside of expected range -#']
                    %crash = yes %May be overkill
                end
                %Check to see if next or prior channel was a flip
                confirmedFlip = 0;
                if prosPolReverseChanMin ~= size(stimHyperData,1) &&  prosPolReverseChanMin ~= 1
                    if ( abs(channelVals(prosPolReverseChanMin-1)) ~= channelVals(prosPolReverseChanMin-1) && abs(channelVals(prosPolReverseChanMin+1)) == channelVals(prosPolReverseChanMin+1) ) || ...
                            ( abs(channelVals(prosPolReverseChanMin-1)) == channelVals(prosPolReverseChanMin-1) && abs(channelVals(prosPolReverseChanMin+1)) ~= channelVals(prosPolReverseChanMin+1) )
                        disp(['-- Polarity flip confirmed around prospective (minima) channel ',num2str(prosPolReverseChanMin),' --'])
                            %Note: Does not check for more than one flip
                        confirmedFlip = 1;
                    else
                        ['#- Warning: No polarity flip detected around prospective (minima) channel ',num2str(prosPolReverseChanMin),' -#']
                        confirmedFlip = 0;
                    end
                end

                %Find literal polarity reversal by iterating through all channels
                prosPolReverseChanFlips = [];
                numPolFlips = 0;
                for i = 2:size(stimHyperData,1)-1
                    if ( abs(channelVals(i-1)) ~= channelVals(i-1) && abs(channelVals(i+1)) == channelVals(i+1) ) || ...
                            ( abs(channelVals(i-1)) == channelVals(i-1) && abs(channelVals(i+1)) ~= channelVals(i+1) )
                            %"Find when prior channel is negative and current is positive OR prior is positive and current is negative"
                        if confirmedFlip ~= 1
                            disp(['-- Polarity flip confirmed around channel ',num2str(i),' --'])
                        end
                        numPolFlips = numPolFlips + 1;
                        %if numPolFlips == 1
                        prosPolReverseChanFlips = [prosPolReverseChanFlips,i];
                        %end
                    else
                        if confirmedFlip ~= 1
                            disp(['#- Warning: No polarity flip detected around channel ',num2str(i),' -#'])
                        end
                    end
                end
                
                %Testatory (Shift with processing writing)
                figure
                for chanInd = 1:size(stimHyperData,1)
                    %plotData = nanmean(stimHyperData(chanInd,:,:),3);
                    plotData = meanHyperData(chanInd,:); %Original
                    %%plotData = plotData - plotData(1); %Correct start to 0
                    %%plotData = abs(plotData); %Abs as used for firstTFs calculations (Optional)
                        %Note: Not completely representative, given lack of 0 correction in firstTFs calculation
                    %Minima result
                    if chanInd == prosPolReverseChanMin
                        plot( plotData , 'LineWidth' , 2 , 'Color', 'k', 'LineStyle', ':' )
                    end
                    %Flip result
                    if isempty(prosPolReverseChanFlips) ~= 1 && chanInd == prosPolReverseChanFlips(1)
                        plot( plotData , 'LineWidth' , 2 , 'Color', 'g', 'LineStyle', '--' )
                    end
                    plot( plotData )
                    hold on
                    %pause(1)
                    %(Attempt to) Add in first TF location
                    try
                        scatter(firstTFs(chanInd),plotData( firstTFs(chanInd) ))
                        text(firstTFs(chanInd),plotData( firstTFs(chanInd) ),num2str(chanInd), 'FontSize', 9, 'Color', 'r')
                    end
                    
                end
                title([datasetname,' - ',allsetname,' - Average ERPs at stimulus onset (G-First flip, K-Minima)'])

                %See if minima method and flip method agree on polarity reversal location
                if isempty(prosPolReverseChanFlips) ~= 1
                    if prosPolReverseChanFlips(1) == prosPolReverseChanMin
                        disp(['-- Polarity reversal (minimum) and first polarity reversal (flip) agree (',num2str(prosPolReverseChanMin),') --'])
                        polReverseChan = prosPolReverseChanMin;
                    else
                        ['#- Polarity reversal (minimum) and Polarity reversals (flip) disagree; Using first flip channel (',num2str(prosPolReverseChanFlips(1)),') -#']
                        polReverseChan = prosPolReverseChanFlips(1);
                    end
                else
                    ['#- Polarity reversals (flip) list is empty; Using minima result (',num2str(prosPolReverseChanMin),') -#']
                     polReverseChan = prosPolReverseChanMin;
                end

                EEG.detectedPolReversalChan = polReverseChan;

                %Pretty plot
                figure
                maxVal = nanmax(nanmax(stimHyperData(:,:,end)));
                for chanInd = 1:size(stimHyperData,1)
                    %plotData = meanHyperData(chanInd,:); %Original
                    plotData = stimHyperData(chanInd,:,end); %Original
                    plotData = plotData - plotData(1) + (size(stimHyperData,1)-chanInd)*maxVal ; %Correct start to 0
                    %%plotData = abs(plotData); %Abs as used for firstTFs calculations (Optional)
                        %Note: Not completely representative, given lack of 0 correction in firstTFs calculation
                    %{
                    %Minima result
                    if chanInd == prosPolReverseChanMin
                        plot( plotData , 'LineWidth' , 2 , 'Color', 'k', 'LineStyle', ':' )
                    end
                    %Flip result
                    if isempty(prosPolReverseChanFlips) ~= 1 && chanInd == prosPolReverseChanFlips(1)
                        plot( plotData , 'LineWidth' , 2 , 'Color', 'g', 'LineStyle', '--' )
                    end
                    %}
                    plot( plotData, 'LineWidth', 1.25 )
                    hold on
                    %pause(1)
                    %{
                    %Add in first TF location
                    scatter(firstTFs(chanInd),plotData( firstTFs(chanInd) ))
                    text(firstTFs(chanInd),plotData( firstTFs(chanInd) ),num2str(chanInd), 'FontSize', 9, 'Color', 'r')
                    %}
                end
                title([datasetname,' - ',allsetname,' - Stacked ERPs at stimulus onset (fs: ',num2str(resamplefreq),')'])
                xlim([0,size(plotData,2)])
                %taishi

            end            
            %#######

            %%Step 5: Create EEGlab file for LfP data..
            S = [];
            %S.eeg_filename = [datasetname  '_chunk_' char(sprintf("%0.2d",chunkidx))];
            S.eeg_filename = [datasetname,'_chunk_',num2str(blockID)];
            try
                S.eeg_filename = [char(horzcat(S.eeg_filename{1:length(S.eeg_filename)}))];
            catch
                S.eeg_filename = [char(horzcat(S.eeg_filename(1:length(S.eeg_filename))))];
            end
            % S.output_folder  = outputfolder; % I had to add char(horzcat( due to the
            % way the new code finds the directory to loop through
            outputfolder = char(horzcat(outputfolder(1:length(outputfolder))));
            S.output_folder  = outputfolder;

            %fprintf(['Saving to %s.set.\n',S.eeg_filename]);
            disp(['-- Attempting to save to ',S.eeg_filename, ' --'])

            new_dir = outputfolder;
             if ~isdir(new_dir)
                 mkdir(new_dir);
             end
             save_flag = 0;

             while save_flag == 0
                 try
                     %  pop_saveset(EEG,'filepath',S.output_folder,'filename',[S.eeg_filename '.set']);
                     varDetails = whos();
                     findEgg = [];
                     for eggInd = 1:size(varDetails,1)
                         if isempty( strfind( varDetails(eggInd).name , 'EEG' ) ) ~= 1
                            findEgg = eggInd; %Technically susceptible to nested names
                            continue
                         end
                     end
                     %if whos('EEG').bytes < 1000 * 1000 * 1000 * 2
                     if varDetails(findEgg).bytes < 1000 * 1000 * 1000 * 2
                         save([new_dir filesep S.eeg_filename '.mat'],['EEG']);
                     else
                         disp(['EEG structure larger than 2GB; Saving with alternative method'])
                         save([new_dir filesep S.eeg_filename '.mat'],['EEG'], '-v7.3');
                     end
                     save_flag = 1;
                     disp('Successfully saved data.');
                 catch
                     pause(5)
                     warning('Save failed. Trying again in 5 seconds...');
                     save_flag = 0;
                 end
             end    

             clear chandata stimdata %data %data not cleared so as to be used slightly lower
             clear chandata_resamp stimdata_resamp
             clear EEG stimdata_resamp data LFP %Added more clearing to improve memory use
             
             processingEndTime = clock;
             disp(['(Chunk processed in ',num2str(etime(processingEndTime,processingStartTime)),'s)'])
             
             %Check if error occurred and iterate value if so
            if thisChunkIsError == 1
                numErrorChunks = numErrorChunks + 1;
            end
        end

        % tempval = [];
        % disp('Turning diary off...'); diary OFF;
        % 
        % if length(fopen('all')) > 5
        %     disp('Using fclose to close open file connections.');
        %     fclose('all')
        % end
        
        if numErrorChunks > 0
            ['## Caution: ',num2str(numErrorChunks),' chunks reported an error during processing ##']
        end

        
    end

end

disp(['Completed extraction of selected data'])

