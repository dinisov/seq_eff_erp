function [PHOT] = ePhysBehavSep(phot, startTime, endTime, fly, block, date, state, homeDirectory, timeThreshold) % and more things (vidDirectory?)
    % behavioural separation function for ePhys data
    %   uses pixel subtraction method to identify periods of inactivity.
    %   removes photodiode activity during these periods using time as anchor
    %
    %   phot: photodiode date [channels x data]
    %   ###photTimes = posix time associated with each element of phot's length (EEG.epoch_times)###
    %   startTime = time of exp start (for phot) [EEG.timestart]
    %   endTime = 
    %   fly: int fly number
    %   block: int block number
    %   date: int date for current fly (in format of ddMMyy)
    %   state: whether to remove photodiode at active (0) or inactive (1) periods

    %presets for testing
    % phot = EEG.PHOT.data;
    % fly = 2;
    % block = 1;
    % state = 1;
    % startTime = EEG.timestart;
    % endTime = EEG.timeend;
    % date = EEG.filename;
    % % homeDirectory = '../../Mae 2026';
    % % homeDirectory = 'I:\BVS2026TWCF-Q9201\Mae 2026';
    % homeDirectory = 'C:\Users\uqaspen5\ePhys';
    % timeThreshold = 20;

    % disp("'WARNING behavSep function currently hardcoded for Mae's file architecture")
    if state == -1
        disp('No behavioural separation conducted because behavState set to -1')
        PHOT = phot;
        return
    end

    disp(['loading files for behavioural separation: fly ', num2str(fly), ' block ', num2str(block)])
    % will need to put in an intelligent method for loading behavData below
    % behavData = readtable('\\data.qbi.uq.edu.au\RFDG2021-Q4413\Andre\ePhys_Data\Fly2\fly2_30_10_25_01__mov.csv');
    date = datetime(date, 'Format', 'ddMMyy');
    flyVidName = ['fly', num2str(fly), char(datetime(date, 'Format', 'dd_MM_yy')), '_'];
    % vidDirectory = 'I:\BVS2026TWCF-Q9201\Mae 2026\Videos\Fly Videos'; %eventually make dynamic (maybe prereq for function?)
    vidDirectory = [homeDirectory, '\Videos'];

    posixStart = posixtime(datetime([char(date), startTime], 'Format', 'ddMMyyHH:mm:ss')); %, 'TimeZone', '+10')); % did you know the capital H is for military time, lowercase h is 12-hour
    posixEnd = posixtime(datetime([char(date), endTime], 'Format', 'ddMMyyHH:mm:ss')); %, 'TimeZone', '+10'));
    photTimes = linspace(posixStart, posixEnd, size(phot, 2)); %maybe do this outside of function??

    
    %copying behavProcess loading

    thisCSVs = dir( [vidDirectory,filesep,'fly',num2str(fly),'*mov.csv'] ); %returns structure
    %QA
    failCount = 0;
    if isempty(thisCSVs)
        ['-# No data found for ',num2str(fly),' #-']
        disp([[vidDirectory,filesep,'fly',num2str(fly),'*mov.csv'] ])
        failCount = failCount + 1;
        crash = yes
    end

    %Collate
    behavData = [];
    sepData = [];
    for IIDN = 1:size(thisCSVs,1)
        %allData = [allData; readtable([thisCSVs(IIDN).folder,filesep,thisCSVs(IIDN).name]) ];
        sepData{IIDN} = readtable([thisCSVs(IIDN).folder,filesep,thisCSVs(IIDN).name]);
        behavData = [behavData;...
            readtable([thisCSVs(IIDN).folder,filesep,thisCSVs(IIDN).name]), ...
            array2table( repmat(IIDN, size(sepData{IIDN},1), 1) ,'VariableNames',{'CSVnum'})];
    end
    % behavData = readtable("I:\BVS2026TWCF-Q9201\Mae 2026\Videos\Fly Videos\fly1_29_10_25_01_.avi");
    % "I:\BVS2026TWCF-Q9201\Mae 2026\Videos\Fly Videos\fly1_29_10_25_01_.avi"
    

    %important variables (consider moving to options in function)
    acMeanThresh = +2;
    minAcTime = 20;
    minSleepBoutTime = timeThreshold;
    
    %% converting behavCSVs into posix
    % currently working under assumption of one csv per fly (may need to adjust later - take inspo from behavProcess)
    %getting posix times for behavData (from behavProcess)
    behavDates = num2str([behavData.Year,behavData.Month,behavData.Date,behavData.Hour,behavData.Mins,behavData.Seconds,behavData.usec]);
    %Clear up potentially missing proximal zeroes
    for wU = 5:-1:0 %uWu
        zeroInds = find(behavDates(:,end-wU) == ' ');
        behavDates(zeroInds,end-wU) = '0';
    end
    behavDatetimes = datetime(behavDates,'Format', 'yyyy      MM      dd      HH      mm      ss  SSSSSS');
    
    behavPosix = posixtime(behavDatetimes); %posix times for elements of behavData (used for synchronisation)
    
    %% behavioural separation
    % (from behavProcess)
    behavFrameRate = 1 / nanmedian(diff([behavData.Seconds+(behavData.usec/1000000)]));
    %Bootleg derive the framerate by finding the median time difference between frames
    
    disp(['-- commencing behavioural separation for fly ', num2str(fly), ' Block ', num2str(block), ' --']) %current fly needs to be fixed
    %Calculate sleep/wake
    movData = behavData.avCntrSize;
    movMean = nanmean(movData);
    movSD = nanstd(movData);
    
    acThreshLevel = movMean + acMeanThresh*movSD;
    
    acUpper = movData > acThreshLevel;
    
    %New, BWLabel based method
    %minActivityTime = minAcTime * BaseFrameRate;
    minActivityTime = minAcTime * behavFrameRate;
    
    %tempUpperBinary = isnan( acUpper ) ~= 1;
    %tempUpperBW = bwlabel( tempUpperBinary );
    tempUpperBW = bwlabel( ~acUpper ); %Label gaps between ac pockets
    %invTempUpperBW = bwlabel( acUpper ); %Label ac pockets
    
    temp = nansum( tempUpperBW == [1:nanmax(tempUpperBW)] , 1); %Find sizes of all gaps
    %temp2 = nansum( invTempUpperBW == [1:nanmax(invTempUpperBW)] , 1); %Find sizes of all pockets
    
    %temp2 = nansum(tempUpperBW == find(temp <= (minAcTime*behavFrameRate)),2);
    
    acUpperProc = acUpper;
    acUpperProc( nansum(tempUpperBW == find(temp <= (minAcTime*behavFrameRate)),2) == 1 ) = 1;
    %Find all gaps separated by <minAcTime, set to 1 in binary data
    
    inacBinaryProc = ~acUpperProc;
    inacBinaryBW = bwlabel( inacBinaryProc );
    temp = nansum( inacBinaryBW == [1:nanmax(inacBinaryBW)] , 1); %Find sizes of all inacBouts
    disp(['There are ',num2str(nanmax(inacBinaryBW)),' prospective inactivity bouts (Av: ',num2str(nanmean(temp/behavFrameRate)),'s)'])
    
    inacBinaryProc( nansum(inacBinaryBW == find(temp < (minSleepBoutTime*behavFrameRate)),2) == 1 ) = 0;
    %Find all prospective sleep bouts < minimum size, flatten
    inacBinaryBW = bwlabel( inacBinaryProc ); %Recalculate
    
    %reporting number of bouts uncovered
    holeSizes = nansum( inacBinaryBW == [1:nanmax(inacBinaryBW)] , 1);
    disp(['Fly ', num2str(fly), ' Block ', num2str(block), ': There were ',num2str(length(holeSizes)),' sleep bouts meeting criteria (Min ',num2str(minSleepBoutTime),'s)']) 
    
    
    %% synchronisation
    % using posix times to find indices of inactivity for phot data
    % NOTE: inac doesn't alway mean inac below this line (will depend on state)
    if state == 1 %inactive
        bwInacBinary = bwlabel(~inacBinaryProc);
    else %state == 0, active
        bwInacBinary = bwlabel(inacBinaryProc);
    end
    inacTimes = nan(max(bwInacBinary), 2); % coordinates for start and end of activity gaps
    for gap = 1:max(bwInacBinary)
        inacTimes(gap, 1) = behavPosix(find(bwInacBinary == gap, 1, 'first')); %start
        inacTimes(gap, 2) = behavPosix(find(bwInacBinary == gap, 1, 'last')); %end
    end
    
    % ablating photData in inactive periods
    wholePhotAvg = nanmean(phot, 2);
    for gap = 1:size(inacTimes, 1)
        inacPhotIndex = find(photTimes > inacTimes(gap,1) & photTimes < inacTimes(gap,2));
        if isempty(inacPhotIndex) % behavioural data outside of experiment time
            continue
        end
    
        % using mean of wholeExpAvg to ablate phot spikes
        for chan = 1:size(phot, 1)
            phot(chan, inacPhotIndex) = wholePhotAvg(chan);
        end
    end
    
    figure
    plot(photTimes - photTimes(1), phot(1,:))
    % plot(photTimes, phot(1,:))
    hold on
    for chan = 1:size(phot, 1)
        plot(photTimes - photTimes(1), phot(chan,:))
    end
    % plot(photTimes - photTimes(1), phot(2,:))
    % plot(photTimes, phot(2,:))
    ys = get(gca, 'ylim');
    % plot(behavPosix, inacBinaryProc*(ys(2)/2) + ys(2), 'LineWidth', 2)    
    plot(behavPosix - photTimes(1), inacBinaryProc*(ys(2)/2) + ys(2), 'LineWidth', 2)
    xlim([0, photTimes(end) - photTimes(1)])
    xlabel('time(posix)')
    title(['separated photodiode for fly ', num2str(fly), ' block ', num2str(block)])

    % resaving photData/returning phot
    PHOT = phot;

end