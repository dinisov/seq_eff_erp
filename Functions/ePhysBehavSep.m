function [PHOT] = ePhysBehavSep(phot, startTime, endTime, fly, block, date, state, homeDirectory, timeThreshold, shuffleMode) % and more things (vidDirectory?)
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
    % fly = 23;
    % block = 1;
    % state = 1;
    % startTime = EEG.timestart;
    % endTime = EEG.timeend;
    % date = EEG.filename;
    % % homeDirectory = '../../Mae 2026';
    % % homeDirectory = 'I:\BVS2026TWCF-Q9201\Mae 2026';
    % homeDirectory = 'C:\Users\uqaspen5\ePhys';
    % timeThreshold = 20;

    %shuffleMode - 0 - Do nothing, 1 - Use random positions within opposite state pair matches, 2 - Use first X seconds based on pair match length, 3 - As with 2, but use last X seconds

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
    % vidDirectory = [homeDirectory, '\Videos\Fly Videos\Mae'];

    posixStart = posixtime(datetime([char(date), startTime], 'Format', 'ddMMyyHH:mm:ss')); %, 'TimeZone', '+10')); % did you know the capital H is for military time, lowercase h is 12-hour
    posixEnd = posixtime(datetime([char(date), endTime], 'Format', 'ddMMyyHH:mm:ss')); %, 'TimeZone', '+10'));
    photTimes = linspace(posixStart, posixEnd, size(phot, 2)); %maybe do this outside of function??

    
    %copying behavProcess loading

    thisCSVs = dir( [vidDirectory,filesep,'fly',num2str(fly),'*mov.csv'] ); %returns structure
    %QA
    if isempty(thisCSVs)
        ['-# No data found for ',num2str(fly),' #-']
        failCount = failCount + 1;
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
        bwInacBinary = bwlabel(~inacBinaryProc); %Empirically confirmed correct orientation
    else %state == 0, active
        bwInacBinary = bwlabel(inacBinaryProc); %Empirically confirmed correct orientation
    end

    % behavData informative plot
    coords = find(inacBinaryProc == state); %Confirmed (empirically) correct
    temp = nan(size(inacBinaryProc));
    temp(coords) = inacBinaryProc(coords);
    figure
    plot(behavPosix - photTimes(1), movData)
    hold on
    plot(behavPosix - photTimes(1), repmat(acThreshLevel, size(behavPosix, 1), size(behavPosix, 2)), 'linewidth', 1.5, 'color', 'k', 'linestyle', '--')
    % line([behavPosix(1), behavPosix(2)], [acThreshLevel, acThreshLevel], 'linewidth', 1.5, 'color', 'k', 'linestyle', '--')
    bYs = get(gca, 'ylim');
    plot(behavPosix - photTimes(1), inacBinaryProc*(bYs(2)/2) + bYs(2), 'linewidth', 2)
    plot(behavPosix - photTimes(1), temp*(bYs(2)/2) + bYs(2), 'Color','r', 'linewidth', 3)
    xlim([0, photTimes(end) - photTimes(1)])
    xlabel('time(seconds)')
    title(['Behavioural Data for fly ', num2str(fly), ' block ', num2str(block), ', timeThres:', num2str(timeThreshold),', Target state:',num2str(state)])
    set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12, 'box', 'off')
    %legend({'behavData', 'activityThreshold', 'inacBinary'})
    legend({'behavData', 'activityThreshold', 'inacBinary','target state'})
    set(gcf, 'name', ['fly ', num2str(fly), ' block ', num2str(block), ', separated photodiode'])  

    %Shuffle, if requested
    if shuffleMode
        bwInacBinaryActive = bwInacBinary; %For editing
        temp = bwlabel( ~bwInacBinary ); %Technically this is bwlabelling a bwlabel product, but empirically this is fine due to how booleans are performed in MATLAB
            %This is just a bwlabel of the inverse (So, activity when looking for inactivity and vice versa)

        %Test figure for proof of inversion
        %{{
        figure
        plot(temp-25)
        hold on
        plot(bwInacBinary+1.2)
        plot(~bwInacBinary)
        xlim([1,length(bwInacBinary)])
        if state == 1
            legend({'Inverted bwlabel [Inactivity (Target)]','Original bwlabel [Labels Activity (Not target)]','Inversion'})
        else
            legend({'Inverted bwlabel [Activity (Target)]','Original bwlabel [Labels Inactivity (Not target)]','Inversion'})
        end
        %}

        %Begin shuffling
        if shuffleMode == 1 || shuffleMode == 2 || shuffleMode == 3 %Random pos, First, Last
            %for gap = 1:nanmax(bwInacBinary) %Reminder: bwInacBinary labels the inverse of the target state
            for gap = 1:nanmax(temp) %Reminder: bwInacBinary labels the inverse of the target state

                %if gap == 18
                %fuentego
                %end

                %origCoords = find( bwInacBinary == gap );
                origCoords = find( temp == gap ); %Swap who looking at
                %matchInv = nansum( temp == gap ); %How large is matched partner (e.g. If looking at first activity bout, how large is first inactivity bout)
                matchInv = nansum( bwInacBinary == gap ); %How large is matched partner (e.g. If looking at first activity bout, how large is first inactivity bout)
                %Trim if match larger than original
                if matchInv > numel(origCoords)
                    disp(['Gap #',num2str(gap),' original length shorter than match partner (',num2str(numel(origCoords)),' vs ',num2str(matchInv),')'])
                    matchInv = numel(origCoords);
                end
                if shuffleMode == 1 %Random
                    pos = randi( numel(origCoords)-matchInv+1 ); %Find a random position at least matchInv from end of original gap
                        %May act funny if numel(origCoords) == matchInv?
                    origCoords = origCoords( pos:pos+matchInv-1 );
                elseif shuffleMode == 2 %First
                    origCoords = origCoords( 1:matchInv );                                   
                elseif shuffleMode == 3 %Last
                    origCoords = origCoords( end-matchInv+1:end );
                end

                %subtractCoords = setdiff( find( bwInacBinary == gap ), origCoords ); %What coords to flatten now
                subtractCoords = setdiff( find( temp == gap ), origCoords ); %What coords to flatten now

                %{
                need to leave bwInacBinary 'intact' but ablate extra according to origCoords
                    ablate = add more coords to bwInacBinary? (because used for flattening later?)
                %}

                %bwInacBinaryActive( subtractCoords ) = 0;
                %bwInacBinaryActive( subtractCoords ) = 1; %NOT CORRECT

                if shuffleMode == 1 %Random
                    bwInacBinaryActive( subtractCoords( find( subtractCoords < origCoords(1) ) ) ) = gap; %Preceding elements
                    bwInacBinaryActive( subtractCoords( find( subtractCoords > origCoords(end) ) ) ) = gap + 1; %Postceding elements
                elseif shuffleMode == 2 %First
                    bwInacBinaryActive( subtractCoords ) = gap + 1; %MAY CAUSE ISSUES IF LAST GAP?                                  
                elseif shuffleMode == 3 %Last
                    bwInacBinaryActive( subtractCoords ) = gap; %subtractCoords here relate to start of bout?
                end


                %Reminder: bwInacBinary will be iterated along later and all labelled regions will be flattened

            end

            %Testatory figure
            figure
            subplot(2,1,1)
            plot(bwInacBinary,'Color','b')
            hold on
            plot( inacBinaryProc- 1.2,'Color','c')
            ylim([-2,nanmax(bwInacBinaryActive)+1])
            xlim([1,numel(bwInacBinaryActive)])
            title(['Original ac/inac binary (Mode ',num2str(shuffleMode),')'])
            legend({'Original target binary'})
            %hold on
            subplot(2,1,2)
            plot(bwInacBinaryActive,'Color','red')
            hold on
            plot( inacBinaryProc- 1.2,'Color','c')
            ylim([-2,nanmax(bwInacBinaryActive)+1])
            xlim([1,numel(bwInacBinaryActive)])
            title(['Adjusted ac/inac binary (Mode ',num2str(shuffleMode),')'])
            legend({'Modified target binary','Raw ac/inac binary'})

        end

        bwInacBinary = bwInacBinaryActive;

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

    % photodiode and binary
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
    %plot(behavPosix - photTimes(1), inacBinaryProc*(ys(2)/2) + ys(2), 'LineWidth', 2) %Andre method
    if shuffleMode
        plot(behavPosix - photTimes(1), (inacBinaryProc*0.5)*(ys(2)/2) + ys(2), 'LineWidth', 2) %Plot original
    end
    plot(behavPosix - photTimes(1), ~[bwInacBinary == 0]*(ys(2)/2) + ys(2), 'LineWidth', 2) %Matt method, uses inversion + boolean to more accurately represent data that was used to flatten phot
    xlim([0, photTimes(end) - photTimes(1)])
    xlabel('time(seconds)')
    title(['separated photodiode for fly ', num2str(fly), ' block ', num2str(block)])
    set(gca, 'LineWidth', 2, 'TickDir', 'out', 'FontWeight', 'bold', 'FontSize', 12, 'box', 'off')
    set(gcf, 'name', ['fly ', num2str(fly), ' block ', num2str(block), ', separated photodiode'])
    
    % resaving photData/returning phot
    PHOT = phot;

end