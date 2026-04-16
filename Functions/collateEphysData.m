function BLOCKS = collateEphysData(fly_record,chosenOnes,focusPeak,timeFrequency,homeDirectory,aux_plots, zeroShiftMode, signalInversion,overrideChannel,rawDataPlot)
%collateEphysData Summary of this function goes here
%   Detailed explanation goes here
    if exist('overrideChannel') && ~isempty(overrideChannel)
        disp(['LFP channel #',num2str(overrideChannel),' will be used for all analysis'])
    end
    if exist('signalInversion') && isequal(signalInversion,'-')
        disp(['LFP data will be inverted'])
    end

    BLOCKS = struct;

    %zeroShiftNotified = 0;
    % loop though all blocks, irrespective of experiment
    for b = find(chosenOnes)

        fly = fly_record.Fly(b);
        date = datestr(fly_record.Date(b),'ddmmyy');
        block = num2str(fly_record.Block(b));

        % load this block's data
        %[homeDirectory '/SEoutput/' date '/LFP/Analyzed_TagTrials_Block' block '/' date '_chunk_0']
        load([homeDirectory '/SEoutput/' date '/LFP/Analyzed_TagTrials_Block' block '/' date '_chunk_0'],'EEG');
        disp(['Loaded ',homeDirectory '/SEoutput/' date '/LFP/Analyzed_TagTrials_Block' block '/' date '_chunk_0'])

        if ~exist('overrideChannel') || isempty(overrideChannel)
            chanToUse = fly_record.LFPChannel(b);
        else
            chanToUse = overrideChannel;
        end
        %QA
        if chanToUse < 1 || chanToUse > size(EEG.LFP1.data,1)
            ['## Alert: Requested channel to use smaller than zero or larger than data ##']
            crash = yes
        end

        % photodiode and lfp data
        %LFP = -EEG.LFP1.data(fly_record.LFPChannel(b),:);
        if ~exist('signalInversion') || isequal(signalInversion,'+')
            LFP = EEG.LFP1.data(chanToUse,:);
        else
            LFP = -EEG.LFP1.data(chanToUse,:);
        end
        PHOT = EEG.PHOT.data;
        rawPHOT = EEG.PHOT.data; %preserve raw unfiltered photodiode data
        TIMES = EEG.times;
        resampleFreq = EEG.srate;

        if rawDataPlot
            %Plot raw data
                %Note: This plot is likely not good for MATLAB speed/memory
            legNames = {};
            figure
            %plot(LFP(1,:),'Color','k')
            if size(PHOT,2) == size(TIMES,2) %See below
                plot(TIMES, LFP(1,:),'Color','k')
            else
                plot(LFP(1,:),'Color','k')
            end
            legNames = [legNames,{'LFP'}];
            temp = range(LFP(1,:));
            hold on
            for pInd = 1:size(PHOT,1)
                %plot(normalize( PHOT(pInd,:) , 'Range', [-0.1*temp, 0.1*temp] ) - 0.5*temp - pInd*0.1*temp)
                if size(PHOT,2) == size(TIMES,2)
                    plot(TIMES, normalize( PHOT(pInd,:) , 'Range', [-0.1*temp, 0.1*temp] ) - 0.5*temp - pInd*0.1*temp)
                    xlabel('Time (s)')
                else
                    plot(normalize( PHOT(pInd,:) , 'Range', [-0.1*temp, 0.1*temp] ) - 0.5*temp - pInd*0.1*temp)
                    if pInd == 1
                    disp(['-# Caution: ',date,' B',block,' size of PHOT (',num2str(size(PHOT,2)),') and TIMES (',num2str(size(TIMES,2)),') data different #-'])
                    end
                end
                legNames = [legNames,{['PHOT',num2str(pInd)]}];
            end
            legend(legNames)
            
            title(['Raw data plot - ',date,' B',block,' (Fly #',num2str(fly),')'])
        end


        dataIsMulti = 0;
        if isfield(EEG,'recordingType') && isequal(EEG.recordingType,'multi')
            dataIsMulti = 1;
            disp(['Data detected as multichannel'])
        end

        % remove first/last n_sec seconds of recording
        n_sec = 1;
        chunk_a_time = n_sec*resampleFreq;
        LFP = LFP(:,chunk_a_time:end-chunk_a_time);
        PHOT = PHOT(:,chunk_a_time:end-chunk_a_time);
        rawPHOT = rawPHOT(:,chunk_a_time:end-chunk_a_time);
        TIMES = TIMES(:,chunk_a_time:end-chunk_a_time);

        %%

        % correct for the fact that the left photodiode is inverted
        % such that peaks are always upward for peak detection
        if ~dataIsMulti %Synapse data most likely
            photOneChannel = 1;
            photTwoChannel = 2;
        else %Multichannel data, probably processed like Synapse data
            photOneChannel = 1;
            photTwoChannel = 4;
            PHOT = [PHOT(photOneChannel,:);...
                PHOT(photTwoChannel,:)]; %Alter phot structure to be 2x1, as Synapse data presumably is
            disp(['(Altered phot channels used for multi analysis)'])
        end
        if ~dataIsMulti %Don't do this inversion with multichannel data (for now)
        if ~contains(fly_record.Comments(b),'dark','IgnoreCase',true)
            PHOT(2,:) = -PHOT(2,:);
            rawPHOT(2,:) = -rawPHOT(2,:);
        else
            PHOT(1,:) = -PHOT(1,:); 
            rawPHOT(1,:) = -rawPHOT(1,:);
        end
        end

        if rawDataPlot
            %BOOTLEG TRANSIENT [P]RECAPTURE
                %Note: Adds at least ~10s to time (per fly-block)
            smoothSpan = 2*(1/60)*resampleFreq; %Smooth phot data for span of at least 2 refresh cycles
                %Note: Assumes 60Hz refresh rate
            capInds = [];
            figure
            for pInd = [photOneChannel,photTwoChannel]
                temp = smooth( PHOT(pInd,:) , smoothSpan );
                temp = highpass( temp, 0.1, resampleFreq );
                tempMed = nanmedian(temp);
                tempSD = nanstd( temp );
                tempInds = find(temp > tempMed + tempSD);
                capInds = [capInds,tempInds'];
    
                subplot(1,2,pInd)
                plot(temp)
                hold on
                line([0,size(temp,1)],[tempMed+tempSD,tempMed+tempSD])
                title(['PHOT chan. ',num2str(pInd),' ', date,' B',block,' (Fly #',num2str(fly),')'])
            end

            capInds = unique(capInds); %Because (theoretically) combining multiple phot channels
            temp = zeros(1,size(PHOT,2));
            temp(capInds) = 1;
            tempBW = bwlabel(temp);
            tempEventNum = nanmax(tempBW);
            disp(['<< Bootleg event finding: ',num2str(tempEventNum),' event/s >>'])
            tempSize = nan(1,tempEventNum);
            for tempI = 1:tempEventNum
                tempSize(tempI) = nansum(tempBW == tempI);
            end
            tempMax = nanmax(tempSize);

            tempTrans = nan(tempEventNum,tempMax);
            tempPhot = nan(tempEventNum,tempMax,size(PHOT,1));
            for tempI = 1:tempEventNum
                temp2 = LFP( tempBW == tempI );
                tempTrans(tempI, [1:size(temp2,2)] ) = temp2;
                for pInd = [photOneChannel,photTwoChannel]
                    temp2Phot = PHOT(pInd, tempBW == tempI );
                    tempPhot(tempI, [1:size(temp2,2)] , pInd ) = temp2Phot;
                end
            end

            tempMean = nanmean(tempTrans,1);
            tempSEM = nanstd(tempTrans,[],1) ./ sqrt( nansum( ~isnan(tempTrans) , 1 ) );
            tempPhotMean = nanmean(tempPhot,[1,3]);
            tempPhotSEM = nanstd(tempPhot,[],[1,3]) ./ sqrt( nansum( ~isnan(tempPhot) , [1,3] ) );
                %Note: SEM difficult to use in following, because of normalization

            figure
            %plot(tempMean);
            errorbar(tempMean,tempSEM,'Color','k')
            temp = range(tempMean);
            hold on
            %for pInd = [photOneChannel,photTwoChannel]
            plot(normalize( tempPhotMean , 'Range', [-0.1*temp, 0.1*temp] ) - 0.5*temp - pInd*0.1*temp)
                %legNames = [legNames,{['PHOT',num2str(pInd)]}];
            %end
            title(['Bootleg transient acquisition plot- ',date,' B',block,' (Fly #',num2str(fly),')'])
            legend({'LFP','PHOT'})        
            xlabel('Time (frames)') %Too lazy atm to acquire/process timebase for all bootleg transient acquisition
        end

    %     PHOT(3,:) = -PHOT(3,:);

        % butterworth and/or notch and/or savitsky-golay
        % check parameters inside
        LFP = filterLFP(LFP, resampleFreq);

        %sometimes the photodiode has hiccups
        PHOT = trim_phot_outliers(PHOT, 12);

        % trim LFP (remove anything from LFP beyond some sd set in fly_record)
        LFP = trimLFP(LFP,fly_record.LFPsd(b),aux_plots);

        %Shift if requested
        if zeroShiftMode %&& zeroShiftNotified ~= 1
            LFP = LFP + abs(nanmin(LFP));
            disp(['Zero shift applied to LFP data'])
            %zeroShiftNotified = 1;
        end

        %% 

    %     figure; plot(rawPHOT(2,:)); hold on; plot(PHOT(2,:));

    %     figure; plot(PHOT(1,:)); hold on; plot(xlim, [photSD photSD]); plot(xlim, [-photSD -photSD]);
    %     figure; plot(PHOT(2,:)); hold on; plot(xlim, [photSD photSD]); plot(xlim, [-photSD -photSD]);

        BLOCKS(b).LFP = LFP;
        BLOCKS(b).PHOT = PHOT;
        BLOCKS(b).rawPHOT = rawPHOT;

        %BLOCKS(b).times = EEG.times;
        BLOCKS(b).times = TIMES; %Use equivalently processed times structure
        BLOCKS(b).resampleFreq = resampleFreq;

        % only for block experiments
        BLOCKS(b).interBlockPeriod = fly_record.InterBlockPeriod(b);
        BLOCKS(b).focusPeak = focusPeak;
        BLOCKS(b).blockLength = fly_record.BlockLength(b);

        % if this is a block experiment and the focus peak is the last in the
        % train, disregard the window setting and look at the entire period
        % between blocks (both for time-domain and time-frequency analyses)
        if ~isnan(fly_record.BlockLength(b)) && (focusPeak == BLOCKS(b).blockLength) && timeFrequency
            BLOCKS(b).window = [0 fly_record.InterBlockPeriod(b)];
        else
            BLOCKS(b).window = [-fly_record.SDT(b).*fly_record.Window1(b) fly_record.ISI(b).*fly_record.Window2(b)];
        end

        BLOCKS(b).ISI = fly_record.ISI(b) + fly_record.SDT(b);
    %     BLOCKS(b).SDT = fly_record.SDT(b);
        BLOCKS(b).peakThreshold = fly_record.Threshold(b);

        % whether to do a 1 or 2 photodiode analysis
        % 2 is better currently but 1 is more universal/convenient
        BLOCKS(b).PHOTType = fly_record.PHOTType(b);

        BLOCKS(b).dataIsMulti = dataIsMulti;

        %Transect stuff
        if ~isempty(fly_record.TransectTime(b)) && ~isnan(fly_record.TransectTime(b))
            BLOCKS(b).transectTime = fly_record.TransectTime(b); %For transect analysis/s
        else
            disp(['-# Caution: ',date,' block ',num2str(block),' lacked specified transect time; Using default 180'])
            BLOCKS(b).transectTime = 180;
        end
        if ~isempty(fly_record.AvTransectWindow1(b)) && ~isnan(fly_record.AvTransectWindow1(b))
            BLOCKS(b).avTransectWindow = [fly_record.AvTransectWindow1(b),fly_record.AvTransectWindow2(b)]; %Use frame timescale as above
        else
            BLOCKS(b).avTransectWindow = []; %Perform no average transect calcs
        end


        BLOCKS(b).fly = fly;
        BLOCKS(b).date = date;
        BLOCKS(b).block = block;

    end
end

