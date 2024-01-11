function BLOCKS = collateEphysData(fly_record,chosenOnes,focusPeak,timeFrequency,homeDirectory,aux_plots)
%collateEphysData Summary of this function goes here
%   Detailed explanation goes here
    BLOCKS = struct;

    % loop though all blocks, irrespective of experiment
    for b = find(chosenOnes)

        date = datestr(fly_record.Date(b),'ddmmyy');
        block = num2str(fly_record.Block(b));

        % load this block's data
        load([homeDirectory '/SEoutput/' date '/LFP/Analyzed_TagTrials_block' block '/' date '_chunk_0'],'EEG');

        % photodiode and lfp data
        LFP = EEG.LFP1.data(fly_record.LFPChannel(b),:);
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
        if ~contains(fly_record.Comments(b),'dark','IgnoreCase',true)
            PHOT(2,:) = -PHOT(2,:);
            rawPHOT(2,:) = -rawPHOT(2,:);
        else
            PHOT(1,:) = -PHOT(1,:); 
            rawPHOT(1,:) = -rawPHOT(1,:);
        end

    %     PHOT(3,:) = -PHOT(3,:);

        % butterworth and/or notch and/or savitsky-golay
        % check parameters inside
        LFP = filterLFP(LFP, resampleFreq);

        %sometimes the photodiode has hiccups
        PHOT = trim_phot_outliers(PHOT, 12);

        % trim LFP (remove anything from LFP beyond some sd set in fly_record)
        LFP = trimLFP(LFP,fly_record.LFPsd(b),aux_plots);

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

    end
end

