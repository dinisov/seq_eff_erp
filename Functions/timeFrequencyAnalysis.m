% analyse in which frequency bands SLRP and LRPR live
function timeFrequencyAnalysis(FLIES, chosenFlies, homeDirectory, plot_individual_flies, plotComponents, customComparisons, n_back, options)

arguments
    FLIES struct
    chosenFlies double
    homeDirectory char
    plot_individual_flies double
    plotComponents double
    customComparisons double
    n_back double
    options.isoMode double = 0 %Whether to analyse time frequency for isomers (1) or merged (0)
end

isoMode = options.isoMode;

load('slrp_lrpr.mat','slrp','lrpr','weird');

% gets rid of empty structures
FLIES = FLIES(~cellfun(@isempty,{FLIES.allERPs}));

resultsDirectory = [homeDirectory '\Results_Time_Frequency\'];
if exist(resultsDirectory) ~= 7
    mkdir(resultsDirectory)
    disp(['Time/Freq results directory had to be made'])
    disp(resultsDirectory)
end

%Matthew system for turning Dinis X labels into something conveniently usable
switch n_back
    case 5
        load('binomial_x_labels_latex_alt_rep.mat','binomial_x_labels_latex');
        labels = binomial_x_labels_latex;
    otherwise
        loadName = [num2str(n_back),'-back_legend.mat'];
        eval(['load ',loadName])
        %labels = anynomial_x_labels_latex; %Old style; Native ordering
        labels = anynomial_x_labels_latex_canonical; %Matches what is applied by seq_eff_order in analyseSequentialEffects 
end
%this is just to help turn horizontal sequences into vertical ones
%ind_horiz = sub2ind(size(binomial_x_labels_latex{1}),1:4,[1 1 1 5]); %Hardcoded n-back of 5
ind_horiz = sub2ind(size(labels{1}),1:n_back-1,[ones(1,n_back-2) 5]); %Dynamic
exLabels = [];
for s = 1:size(labels,2)
    %exLabels{s} = binomial_x_labels_latex{s}(ind_horiz);
    exLabels{s} = labels{s}(ind_horiz);
end

% individualFlies = false;

n_flies = length(FLIES);

time_bounds = FLIES(1).window*1000;

% options = optimset('Algorithm','interior-point','FinDiffType','central');

f_limits = [0 300];

f = FLIES(1).f;

% f_reduced = f(f > f_limits(1) & f < f_limits(2));

%for plots
y_ticks = 1:2:length(f);
y_tick_labels = f(y_ticks);

%for preallocation
fSize = length(f(f > f_limits(1) & f < f_limits(2)));
tSize = size(FLIES(1).magnitudeSEs,2);
%Matthew QA
temp = [];
for fly = 1:size(FLIES,2)
    temp(fly) = size(FLIES(1).magnitudeSEs,3);
end
if size(unique(temp),2) > 1
    ['## Alert: Dissimilar window sizes used for multi-fly time/freq analysis; Cannot proceed ##']
    return
end
zSize = unique(temp); %If we've reached here, this is a singular value

%% calculate r and p values (r plots)

if n_back == 5
    r_slrp = zeros(fSize,tSize,n_flies);
    r_lrpr = zeros(fSize,tSize,n_flies);
    r_weird = zeros(fSize,tSize,n_flies);
    r_td = zeros(fSize,tSize,n_flies);
    
    p_slrp = zeros(fSize,tSize,n_flies);
    p_lrpr = zeros(fSize,tSize,n_flies);
    p_weird = zeros(fSize,tSize,n_flies);
    p_td = zeros(fSize,tSize,n_flies);
end

%Prepare isomer stuff (or not)
if isoMode
    %QA
    if ~isfield(FLIES(1),'ISOMER')
        ['-# Isomer analysis for frequency requested, but no isomer data exists #-']
        crash = yes
    end
    isoList = fieldnames(FLIES(1).ISOMER);
    isoN = size(isoList,1);
    disp(['-- Time/Frequency will be analysed for ',num2str(isoN),' isomers --'])
    disp(isoList)
else
    isoN = 1;
    disp(['--Time/Frequency will be analysed for merged data --'])
end

for iso = 1:isoN %Do even if using merged, for reasons of simplicity
    %Quick preparation if isomer
    if isoMode
        thisIsom = isoList{iso}; %Will only be enacted in isomer mode, so no risk when isoN == 1 but isoMode == 0
        figAdd = ['-isomer',thisIsom];
    else
        figAdd = ['-merged'];
    end
    for fly = 1:n_flies
    
        thisFly = chosenFlies(fly);
        
        if ~isoMode
            td_profile = FLIES(fly).PROFILE.amplitude;
            magSEs = FLIES(fly).magnitudeSEs; 
            magSEs = magSEs(f > f_limits(1) & f < f_limits(2),:,:);
        else
            td_profile = FLIES(fly).ISOMER.(thisIsom).PROFILE.amplitude;
            magSEs = FLIES(fly).ISOMER.(thisIsom).magnitudeSEs; 
            magSEs = magSEs(f > f_limits(1) & f < f_limits(2),:,:);
        end
       
        %magSEs = FLIES(fly).magnitudeSEs; magSEs = magSEs(f > f_limits(1) & f < f_limits(2),:,:); %Moved slightly above
    
        f_size = size(magSEs,1); t_size = size(magSEs,3);
    
        if n_back == 5
            parfor f_ind = 1:f_size% frequency index
        
                for t_ind = 1:t_size% time index
                    
                   seq_eff_pattern = magSEs(f_ind,:,t_ind).';
        
                   [r,p] = corrcoef(seq_eff_pattern,slrp);
                   r_slrp(f_ind,t_ind,fly) = r(2); p_slrp(f_ind,t_ind,fly) = p(2); %#ok<PFOUS>
        
                   [r,p] = corrcoef(seq_eff_pattern,-lrpr);
                   r_lrpr(f_ind,t_ind,fly) = r(2); p_lrpr(f_ind,t_ind,fly) = p(2); %#ok<PFOUS>
        
                   [r,p] = corrcoef(seq_eff_pattern,weird);
                   r_weird(f_ind,t_ind,fly) = r(2); p_weird(f_ind,t_ind,fly) = p(2); %#ok<PFOUS>
                   
                  [r,p] = corrcoef(seq_eff_pattern,td_profile);
                   r_td(f_ind,t_ind,fly) = r(2); p_td(f_ind,t_ind,fly) = p(2); %#ok<PFOUS>
        
                end
                
            end
        end
        
        if plot_individual_flies
        
            if plotComponents && n_back == 5
                figure; 
                imagesc(r_slrp(:,:,fly),'xdata',time_bounds); colorbar; colormap('jet');
                set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
                xlabel('time (ms)'); ylabel('Frequency (Hz)');
                title(['SLRP Fly ' num2str(thisFly) figAdd]);
                %saveas(gcf,[resultsDirectory 'slrp_fly_' num2str(fly) '.png']);
                %saveas(gcf,[resultsDirectory 'slrp_fly_' num2str(thisFly) '.png']);
                saveas(gcf,[resultsDirectory 'slrp_fly_' num2str(thisFly) figAdd '.png']);
        
                figure; 
                imagesc(r_lrpr(:,:,fly),'xdata',time_bounds); colorbar; colormap('jet');
                set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
                xlabel('time (ms)'); ylabel('Frequency (Hz)');
                title(['LRPR Fly ' num2str(thisFly) figAdd]);
                %saveas(gcf,[resultsDirectory 'lrpr_fly_' num2str(fly) '.png']);
                %saveas(gcf,[resultsDirectory 'lrpr_fly_' num2str(thisFly) '.png']);
                saveas(gcf,[resultsDirectory 'lrpr_fly_' num2str(thisFly) figAdd '.png']);
        
                figure; 
                imagesc(r_weird(:,:,fly),'xdata',time_bounds); colorbar; colormap('jet');
                set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
                xlabel('time (ms)'); ylabel('Frequency (Hz)');
                title(['WEIRD Fly ' num2str(thisFly) figAdd]);
                %saveas(gcf,[resultsDirectory 'weird_fly_' num2str(fly) '.png']);
                saveas(gcf,[resultsDirectory 'weird_fly_' num2str(thisFly) figAdd '.png']);
                
                figure; 
                imagesc(r_td(:,:,fly),'xdata',time_bounds); colorbar; colormap('jet');
                set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
                xlabel('time (ms)'); ylabel('Frequency (Hz)');
                title(['TD PROFILE Fly ' num2str(thisFly) figAdd]);
                %saveas(gcf,[resultsDirectory 'td_profile_fly_' num2str(fly) '.png']);
                saveas(gcf,[resultsDirectory 'td_profile_fly_' num2str(thisFly) figAdd '.png']);
            end
    
            % spectrogram for AAAA minus AAAR
                %Replaced with custom comparisons
            %{
            figure; 
            imagesc(squeeze(magSEs(:,16,:)) - squeeze(magSEs(:,8,:)),'xdata',time_bounds); colorbar;
            set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
            xlabel('time (ms)'); ylabel('Frequency (Hz)');
            title('AAAA minus AAAR');
            saveas(gcf,[resultsDirectory 'AAAA_minus_AAAR_fly' num2str(fly) '.png']);
    
            % spectrogram for RRRR minus RRRA
            figure; imagesc(squeeze(magSEs(:,9,:)) - squeeze(magSEs(:,1,:)),'xdata',time_bounds); colorbar;
            set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
            xlabel('time (ms)'); ylabel('Frequency (Hz)');
            title('RRRR minus RRRA');
            saveas(gcf,[resultsDirectory 'RRRR_minus_RRRA_fly' num2str(fly) '.png']);
    
            % spectrogram for RRRA - AAAR
            figure; imagesc(squeeze(magSEs(:,9,:)) - squeeze(magSEs(:,8,:)),'xdata',time_bounds); colorbar;
            set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
            xlabel('time (ms)'); ylabel('Frequency (Hz)');
            title('RRRR minus RRRA');
            saveas(gcf,[resultsDirectory 'RRRR_minus_AAAR_fly' num2str(fly) '.png']);
            %}
    
            %Custom comparisons
            %{
            for comp = 1:size(customComparisons,1)
                figure; 
                c1 = customComparisons(comp,1);
                c2 = customComparisons(comp,2);
                imagesc(squeeze(magSEs(:,c1,:)) - squeeze(magSEs(:,c2,:)),'xdata',time_bounds); colorbar;
                set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
                xlabel('time (ms)'); ylabel('Frequency (Hz)');
                %title('AAAA minus AAAR');
                title(['Fly #', num2str(thisFly),' ',exLabels{c1},' minus ',exLabels{c2}]);
                %saveas(gcf,[resultsDirectory 'AAAA_minus_AAAR_fly' num2str(fly) '.png']);
                saveas(gcf,[resultsDirectory exLabels{c1} '_minus_' exLabels{c2} '_fly' num2str(thisFly) '.png']);
            end
            %}
            for comp = 1:size(customComparisons,1)
                figure; 
                c1 = customComparisons(comp,1);
                c2 = customComparisons(comp,2);
                if any([c1,c2] == -1) %Special case average requested
                    if c1 == -1
                        c1Actual = [1:size(magSEs,2)]; %All chans
                        c1ActualLabel = 'all seq av.';
                    else
                        c1Actual = c1;
                        c1ActualLabel = exLabels{c1Actual};
                    end
                    if c2 == -1
                        c2Actual = [1:size(magSEs,2)]; %All chans
                        c2ActualLabel = 'all seq av';
                    else
                        c2Actual = c2;
                        c2ActualLabel = exLabels{c2Actual};
                    end
                else %Normal case/s
                    c1Actual = c1;
                    c2Actual = c2;
                    c1ActualLabel = exLabels{c1Actual};
                    c2ActualLabel = exLabels{c2Actual};
                end
                %imagesc(squeeze(meanMagSEs(:,c1,:)) - squeeze(meanMagSEs(:,c2,:)),'xdata',time_bounds); colorbar;
                imagesc(squeeze( nanmean(magSEs(:,c1Actual,:),2) ) - squeeze( nanmean(magSEs(:,c2Actual,:),2) ),'xdata',time_bounds); colorbar; %Add in nanmean to allow for averaging
                set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
                xlabel('time (ms)'); ylabel('Frequency (Hz)');
                title(['Fly #', num2str(thisFly),' ',c1ActualLabel,' minus ',c2ActualLabel figAdd]);
                saveas(gcf,[resultsDirectory c1ActualLabel '_minus_' c2ActualLabel '_fly' num2str(thisFly) figAdd '.png']);
            end
    
        
        end
        
    end

    %Plots
    if plotComponents && n_back == 5
        figure; 
        imagesc(mean(r_slrp,3),'xdata',time_bounds); colorbar; colormap('jet');
        set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
        xlabel('time (ms)'); ylabel('Frequency (Hz)');
        title(['Mean SLRP' figAdd ]);
        saveas(gcf,[resultsDirectory 'mean_slrp' figAdd '.png']);
        
        figure;
        imagesc(mean(r_lrpr,3),'xdata',time_bounds); colorbar; colormap('jet');
        set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
        xlabel('time (ms)'); ylabel('Frequency (Hz)');
        title(['Mean LRPR' figAdd ]);
        saveas(gcf,[resultsDirectory 'mean_lrpr' figAdd '.png']);
        
        figure; 
        imagesc(mean(r_weird,3),'xdata',time_bounds); colorbar; colormap('jet');
        set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
        xlabel('time (ms)'); ylabel('Frequency (Hz)');
        title(['Mean WEIRD' figAdd]);
        saveas(gcf,[resultsDirectory 'mean_weird' figAdd '.png']);
        
        figure; 
        imagesc(mean(r_td,3),'xdata',time_bounds); colorbar; colormap('jet');
        set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
        xlabel('time (ms)'); ylabel('Frequency (Hz)');
        title(['Mean TD PROFILE' figAdd]);
        saveas(gcf,[resultsDirectory 'mean_td' figAdd '.png']);
    elseif plotComponents && n_back ~= 5
        disp(['Time frequency plotting of components requested but n_back ~= 5, and calculated components do not exist'])
    end

end %isoN end

%% plot mean oddballs

%magSEs = FLIES(1).magnitudeSEs; %Moved below because of isomer reasons
%magSEs = magSEs(f > f_limits(1) & f < f_limits(2),:,:);

for iso = 1:isoN
    %Quick preparation if isomer
    if isoMode
        thisIsom = isoList{iso}; %Will only be enacted in isomer mode, so no risk when isoN == 1 but isoMode == 0
        figAdd = ['-isomer',thisIsom];
        magSEs = FLIES(1).ISOMER.(thisIsom).magnitudeSEs; %Moved below because of isomer reasons
        magSEs = magSEs(f > f_limits(1) & f < f_limits(2),:,:);
    else
        figAdd = ['-merged'];
        magSEs = FLIES(1).magnitudeSEs; %Moved below because of isomer reasons
        magSEs = magSEs(f > f_limits(1) & f < f_limits(2),:,:);
    end

    %calculate mean spectra
    meanMagSEs = zeros(size(magSEs));
    
    for fly = 1:n_flies
    
        if ~isoMode
            magSEs = FLIES(fly).magnitudeSEs; magSEs = magSEs(f > f_limits(1) & f < f_limits(2),:,:);
        else
            magSEs = FLIES(fly).ISOMER.(thisIsom).magnitudeSEs; 
            magSEs = magSEs(f > f_limits(1) & f < f_limits(2),:,:);
        end

        meanMagSEs = meanMagSEs + magSEs;
    
    end
    
    meanMagSEs = meanMagSEs/n_flies;
    
    %{
    % spectrogram for AAAA minus AAAR
    figure; 
    imagesc(squeeze(meanMagSEs(:,16,:)) - squeeze(meanMagSEs(:,8,:)),'xdata',time_bounds); colorbar;
    set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
    xlabel('time (ms)'); ylabel('Frequency (Hz)');
    title('Mean AAAA minus AAAR');
    saveas(gcf,[resultsDirectory 'mean_AAAA_minus_AAAR.png']);
    
    % spectrogram for RRRR minus RRRA
    figure; imagesc(squeeze(meanMagSEs(:,9,:)) - squeeze(meanMagSEs(:,1,:)),'xdata',time_bounds); colorbar;
    set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
    xlabel('time (ms)'); ylabel('Frequency (Hz)');
    title('Mean RRRR minus RRRA');
    saveas(gcf,[resultsDirectory 'mean_RRRR_minus_RRRA.png']);
    %}
    
    %Custom comparisons
    for comp = 1:size(customComparisons,1)
        figure; 
        c1 = customComparisons(comp,1);
        c2 = customComparisons(comp,2);
        if any([c1,c2] == -1) %Special case average requested
            if c1 == -1
                c1Actual = [1:size(meanMagSEs,2)]; %All chans
                c1ActualLabel = 'all seq av.';
            else
                c1Actual = c1;
                c1ActualLabel = exLabels{c1Actual};
            end
            if c2 == -1
                c2Actual = [1:size(meanMagSEs,2)]; %All chans
                c2ActualLabel = 'all seq av';
            else
                c2Actual = c2;
                c2ActualLabel = exLabels{c2Actual};
            end
        else %Normal case/s
            c1Actual = c1;
            c2Actual = c2;
            c1ActualLabel = exLabels{c1Actual};
            c2ActualLabel = exLabels{c2Actual};
        end
        %imagesc(squeeze(meanMagSEs(:,c1,:)) - squeeze(meanMagSEs(:,c2,:)),'xdata',time_bounds); colorbar;
        imagesc(squeeze( nanmean(meanMagSEs(:,c1Actual,:),2) ) - squeeze( nanmean(meanMagSEs(:,c2Actual,:),2) ),'xdata',time_bounds); colorbar; %Add in nanmean to allow for averaging
        set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
        xlabel('time (ms)'); ylabel('Frequency (Hz)');
        %title('AAAA minus AAAR');
        %title(['Mean ', exLabels{c1},' minus ',exLabels{c2}]);
        title(['Mean ', c1ActualLabel,' minus ',c2ActualLabel figAdd]);
        %saveas(gcf,[resultsDirectory 'AAAA_minus_AAAR_fly' num2str(fly) '.png']);
        saveas(gcf,[resultsDirectory 'mean_' c1ActualLabel '_minus_' c2ActualLabel figAdd '.png']);
    end

end

%% plot mean r for all flies
    %Moved above to nullify issues with correct data selection
%{
if plotComponents && n_back == 5
    figure; 
    imagesc(mean(r_slrp,3),'xdata',time_bounds); colorbar; colormap('jet');
    set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
    xlabel('time (ms)'); ylabel('Frequency (Hz)');
    title('Mean SLRP');
    saveas(gcf,[resultsDirectory 'mean_slrp.png']);
    
    figure;
    imagesc(mean(r_lrpr,3),'xdata',time_bounds); colorbar; colormap('jet');
    set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
    xlabel('time (ms)'); ylabel('Frequency (Hz)');
    title('Mean LRPR');
    saveas(gcf,[resultsDirectory 'mean_lrpr.png']);
    
    figure; 
    imagesc(mean(r_weird,3),'xdata',time_bounds); colorbar; colormap('jet');
    set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
    xlabel('time (ms)'); ylabel('Frequency (Hz)');
    title('Mean WEIRD');
    saveas(gcf,[resultsDirectory 'mean_weird.png']);
    
    figure; 
    imagesc(mean(r_td,3),'xdata',time_bounds); colorbar; colormap('jet');
    set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
    xlabel('time (ms)'); ylabel('Frequency (Hz)');
    title('Mean TD PROFILE');
    saveas(gcf,[resultsDirectory 'mean_td.png']);
elseif plotComponents && n_back ~= 5
    disp(['Time frequency plotting of components requested but n_back ~= 5, and calculated components do not exist'])
end
%}
