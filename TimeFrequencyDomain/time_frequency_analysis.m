% analyse in which frequency bands SLRP and LRPR live

clear; close all;

load slrp_lrpr.mat

load data_flies_time_frequency_wavelet_12dot5.mat

FLIES = FLIES(~cellfun(@isempty,struct2cell(FLIES)));

homeDirectory = 'D:\group_swinderen\Dinis';

resultsDirectory = [homeDirectory '\Results_Time_Frequency\Wavelets\'];

n_flies = length(FLIES);

% frequency_bounds = [0 60];

time_bounds = [-10 70];

options = optimset('Algorithm','interior-point','FinDiffType','central');

resampleFreq = 3000;%CHECK

fly_index = 1;

f_limits = [0 300];

f = FLIES(1).LIT.f;

f_reduced = f(f > f_limits(1) & f < f_limits(2));

%for plots
y_data = [max(f_reduced) min(f_reduced)];

%for preallocation
fSize = length(f(f > f_limits(1) & f < f_limits(2)));
tSize = size(FLIES(1).LIT.magnitudeSEs,2);

%% calculate variance across sequences to identify regions of interest for
% sequential dependencies
% magnitudeVar = [];
% phaseVar = [];
% 
% for fly = 1:n_flies
%         
%         t = FLIES(fly).LIT.t*1000;
%         
%         magVar = squeeze(var(FLIES(fly).LIT.magnitudeSEs,0,2));
%         magVar = magVar(f > f_limits(1) & f < f_limits(2),:);
%         magnitudeVar(:,:,fly) = magVar;
% %         figure; imagesc(magVar,'xdata',[min(t) max(t)],'ydata',[min(f_reduced) max(f_reduced)]);
% %         title(['Variance magnitude fly' num2str(fly)]);
% 
%         phaVar = squeeze(var(FLIES(fly).LIT.phaseSEs,0,2));
%         phaVar = phaVar(f > f_limits(1) & f < f_limits(2),:);
%         phaseVar(:,:,fly) = phaVar; %#ok<*SAGROW>
% %         figure; imagesc(phaVar,'xdata',[min(t) max(t)],'ydata',[min(f_reduced) max(f_reduced)]);
% %         title(['Variance phase fly' num2str(fly)]);
%     
% end
% 
% figure; imagesc(mean(magnitudeVar, 3),'xdata',[min(t) max(t)],'ydata',[min(f_reduced) max(f_reduced)]); title('magnitude var');
% figure; imagesc(mean(phaseVar,3),'xdata',[min(t) max(t)],'ydata',[min(f_reduced) max(f_reduced)]); title('phase var');

%% calculate fit to SLRP and LRPR to individual flies (this will take a loooong time)

% r_squared_slrp = zeros(fSize,tSize,n_flies);
% r_squared_lrpr = zeros(fSize,tSize,n_flies);
% r_squared_overall = zeros(fSize,tSize,n_flies);
% r_squared_weird = zeros(fSize,tSize,n_flies);
% 
% for fly = 1:n_flies
%         
%         magSEs = FLIES(fly).LIT.magnitudeSEs; magSEs = magSEs(f > f_limits(1) & f < f_limits(2),:,:);
% 
%         for f_ind = 1:size(magSEs,1)% frequency index
% 
%             for t_ind = 1:size(magSEs,3)% time index
% 
%                 seq_eff_pattern = magSEs(f_ind,:,t_ind).';
% 
%                 %fit only to slrp or lrpr
%                [x_slrp,sse_slrp] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[-inf  0 0 -inf],[inf 0 0 inf],[],options);
%                [x_lrpr,sse_lrpr] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 1 1 0],[],[],[],[],[0  -inf 0 -inf],[0 inf 0 inf],[],options);
%                [x_weird,sse_weird] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);
% 
%                %fit overall
%                [x_overall,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 1 0 0],[],[],[],[],[-inf  -inf 0 -inf],[inf inf 0 inf],[],options);
% 
%                sse_total = sum((seq_eff_pattern-mean(seq_eff_pattern)).^2);
% 
%                r_squared_slrp(f_ind,t_ind,fly) = 1-(sse_slrp/sse_total);
%                r_squared_lrpr(f_ind,t_ind,fly) = 1-(sse_lrpr/sse_total);
%                r_squared_weird(f_ind,t_ind,fly) = 1-(sse_weird/sse_total);
%                r_squared_overall(f_ind,t_ind,fly) = 1-(sse_overall/sse_total);
% 
%             end
% 
%         end
%         
%         figure; title(['SLRP Fly ' num2str(fly)]);
%         imagesc(r_squared_slrp(:,:,fly),'xdata',time_bounds,'ydata',[min(f_reduced) max(f_reduced)]); colorbar;
%         xlabel('time (ms)'); ylabel('Frequency (Hz)');
%         saveas(gcf,[resultsDirectory 'slrp_fly_' num2str(fly) '.png']);
%         
%         figure; title(['LRPR Fly ' num2str(fly)]);
%         imagesc(r_squared_lrpr(:,:,fly),'xdata',time_bounds,'ydata',[min(f_reduced) max(f_reduced)]); colorbar;
%         xlabel('time (ms)'); ylabel('Frequency (Hz)');
%         saveas(gcf,[resultsDirectory 'lrpr_fly_' num2str(fly) '.png']);
%         
%         figure; title(['WEIRD Fly ' num2str(fly)]);
%         imagesc(r_squared_weird(:,:,fly),'xdata',time_bounds,'ydata',[min(f_reduced) max(f_reduced)]); colorbar;
%         xlabel('time (ms)'); ylabel('Frequency (Hz)');
%         saveas(gcf,[resultsDirectory 'weird_fly_' num2str(fly) '.png']);
%         
%         figure; title(['COMBINED Fly ' num2str(fly)]);
%         imagesc(r_squared_overall(:,:,fly),'xdata',time_bounds,'ydata',[min(f_reduced) max(f_reduced)]); colorbar;
%         xlabel('time (ms)'); ylabel('Frequency (Hz)');
%         saveas(gcf,[resultsDirectory 'combined_fly_' num2str(fly) '.png']);
%         
%         close all;
%     
% end
 
%% mean magnitude SEs across all flies
meanMagSEs = FLIES(1).LIT.magnitudeSEs;

for fly = 2:n_flies
   
    meanMagSEs = meanMagSEs + FLIES(fly).LIT.magnitudeSEs;
    
end

meanMagSEs = meanMagSEs/n_flies;

meanMagSEs = meanMagSEs(f > f_limits(1) & f < f_limits(2),:,:);

%% phase SEs

% meanMagSEs = FLIES(1).LIT.phaseSEs;
% 
% for fly = 2:n_flies
%    
%     meanMagSEs = meanMagSEs + FLIES(fly).LIT.phaseSEs;
%     
% end
% 
% meanMagSEs = meanMagSEs/n_flies;
% 
% meanMagSEs = meanMagSEs(f > f_limits(1) & f < f_limits(2),:,:);

%% calculate overall fit to SLRP and LRPR

r_squared_slrp = zeros(fSize,tSize);
r_squared_lrpr = zeros(fSize,tSize);
r_squared_overall = zeros(fSize,tSize);
r_squared_weird = zeros(fSize,tSize);

for f_ind = 1:size(meanMagSEs,1)% frequency index

    for t_ind = 1:size(meanMagSEs,3)% time index

        seq_eff_pattern = meanMagSEs(f_ind,:,t_ind).';

        %fit only to slrp or lrpr
       [x_slrp,sse_slrp] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[-inf  0 0 -inf],[inf 0 0 inf],[],options);
       [x_lrpr,sse_lrpr] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 1 1 0],[],[],[],[],[0  -inf 0 -inf],[0 inf 0 inf],[],options);
       [x_weird,sse_weird] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);

       %fit overall (no weird for 1.25 Hz)
       [x_overall,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 1 0 0],[],[],[],[],[-inf  -inf 0 -inf],[inf inf 0 inf],[],options);

       sse_total = sum((seq_eff_pattern-mean(seq_eff_pattern)).^2);

       r_squared_slrp(f_ind,t_ind) = 1-(sse_slrp/sse_total);
       r_squared_lrpr(f_ind,t_ind) = 1-(sse_lrpr/sse_total);
       r_squared_weird(f_ind,t_ind) = 1-(sse_weird/sse_total);
       r_squared_overall(f_ind,t_ind) = 1-(sse_overall/sse_total);

    end

end

%% plot r_squared

%make a tick every x steps
y_ticks = 1:2:length(f);
y_tick_labels = f(y_ticks);

figure; title('SLRP All');
imagesc(r_squared_slrp(:,:),'xdata',time_bounds); colorbar; colormap('hot');
set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
xlabel('time (ms)'); ylabel('Frequency (Hz)');
saveas(gcf,[resultsDirectory 'slrp_all.png']);

figure; title('LRPR All');
imagesc(r_squared_lrpr,'xdata',time_bounds); colorbar; colormap('hot');
set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
xlabel('time (ms)'); ylabel('Frequency (Hz)');
saveas(gcf,[resultsDirectory 'lrpr_all.png']);

figure; title('WEIRD All');
imagesc(r_squared_weird,'xdata',time_bounds); colorbar; colormap('hot');
set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
xlabel('time (ms)'); ylabel('Frequency (Hz)');
saveas(gcf,[resultsDirectory 'weird_all.png']);

figure; title('COMBINED All ');
imagesc(r_squared_overall,'xdata',time_bounds); colorbar; colormap('hot');
set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
xlabel('time (ms)'); ylabel('Frequency (Hz)');
saveas(gcf,[resultsDirectory 'combined_all.png']);

%% create plots for difference between combined and slrp/lrpr

figure; title(['Overall minus SLRP Fly ' num2str(fly)]);
imagesc(r_squared_overall-r_squared_slrp,'xdata',time_bounds); colorbar; colormap('hot');
set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
xlabel('time (ms)'); ylabel('Frequency (Hz)');
saveas(gcf,[resultsDirectory 'combined_minus_slrp.png']);

figure; title(['Overall minus LRPR Fly ' num2str(fly)]);
imagesc(r_squared_overall-r_squared_lrpr,'xdata',time_bounds); colorbar; colormap('hot');
set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
xlabel('time (ms)'); ylabel('Frequency (Hz)');
saveas(gcf,[resultsDirectory 'combined_minus_lrpr.png']);

%% plot spectrograms

% for i = 1:16 
% 
%     figure;
%     imagesc(squeeze(meanMagSEs(:,i,:)),'xdata',time_bounds); colorbar;
%     set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
%     xlabel('time (ms)'); ylabel('Frequency (Hz)');
%     saveas(gcf,[resultsDirectory 'spectrogram_' num2str(i) '.png']);
%     
% end

%% fit and plot some seq eff profiles of interest

% 1.25 Hz cases
% profiles = [32 410; 37 1776; 15 365; 41 611; 51 711; 11 428; 24 558; 11 34; 32 837];
% 
% type = [1 1 2 2 2 3 3 4 4];

% 6.25 Hz case
% profiles = [34 223;38 303;13 21;19 377];
% 
% type = [1 2 3 4];

% 12.5 Hz
profiles = [22 175;28 238;28 98;12 210];

type = [1 2 3 4];

% 25 Hz
% profiles = [16 49;6 60; 10 85; 7 92; 9 80];
% 
% type = [1 2 2 3 4];

%25 Hz phase
% profiles = [8 69];
% 
% type = [2];

names = {'slrp','lrpr','weird','overall'};

for i = 1:size(profiles,1)
    
    seq_eff_pattern = meanMagSEs(profiles(i,1),:,profiles(i,2)).';
   
    switch type(i)
        
        case 1
  
           %fit only to slrp or lrpr
           [x,sse_slrp] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[-inf  0 0 -inf],[inf 0 0 inf],[],options);
         
        case 2
            
           [x,sse_lrpr] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 1 1 0],[],[],[],[],[0  -inf 0 -inf],[0 inf 0 inf],[],options);
            
        case 3
            
           [x,sse_weird] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);
            
        case 4
            
           %fit overall (no weird for 1.25 Hz)
            [x,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 1 0 0],[],[],[],[],[-inf  -inf 0 -inf],[inf inf 0 inf],[],options);

    end
    
    seq_eff_fit = x(1)*slrp + x(2)*lrpr + x(3)*weird + x(4);
    
    figure; create_seq_eff_plot(seq_eff_pattern,seq_eff_fit);
    ylabel('Magnitude');
    saveas(gcf,[ resultsDirectory '/' names{type(i)} '_' num2str(profiles(i,1)) '_' num2str(profiles(i,2)) '.png'])
    
end

%% plot RRRR-RRRA and AAAA-AAAR

figure; imagesc(squeeze(meanMagSEs(:,16,:))-squeeze(meanMagSEs(:,8,:)),'xdata',time_bounds); colorbar;
set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
xlabel('time (ms)'); ylabel('Frequency (Hz)');
saveas(gcf,[ resultsDirectory '/' 'AAAA_minus_AAAR.png']);

figure; imagesc(squeeze(meanMagSEs(:,9,:))-squeeze(meanMagSEs(:,1,:)),'xdata',time_bounds); colorbar;
set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
xlabel('time (ms)'); ylabel('Frequency (Hz)');
saveas(gcf,[ resultsDirectory '/' 'RRRR_minus_RRRA.png']);

%% plot overall ERP

overallERP = zeros(length(FLIES(1).LIT.meanERPs),n_flies);

for i = 1:n_flies
    
    overallERP(:,i) = mean(FLIES(i).LIT.meanERPs,2);
    
    figure; plot(smoothdata(overallERP(:,i),'sgolay'),'r','linewidth',2); ylabel('Amplitude (ms)');
    set(gca,'xtick',.1*3000,'xticklabel',0);
    saveas(gcf,[resultsDirectory 'ERP_fly ' num2str(i) '.png']);
    
end

plot(smoothdata(mean(overallERP,2),'sgolay'),'r','linewidth',2); ylabel('Amplitude (ms)');
set(gca,'xtick',.1*3000,'xticklabel',0);
saveas(gcf,[resultsDirectory 'ERP_all_flies.png']);