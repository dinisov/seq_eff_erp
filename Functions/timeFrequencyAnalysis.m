% analyse in which frequency bands SLRP and LRPR live
function timeFrequencyAnalysis(FLIES, homeDirectory, plot_individual_flies)

load('slrp_lrpr.mat','slrp','lrpr','weird');

% gets rid of empty structures
FLIES = FLIES(~cellfun(@isempty,{FLIES.allERPs}));

resultsDirectory = [homeDirectory '\Results_Time_Frequency\'];

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

%% calculate fit to SLRP and LRPR to individual flies (r^2 plots)

% r_squared_slrp = zeros(fSize,tSize,n_flies);
% r_squared_lrpr = zeros(fSize,tSize,n_flies);
% r_squared_overall = zeros(fSize,tSize,n_flies);
% r_squared_weird = zeros(fSize,tSize,n_flies);
% 
% for fly = 1:n_flies
% 
%         magSEs = FLIES(fly).magnitudeSEs; magSEs = magSEs(f > f_limits(1) & f < f_limits(2),:,:);
% 
%         f_size = size(magSEs,1); t_size = size(magSEs,3);
%         
%         parfor f_ind = 1:f_size% frequency index
% 
%             for t_ind = 1:t_size% time index
% 
%                 seq_eff_pattern = magSEs(f_ind,:,t_ind).';
% 
%                 %fit only to slrp or lrpr
%                [~,sse_slrp] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[-inf  0 0 -inf],[inf 0 0 inf],[],options);
%                [~,sse_lrpr] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 1 1 0],[],[],[],[],[0  -inf 0 -inf],[0 inf 0 inf],[],options);
%                [~,sse_weird] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);
% 
%                %fit overall
%                [~,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 1 0 0],[],[],[],[],[-inf  -inf 0 -inf],[inf inf 0 inf],[],options);
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
%         imagesc(r_squared_slrp,'xdata',time_bounds); colorbar; colormap('hot');
%         set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
%         xlabel('time (ms)'); ylabel('Frequency (Hz)');
%         saveas(gcf,[resultsDirectory 'slrp_fly_' num2str(fly) '.png']);
% 
%         figure; title(['LRPR Fly ' num2str(fly)]);
%         imagesc(r_squared_lrpr,'xdata',time_bounds); colorbar; colormap('hot');
%         set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
%         xlabel('time (ms)'); ylabel('Frequency (Hz)');
%         saveas(gcf,[resultsDirectory 'lrpr_fly_' num2str(fly) '.png']);
% 
%         figure; title(['WEIRD Fly ' num2str(fly)]);
%         imagesc(r_squared_weird,'xdata',time_bounds); colorbar; colormap('hot');
%         set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
%         xlabel('time (ms)'); ylabel('Frequency (Hz)');
%         saveas(gcf,[resultsDirectory 'weird_fly_' num2str(fly) '.png']);
% 
%         figure; title(['COMBINED Fly ' num2str(fly)]);
%         imagesc(r_squared_overall,'xdata',time_bounds); colorbar; colormap('hot');
%         set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
%         xlabel('time (ms)'); ylabel('Frequency (Hz)');
%         saveas(gcf,[resultsDirectory 'combined_fly_' num2str(fly) '.png']);
% 
% %         close all;
% 
% end

%% calculate r and p values (r plots)

r_slrp = zeros(fSize,tSize,n_flies);
r_lrpr = zeros(fSize,tSize,n_flies);
r_weird = zeros(fSize,tSize,n_flies);
r_td = zeros(fSize,tSize,n_flies);

p_slrp = zeros(fSize,tSize,n_flies);
p_lrpr = zeros(fSize,tSize,n_flies);
p_weird = zeros(fSize,tSize,n_flies);
p_td = zeros(fSize,tSize,n_flies);

for fly = 1:n_flies
    
    td_profile = FLIES(fly).PROFILE.amplitude;
   
    magSEs = FLIES(fly).magnitudeSEs; magSEs = magSEs(f > f_limits(1) & f < f_limits(2),:,:);

    f_size = size(magSEs,1); t_size = size(magSEs,3);

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
    
    if plot_individual_flies
    
        figure; title(['SLRP Fly ' num2str(fly)]);
        imagesc(r_slrp(:,:,fly),'xdata',time_bounds); colorbar; colormap('jet');
        set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
        xlabel('time (ms)'); ylabel('Frequency (Hz)');
        saveas(gcf,[resultsDirectory 'slrp_fly_' num2str(fly) '.png']);

        figure; title(['LRPR Fly ' num2str(fly)]);
        imagesc(r_lrpr(:,:,fly),'xdata',time_bounds); colorbar; colormap('jet');
        set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
        xlabel('time (ms)'); ylabel('Frequency (Hz)');
        saveas(gcf,[resultsDirectory 'lrpr_fly_' num2str(fly) '.png']);

        figure; title(['WEIRD Fly ' num2str(fly)]);
        imagesc(r_weird(:,:,fly),'xdata',time_bounds); colorbar; colormap('jet');
        set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
        xlabel('time (ms)'); ylabel('Frequency (Hz)');
        saveas(gcf,[resultsDirectory 'weird_fly_' num2str(fly) '.png']);

        figure; title(['TD PROFILE Fly ' num2str(fly)]);
        imagesc(r_td(:,:,fly),'xdata',time_bounds); colorbar; colormap('jet');
        set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
        xlabel('time (ms)'); ylabel('Frequency (Hz)');
        saveas(gcf,[resultsDirectory 'td_profile_fly_' num2str(fly) '.png']);
    
    end
    
end
 
%% plot mean r for all flies

figure; title('Mean SLRP');
imagesc(mean(r_slrp,3),'xdata',time_bounds); colorbar; colormap('jet');
set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
xlabel('time (ms)'); ylabel('Frequency (Hz)');
saveas(gcf,[resultsDirectory 'mean_slrp.png']);

figure; title('Mean LRPR');
imagesc(mean(r_lrpr,3),'xdata',time_bounds); colorbar; colormap('jet');
set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
xlabel('time (ms)'); ylabel('Frequency (Hz)');
saveas(gcf,[resultsDirectory 'mean_lrpr.png']);

figure; title('Mean WEIRD');
imagesc(mean(r_weird,3),'xdata',time_bounds); colorbar; colormap('jet');
set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
xlabel('time (ms)'); ylabel('Frequency (Hz)');
saveas(gcf,[resultsDirectory 'mean_weird.png']);

figure; title('Mean TD PROFILE');
imagesc(mean(r_td,3),'xdata',time_bounds); colorbar; colormap('jet');
set(gca,'ytick',y_ticks,'yticklabel',round(y_tick_labels,1));
xlabel('time (ms)'); ylabel('Frequency (Hz)');
saveas(gcf,[resultsDirectory 'mean_td.png']);

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

% r_squared_slrp = zeros(fSize,tSize);
% r_squared_lrpr = zeros(fSize,tSize);
% r_squared_overall = zeros(fSize,tSize);
% r_squared_weird = zeros(fSize,tSize);
% 
% for f_ind = 1:size(meanMagSEs,1)% frequency index
% 
%     for t_ind = 1:size(meanMagSEs,3)% time index
% 
%         seq_eff_pattern = meanMagSEs(f_ind,:,t_ind).';
% 
%         %fit only to slrp or lrpr
%        [x_slrp,sse_slrp] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[-inf  0 0 -inf],[inf 0 0 inf],[],options);
%        [x_lrpr,sse_lrpr] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 1 1 0],[],[],[],[],[0  -inf 0 -inf],[0 inf 0 inf],[],options);
%        [x_weird,sse_weird] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);
% 
%        %fit overall (no weird for 1.25 Hz)
%        [x_overall,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 1 0 0],[],[],[],[],[-inf  -inf 0 -inf],[inf inf 0 inf],[],options);
% 
%        sse_total = sum((seq_eff_pattern-mean(seq_eff_pattern)).^2);
% 
%        r_squared_slrp(f_ind,t_ind) = 1-(sse_slrp/sse_total);
%        r_squared_lrpr(f_ind,t_ind) = 1-(sse_lrpr/sse_total);
%        r_squared_weird(f_ind,t_ind) = 1-(sse_weird/sse_total);
%        r_squared_overall(f_ind,t_ind) = 1-(sse_overall/sse_total);
% 
%     end
% 
% end

%% plot r_squared

% %make a tick every x steps
% y_ticks = 1:2:length(f);
% y_tick_labels = f(y_ticks);
% 
% figure; title('SLRP All');
% imagesc(r_squared_slrp(:,:),'xdata',time_bounds); colorbar; colormap('hot');
% set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
% xlabel('time (ms)'); ylabel('Frequency (Hz)');
% saveas(gcf,[resultsDirectory 'slrp_all.png']);
% 
% figure; title('LRPR All');
% imagesc(r_squared_lrpr,'xdata',time_bounds); colorbar; colormap('hot');
% set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
% xlabel('time (ms)'); ylabel('Frequency (Hz)');
% saveas(gcf,[resultsDirectory 'lrpr_all.png']);
% 
% figure; title('WEIRD All');
% imagesc(r_squared_weird,'xdata',time_bounds); colorbar; colormap('hot');
% set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
% xlabel('time (ms)'); ylabel('Frequency (Hz)');
% saveas(gcf,[resultsDirectory 'weird_all.png']);
% 
% figure; title('COMBINED All ');
% imagesc(r_squared_overall,'xdata',time_bounds); colorbar; colormap('hot');
% set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
% xlabel('time (ms)'); ylabel('Frequency (Hz)');
% saveas(gcf,[resultsDirectory 'combined_all.png']);

%% create plots for difference between combined and slrp/lrpr

% figure; title(['Overall minus SLRP Fly ' num2str(fly)]);
% imagesc(r_squared_overall-r_squared_slrp,'xdata',time_bounds); colorbar; colormap('hot');
% set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
% xlabel('time (ms)'); ylabel('Frequency (Hz)');
% saveas(gcf,[resultsDirectory 'combined_minus_slrp.png']);
% 
% figure; title(['Overall minus LRPR Fly ' num2str(fly)]);
% imagesc(r_squared_overall-r_squared_lrpr,'xdata',time_bounds); colorbar; colormap('hot');
% set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
% xlabel('time (ms)'); ylabel('Frequency (Hz)');
% saveas(gcf,[resultsDirectory 'combined_minus_lrpr.png']);

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

% % 1.25 Hz cases
% % profiles = [32 411;37 1774;41 610;50 756;24 558];
% % 
% % type = [1 1 2 2 3];
% 
% % 6.25 Hz case
% profiles = [34 223;38 303;13 21;19 377];
% 
% type = [1 2 3 4];
% 
% % 12.5 Hz
% % profiles = [22 175;28 238;28 98;12 210];
% % 
% % type = [1 2 3 4];
% 
% % 25 Hz
% % profiles = [16 49;6 60; 10 85; 7 92; 9 80];
% % 
% % type = [1 2 2 3 4];
% 
% %25 Hz phase
% % profiles = [8 69];
% % 
% % type = [2];
% 
% names = {'slrp','lrpr','weird','overall'};
% 
% for i = 1:size(profiles,1)
%     
%     seq_eff_pattern = meanMagSEs(profiles(i,1),:,profiles(i,2)).';
%    
%     switch type(i)
%         
%         case 1
%   
%            %fit only to slrp or lrpr
%            [x,sse_slrp] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[-inf  0 0 -inf],[inf 0 0 inf],[],options);
%          
%         case 2
%             
%            [x,sse_lrpr] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[0 1 1 0],[],[],[],[],[0  -inf 0 -inf],[0 inf 0 inf],[],options);
%             
%         case 3
%             
%            [x,sse_weird] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 0 1 0],[],[],[],[],[0  0 -inf -inf],[0 0 inf inf],[],options);
%             
%         case 4
%             
%            %fit overall (no weird for 1.25 Hz)
%             [x,sse_overall] = fmincon(@(x) least_squares_slrp_lrpr_weird(x(1),x(2),x(3),x(4),slrp,lrpr,weird,seq_eff_pattern),[1 1 0 0],[],[],[],[],[-inf  -inf 0 -inf],[inf inf 0 inf],[],options);
% 
%     end
%     
%     seq_eff_fit = x(1)*slrp + x(2)*lrpr + x(3)*weird + x(4);
%     
%     figure; create_seq_eff_plot(seq_eff_pattern,seq_eff_fit);
%     ylabel('Magnitude');
%     saveas(gcf,[ resultsDirectory '/' names{type(i)} '_' num2str(profiles(i,1)) '_' num2str(profiles(i,2)) '.png'])
%     
% end

%% plot RRRR-RRRA and AAAA-AAAR

% figure; imagesc(squeeze(meanMagSEs(:,16,:))-squeeze(meanMagSEs(:,8,:)),'xdata',time_bounds); colorbar;
% set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
% xlabel('time (ms)'); ylabel('Frequency (Hz)');
% saveas(gcf,[ resultsDirectory '/' 'AAAA_minus_AAAR.png']);
% 
% figure; imagesc(squeeze(meanMagSEs(:,9,:))-squeeze(meanMagSEs(:,1,:)),'xdata',time_bounds); colorbar;
% set(gca,'ytick',y_ticks,'yticklabel',floor(y_tick_labels));
% xlabel('time (ms)'); ylabel('Frequency (Hz)');
% saveas(gcf,[ resultsDirectory '/' 'RRRR_minus_RRRA.png']);

%% plot overall ERP

% overallERP = zeros(length(FLIES(1).LIT.meanERPs),n_flies);
% 
% for i = 1:n_flies
%     
%     overallERP(:,i) = mean(FLIES(i).LIT.meanERPs,2);
%     
%     figure; plot(smoothdata(overallERP(:,i),'sgolay'),'r','linewidth',2); ylabel('Amplitude (ms)');
%     set(gca,'xtick',.1*3000,'xticklabel',0);
%     saveas(gcf,[resultsDirectory 'ERP_fly ' num2str(i) '.png']);
%     
% end
% 
% plot(smoothdata(mean(overallERP,2),'sgolay'),'r','linewidth',2); ylabel('Amplitude (ms)');
% set(gca,'xtick',.1*3000,'xticklabel',0);
% saveas(gcf,[resultsDirectory 'ERP_all_flies.png']);
