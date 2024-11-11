% Plotting lagged correlation between MJO and EP SSTA or CP SSTA
% conditioned on CP events, EP events, and El Nino events. We begin with
% the model simulation and then do the same for observational data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Conditioned lagged correlation from model %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%counting ENSO events occuring in each month since MJO is on a
%shorter time scale
window = 30;
total_loop=floor(length(T_E_3R_daily)/window)-1;

EN_monthly=zeros(total_loop,1);
EPEN_monthly=zeros(total_loop,1);
CPEN_monthly=zeros(total_loop,1);
LN_monthly=zeros(total_loop,1);

for k = 1:total_loop
    if mean(T_E_3R_daily(-0+k*window:+2+k*window)) > 0.5 || mean(T_C_3R_daily(-0+k*window:+2+k*window)) > 0.5
        EN_monthly(k)=1;
        if mean(T_E_3R_daily(-0+k*window:+2+k*window))>mean(T_C_3R_daily(-0+k*window:+2+k*window))
            EPEN_monthly(k)=1;
        else
            CPEN_monthly(k)=1;
        end
    elseif mean(T_E_3R_daily(-0+k*window:+2+k*window)) < -0.5 || mean(T_C_3R_daily(-0+k*window:+2+k*window)) < -0.5
        LN_monthly(k)=1;
    end
end



% subsampling of the intraseasonal wind to get daily data
intra_wind2 = movmean(intra_wind,k_dt_daily);  intra_wind2 = intra_wind2(1:k_dt_daily:end);
% variance of the intraseasonal wind
intra_wind3 = abs(intra_wind2);
%subsampling tau to get daily data
all_wind = movmean(tau,k_dt_daily); all_wind = all_wind(1:k_dt_daily:end);
%variance of tau
all_wind2 = abs(all_wind);

% MJO and MJO variances in the Indian Ocean
MJO_index = mean(MJO(1:3,1:k_dt_daily:end));
MJO_index2 = abs(MJO_index);
% MJO and MJO variances in the Pacific Ocean
MJO_index3 = mean(MJO(3:5,1:k_dt_daily:end));
MJO_index4 = abs(MJO_index3);


% Extract MJO segments for each ENSO type
window_months = 6;

lead_lag = 360;
window_lag_lead = lead_lag*2+1; % lead and lag are both 24 months

L2 = 13920;%number of days in each set that we find the correlation coefficent for so that we can calculate a confidence interval
times = floor(length(intra_wind2)/L2); % number of sets of L2 days in in the simulation
corr_coef1 = zeros(times,window_lag_lead); % MJO in Pacific ocean
corr_coef2 = zeros(times,window_lag_lead); % variance of MJO in Pacific ocean
corr_coef3 = zeros(times,window_lag_lead); % MJO in Pacific ocean
corr_coef4 = zeros(times,window_lag_lead); % variance of MJO in Pacific ocean
corr_coef5 = zeros(times,window_lag_lead); % MJO in Pacific ocean
corr_coef6 = zeros(times,window_lag_lead); % variance of MJO in Pacific ocean
corr_coef7 = zeros(times,window_lag_lead); % MJO in Pacific ocean
corr_coef8 = zeros(times,window_lag_lead); % variance of MJO in Pacific ocean


for i=1:times
    %extracting the time segments of the mean MJO and variance of MJO in the pacific
    % and EP SST when EP events are occuring
    EP_mjo1 = extract_mjo_segments(MJO_index3(1+(i-1)*L2:L2+(i-1)*L2), EPEN_monthly(1+(i-1)*L2/30:i*L2/30), window_months);
    EP_mjo2 = extract_mjo_segments(MJO_index4(1+(i-1)*L2:L2+(i-1)*L2), EPEN_monthly(1+(i-1)*L2/30:i*L2/30), window_months);
    EP_sst1 = extract_mjo_segments(T_E_3R_daily(1+(i-1)*L2:L2+(i-1)*L2)', EPEN_monthly(1+(i-1)*L2/30:i*L2/30), window_months);
    
    %extracting the time segments of the mean MJO and variance of MJO in the pacific
    % and EP SST when El Nino events are occuring
    EP_mjo3 = extract_mjo_segments(MJO_index3(1+(i-1)*L2:L2+(i-1)*L2), EN_monthly(1+(i-1)*L2/30:i*L2/30), window_months);
    EP_mjo4 = extract_mjo_segments(MJO_index4(1+(i-1)*L2:L2+(i-1)*L2), EN_monthly(1+(i-1)*L2/30:i*L2/30), window_months);
    EP_sst2 = extract_mjo_segments(T_E_3R_daily(1+(i-1)*L2:L2+(i-1)*L2)', EN_monthly(1+(i-1)*L2/30:i*L2/30), window_months);
    
    
    
    % extracting the time segments of the mean MJO and variance of MJO in the pacific
    % and CP SST when CP events are occuring
    CP_mjo1 = extract_mjo_segments(MJO_index3(1+(i-1)*L2:L2+(i-1)*L2), CPEN_monthly(1+(i-1)*L2/30:i*L2/30), window_months);
    CP_mjo2 = extract_mjo_segments(MJO_index4(1+(i-1)*L2:L2+(i-1)*L2), CPEN_monthly(1+(i-1)*L2/30:i*L2/30), window_months);
    CP_sst1 = extract_mjo_segments(T_C_3R_daily(1+(i-1)*L2:L2+(i-1)*L2)', CPEN_monthly(1+(i-1)*L2/30:i*L2/30), window_months);
    
    %extracting the time segments of the mean MJO and variance of MJO in the pacific
    % and CP SST when El Nino events are occuring
    CP_mjo3 = extract_mjo_segments(MJO_index3(1+(i-1)*L2:L2+(i-1)*L2), EN_monthly(1+(i-1)*L2/30:i*L2/30), window_months);
    CP_mjo4 = extract_mjo_segments(MJO_index4(1+(i-1)*L2:L2+(i-1)*L2), EN_monthly(1+(i-1)*L2/30:i*L2/30), window_months);
    CP_sst2 = extract_mjo_segments(T_C_3R_daily(1+(i-1)*L2:L2+(i-1)*L2)', EN_monthly(1+(i-1)*L2/30:i*L2/30), window_months);


    k = 1;% to index for different lags/leads
    % computing the lagged correlation at different lag/leads j
    for j = -lead_lag:lead_lag
        %correlation between WP mean MJO and TE during EP events
        temp=corrcoef(EP_mjo1(window_lag_lead+j:end-window_lag_lead+j), EP_sst1(window_lag_lead:end-window_lag_lead));
        corr_coef1(i,k) = temp(1,2);
        %correlation between WP MJO variance and TE during EP events
        temp=corrcoef(EP_mjo2(window_lag_lead+j:end-window_lag_lead+j), EP_sst1(window_lag_lead:end-window_lag_lead));
        corr_coef2(i,k) = temp(1,2);
        %correlation between WP mean MJO and TE during El Nino events
        temp=corrcoef(EP_mjo3(window_lag_lead+j:end-window_lag_lead+j), EP_sst2(window_lag_lead:end-window_lag_lead));
        corr_coef3(i,k) = temp(1,2);
        %correlation between WP MJO variance and TE during El Nino events
        temp=corrcoef(EP_mjo4(window_lag_lead+j:end-window_lag_lead+j), EP_sst2(window_lag_lead:end-window_lag_lead));
        corr_coef4(i,k) = temp(1,2);
        %correlation between WP mean MJO and TC during CP events
        temp=corrcoef(CP_mjo1(window_lag_lead+j:end-window_lag_lead+j), CP_sst1(window_lag_lead:end-window_lag_lead));
        corr_coef5(i,k) = temp(1,2);
        %correlation between WP MJO variance and TC during CP events
        temp=corrcoef(CP_mjo2(window_lag_lead+j:end-window_lag_lead+j), CP_sst1(window_lag_lead:end-window_lag_lead));
        corr_coef6(i,k) = temp(1,2);
        %correlation between WP mean MJO and TC during El Nino events
        temp=corrcoef(CP_mjo3(window_lag_lead+j:end-window_lag_lead+j), CP_sst2(window_lag_lead:end-window_lag_lead));
        corr_coef7(i,k) = temp(1,2);
        %correlation between WP MJO variance and TC during El Nino events
        temp=corrcoef(CP_mjo4(window_lag_lead+j:end-window_lag_lead+j), CP_sst2(window_lag_lead:end-window_lag_lead));
        corr_coef8(i,k) = temp(1,2);
        k = k + 1;
    end
end

%take the mean and standard devation of each correlation coefficient
%accross the sets of L2 days
corr_coef1_m=nanmean(corr_coef1,1);
corr_coef1_std=nanstd(corr_coef1,0,1);
corr_coef2_m=nanmean(corr_coef2,1);
corr_coef2_std=nanstd(corr_coef2,0,1);
corr_coef3_m=nanmean(corr_coef3,1);
corr_coef3_std=nanstd(corr_coef3,0,1);
corr_coef4_m=nanmean(corr_coef4,1);
corr_coef4_std=nanstd(corr_coef4,0,1);
corr_coef5_m=nanmean(corr_coef5,1);
corr_coef5_std=nanstd(corr_coef5,0,1);
corr_coef6_m=nanmean(corr_coef6,1);
corr_coef6_std=nanstd(corr_coef6,0,1);
corr_coef7_m=nanmean(corr_coef7,1);
corr_coef7_std=nanstd(corr_coef7,0,1);
corr_coef8_m=nanmean(corr_coef8,1);
corr_coef8_std=nanstd(corr_coef8,0,1);


xx = (-lead_lag:lead_lag)/30;%x axis in months

figure
subplot(2,4,6)
facealpha = 0.2;
hold on

%plot correlation between WP MJO variance and TE during EP events
hh2 = fill([xx fliplr(xx)],...
    [corr_coef2_m-corr_coef2_std fliplr(corr_coef2_m+corr_coef2_std)],'b','facealpha',facealpha);
hh2.EdgeColor = 'none';
plot((-lead_lag:lead_lag)/30,corr_coef2_m,'b','linewidth',2)

%plot correlation between WP mean MJO and TE during EP events
hh1 = fill([xx fliplr(xx)],...
    [corr_coef1_m-corr_coef1_std fliplr(corr_coef1_m+corr_coef1_std)],'r','facealpha',facealpha);
hh1.EdgeColor = 'none';
plot((-lead_lag:lead_lag)/30,corr_coef1_m,'r','linewidth',2)
box on
set(gca,'fontsize',14)
xlabel('(MJOI leads)  <-- lag -->  (T_E leads)')
ylabel('Corr')
title('(f) EP El Nino events')
ylim([-0.2,0.4])
xlim([-12,12])


subplot(2,4,2)
facealpha = 0.2;
hold on

%plot correlation between WP MJO variance and TE during El Nino events
hh2 = fill([xx fliplr(xx)],...
    [corr_coef4_m-corr_coef4_std fliplr(corr_coef4_m+corr_coef4_std)],'b','facealpha',facealpha);
hh2.EdgeColor = 'none';
plot((-lead_lag:lead_lag)/30,corr_coef4_m,'b','linewidth',2)

%plot correlation between WP MJO variance and TE during El Nino events
hh1 = fill([xx fliplr(xx)],...
    [corr_coef3_m-corr_coef3_std fliplr(corr_coef3_m+corr_coef3_std)],'r','facealpha',facealpha);
hh1.EdgeColor = 'none';
plot((-lead_lag:lead_lag)/30,corr_coef3_m,'r','linewidth',2)
box on
set(gca,'fontsize',14)
xlabel('(MJOI leads)  <-- lag -->  (T_E leads)')
ylabel('Corr')
title({'Model results';'Lagged correlation with T_E';'(a) All El Nino events'})
ylim([-0.2,0.4])
xlim([-12,12])


subplot(2,4,8)
facealpha = 0.2;
hold on 

% plot correlation between WP MJO variance and TC during CP events
hh2 = fill([xx fliplr(xx)],...
    [corr_coef6_m-corr_coef6_std fliplr(corr_coef6_m+corr_coef6_std)],'b','facealpha',facealpha);
hh2.EdgeColor = 'none';
plot((-lead_lag:lead_lag)/30,corr_coef6_m,'b','linewidth',2)

%plot correlation between WP MJO variance and TC during CP events
hh1 = fill([xx fliplr(xx)],...
    [corr_coef5_m-corr_coef1_std fliplr(corr_coef5_m+corr_coef5_std)],'r','facealpha',facealpha);
hh1.EdgeColor = 'none';
plot((-lead_lag:lead_lag)/30,corr_coef5_m,'r','linewidth',2)
box on
set(gca,'fontsize',14)
xlabel('(MJOI leads)  <-- lag -->  (T_C leads)')
ylabel('Corr')
title('(h) CP El Nino events')
ylim([-0.2,0.4])
xlim([-12,12])


subplot(2,4,4)
facealpha = 0.2;
hold on

%plot correlation between WP MJO variance and TC during El Nino events
hh2 = fill([xx fliplr(xx)],...
    [corr_coef8_m-corr_coef8_std fliplr(corr_coef8_m+corr_coef8_std)],'b','facealpha',facealpha);
hh2.EdgeColor = 'none';
plot((-lead_lag:lead_lag)/30,corr_coef8_m,'b','linewidth',2)

%plot correlation between WP mean MJO and TC during El Nino events
hh1 = fill([xx fliplr(xx)],...
    [corr_coef7_m-corr_coef7_std fliplr(corr_coef7_m+corr_coef7_std)],'r','facealpha',facealpha);
hh1.EdgeColor = 'none';
plot((-lead_lag:lead_lag)/30,corr_coef7_m,'r','linewidth',2)
box on
set(gca,'fontsize',14)
xlabel('(MJOI leads)  <-- lag -->  (T_C leads)')
ylabel('Corr')
title({'Model results';'Lagged correlation with T_C';'(a) All El Nino events'})
ylim([-0.2,0.4])
xlim([-12,12])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Observed lagged correlation %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Obs_fine_2020 %observed wind thermocline depth and sst

%subset of times we use from observational data
time2 = 1:13929; %time indicies for thermocline depth and sst
time3 = 1095+time2; %time indicies for MJO

% Finding the indices of the Pacific ocean domain in observational data sets
Left1_end = (120+2.5) / 2.5 * 10; % Pacific ocean left boundary 120E for SST
Right1_end = 280 / 2.5 * 10; % Pacific ocean right boundary 80W for SST
Middle1_end = 200 / 2.5 * 10; % Middle of Pacific ocean 160W for SST

Left1_C = 160/2.5*10; % CP left boundary for SST 160E
Right1_C = 210/2.5*10; % CP right boundary for SST 150W

Left2_end = (120+1) ; % Pacific ocean left boundary 120E for thermocline
Right2_end = 280; % Pacific ocean right boundary 80W for thermocline
Middle2_end = 200; % Middle of Pacific ocean 160W for thermocline

Left3_end = (120+2.5) / 2.5; % Pacific ocean left boundary 120E for wind bursts
Right3_end = 280 / 2.5; % Pacific ocean right boundary 80W for wind bursts
Middle3_end = 200 / 2.5; % Middle of Pacific ocean 160W for wind bursts

Left3_1 = (25+2.5)/2.5; % Left boundary of the Indian Ocean for MJO 25E


% SSTA Time series 
T_E_3R_obs = mean(sst_a_fine(time2,Middle1_end:Right1_end),2);% SST averaged over the eastern Pacific
T_C_3R_obs = mean(sst_a_fine(time2,Left1_C:Right1_C),2);% SST averaged over the Central Pacific

window = 30; %montly count for observed events
total_loop=floor(length(T_E_3R_obs)/window); %number of months in the dataset

% initializing count arrays for the different types of events
EN_obs=zeros(total_loop,1); % El Ninos
EPEN_obs=zeros(total_loop,1); % EP El Ninos
CPEN_obs=zeros(total_loop,1); % CP El Ninos
LN_obs=zeros(total_loop,1); %La Ninas

%going through each month in the dataset and populating a 1 if there was an
%event for each respective event
for k = 1:total_loop
    if mean(T_E_3R_obs(-0+k*window:+2+k*window)) > 0.5 || mean(T_C_3R_obs(-0+k*window:+2+k*window)) > 0.5
        EN_obs(k)=1;
        if mean(T_E_3R_obs(-0+k*window:+2+k*window))>mean(T_C_3R_obs(-0+k*window:+2+k*window))
            EPEN_obs(k)=1;
        else
            CPEN_obs(k)=1;
        end
    elseif mean(T_E_3R_obs(-0+k*window:+2+k*window)) < -0.5 || mean(T_C_3R_obs(-0+k*window:+2+k*window)) < -0.5
        LN_obs(k)=1;
    end
end


MJO_new2 = mean(MJO_obs(Left3_end:Middle3_end,time3),1);%MJO in the Western Pacific

MJO_index3 = MJO_new2; %MJO in the WP
MJO_index4 = abs(MJO_index3); %amplitude of MJO in the WP


% Extract MJO segments for each ENSO type
window_months = 6;%number of months we collect around the start and end of an event

EP_mjo1 = extract_mjo_segments(MJO_index3, EPEN_obs, window_months); %extract MJO in pacific during EP events
EP_mjo2 = extract_mjo_segments(MJO_index4, EPEN_obs, window_months); %extract amplitude of MJO in pacific during EP events
EP_sst1 = extract_mjo_segments(T_E_3R_obs', EPEN_obs, window_months); %extract SSTA in EP during EP events

EP_mjo3 = extract_mjo_segments(MJO_index3, EN_obs, window_months); %extract MJO in pacific during El Nino events
EP_mjo4 = extract_mjo_segments(MJO_index4, EN_obs, window_months); %extract amplitude of MJO in pacific during El Nino events
EP_sst2 = extract_mjo_segments(T_E_3R_obs', EN_obs, window_months); %extract SSTA in EP during El Nino events

CP_mjo1 = extract_mjo_segments(MJO_index3, CPEN_obs, window_months); %extract MJO in pacific during CP events
CP_mjo2 = extract_mjo_segments(MJO_index4, CPEN_obs, window_months); %extract amplitude of MJO in pacific during CP events
CP_sst1 = extract_mjo_segments(T_C_3R_obs', CPEN_obs, window_months); %extract SSTA in CP during CP events

CP_mjo3 = extract_mjo_segments(MJO_index3, EN_obs, window_months); %extract MJO in pacific during El Nino events
CP_mjo4 = extract_mjo_segments(MJO_index4, EN_obs, window_months); %extract amplitude of MJO in pacific during El Nino events
CP_sst2 = extract_mjo_segments(T_C_3R_obs', EN_obs, window_months); %extract SSTA in CP during El Nino events



lead_lag = 360;
window_lag_lead = lead_lag*2+1; % lead and lag are both 24 months

corr_coef1_obs = zeros(1,window_lag_lead); % MJO in Pacific ocean
corr_coef2_obs = zeros(1,window_lag_lead); % variance of MJO in Pacific ocean
corr_coef3_obs = zeros(1,window_lag_lead); % MJO in Pacific ocean
corr_coef4_obs = zeros(1,window_lag_lead); % variance of MJO in Pacific ocean
corr_coef5_obs = zeros(1,window_lag_lead); % MJO in Pacific ocean
corr_coef6_obs = zeros(1,window_lag_lead); % variance of MJO in Pacific ocean
corr_coef7_obs = zeros(1,window_lag_lead); % MJO in Pacific ocean
corr_coef8_obs = zeros(1,window_lag_lead); % variance of MJO in Pacific ocean

% computing the lagged correlation
j = 1;
for i = -lead_lag:lead_lag
    %MJO in Pacific ocean with EP SSTA during EP events
    temp=corrcoef(EP_mjo1(window_lag_lead+i:end-window_lag_lead+i), EP_sst1(window_lag_lead:end-window_lag_lead));
    corr_coef1_obs(j) = temp(1,2);
    %MJO variance in Pacific ocean with EP SSTA during EP events
    temp=corrcoef(EP_mjo2(window_lag_lead+i:end-window_lag_lead+i), EP_sst1(window_lag_lead:end-window_lag_lead));
    corr_coef2_obs(j) = temp(1,2);
    %MJO in Pacific ocean with EP SSTA during EN events
    temp=corrcoef(EP_mjo3(window_lag_lead+i:end-window_lag_lead+i), EP_sst2(window_lag_lead:end-window_lag_lead));
    corr_coef3_obs(j) = temp(1,2);
    %MJO variance in Pacific ocean with EP SSTA during EN events
    temp=corrcoef(EP_mjo4(window_lag_lead+i:end-window_lag_lead+i), EP_sst2(window_lag_lead:end-window_lag_lead));
    corr_coef4_obs(j) = temp(1,2);
    %MJO in Pacific ocean with CP SSTA during CP events
    temp=corrcoef(CP_mjo1(window_lag_lead+i:end-window_lag_lead+i), CP_sst1(window_lag_lead:end-window_lag_lead));
    corr_coef5_obs(j) = temp(1,2);
    %MJO variance in Pacific ocean with CP SSTA during CP events
    temp=corrcoef(CP_mjo2(window_lag_lead+i:end-window_lag_lead+i), CP_sst1(window_lag_lead:end-window_lag_lead));
    corr_coef6_obs(j) = temp(1,2);
    %MJO in Pacific ocean with CP SSTA during EN events
    temp=corrcoef(CP_mjo3(window_lag_lead+i:end-window_lag_lead+i), CP_sst2(window_lag_lead:end-window_lag_lead));
    corr_coef7_obs(j) = temp(1,2);
    %MJO varinace in Pacific ocean with CP SSTA during EN events
    temp=corrcoef(CP_mjo4(window_lag_lead+i:end-window_lag_lead+i), CP_sst2(window_lag_lead:end-window_lag_lead));
    corr_coef8_obs(j) = temp(1,2);
    j = j + 1;
end


%Plot observed on the same figure as the model
subplot(2,4,5)
hold on
plot((-lead_lag:lead_lag)/30,corr_coef2_obs,'b','linewidth',2)%MJO variance in Pacific ocean with EP SSTA during EP events
plot((-lead_lag:lead_lag)/30,corr_coef1_obs,'r','linewidth',1)%MJO in Pacific ocean with EP SSTA during EP events
box on
set(gca,'fontsize',14)
xlabel('(MJOI leads)  <-- lag -->  (T_E leads)')
ylabel('Corr')
title('(e) EP El Nino events')
ylim([-0.2,0.4])
xlim([-12,12])

subplot(2,4,1)
hold on
plot((-lead_lag:lead_lag)/30,corr_coef4_obs,'b','linewidth',2)%MJO variance in Pacific ocean with EP SSTA during EN events
plot((-lead_lag:lead_lag)/30,corr_coef3_obs,'r','linewidth',1)%MJO in Pacific ocean with EP SSTA during EN events
legend('Correlation using |MJOI|','Correlation using MJOI')
box on
set(gca,'fontsize',14)
xlabel('(MJOI leads)  <-- lag -->  (T_E leads)')
ylabel('Corr')
title({'Observations';'Lagged correlation with T_E';'(a) All El Nino events'})
ylim([-0.2,0.4])
xlim([-12,12])


subplot(2,4,7)
hold on
plot((-lead_lag:lead_lag)/30,corr_coef6_obs,'b','linewidth',2)%MJO variance in Pacific ocean with CP SSTA during CP events
plot((-lead_lag:lead_lag)/30,corr_coef5_obs,'r','linewidth',1)%MJO in Pacific ocean with CP SSTA during CP events
box on
set(gca,'fontsize',14)
xlabel('(MJOI leads)  <-- lag -->  (T_C leads)')
ylabel('Corr')
title('(g) CP El Nino events')
ylim([-0.2,0.4])
xlim([-12,12])


subplot(2,4,3)
hold on
plot((-lead_lag:lead_lag)/30,corr_coef8_obs,'b','linewidth',2)%MJO varinace in Pacific ocean with CP SSTA during EN events
plot((-lead_lag:lead_lag)/30,corr_coef7_obs,'r','linewidth',1)%MJO in Pacific ocean with CP SSTA during EN events
box on
set(gca,'fontsize',14)
xlabel('(MJOI leads)  <-- lag -->  (T_C leads)')
ylabel('Corr')
title({'Observations';'Lagged correlation with T_C';'(c) All El Nino events'})
ylim([-0.2,0.4])
xlim([-12,12])


%Function that takes a chunk of MJO or SST time series, a chunk of one of
%one of the event counters, and the number of months you want to capture
%around the start and end of an event. We use 6 months
function mjo_segments = extract_mjo_segments(MJO, event_index, window_months)
    event_starts = find(diff([0; event_index; 0]) == 1); %makes an array of each index where an event begins by testing where it goes from 0 to 1
    event_ends = find(diff([0; event_index; 0]) == -1) - 1; %makes an array of each index where an event ends by seeing where it goes from 1 to 0
    
    mjo_segments = [];%empty output vector to be populated
    for i = 1:length(event_starts) %goes through each event
        start_month = max(1, event_starts(i) - window_months);%takes 6 months before the start month unless an invalid index
        end_month = min(length(event_index), event_ends(i) + window_months);%takes 6 months after the end month unless and invalid index
        
        start_day = (start_month-1) * 30 + 1;%calculates the start day
        end_day = min(end_month * 30, size(MJO, 2));%calculates the end day
        
        mjo_segments = [mjo_segments, MJO(:, start_day:end_day)];%adds the new event segment to the previous
    end
end