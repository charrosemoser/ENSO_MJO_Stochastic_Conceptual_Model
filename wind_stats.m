%Plotting wind statistics from the model simulation and observation. This
%file contains code to reconstruct the interseasonal and total wind over
%the WP, to plot the PDFs of total wind and intereasonal wind for model
%simulation and observation, and to plot the time series comparing
%interseasonal wind interannual wind and SST in the EP from model
%simulation.

%load data
load uwnd_new_data.mat % wind data
load hgt_new_data.mat % geopotential height data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Reconstruct wind from observations %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Finding the indices of the Pacific ocean domain in observational data sets
Left3_end2 = (140+2.5) / 2.5; % Western Pacific ocean left boundary 140E for wind bursts
Middle3_end2 = 180 / 2.5; % Western Pacific ocean right boundary 180W for wind bursts

%Calculate K and R to reconstruct the observational intraseasonal wind and
%total wind
psi_0 = sqrt(2) * pi^(-1/4); % meridional basis psi_0 at equator
psi_2 = -(4*pi)^(-1/4); % meridional basis psi_2 at equator

%total observed Kelvin and Rossby waves
K_total = (uwnd_mode_0_rmmean_3modes - hgt_mode_0_rmmean_3modes)/2;
R_total = -(uwnd_mode_0_rmmean_3modes + hgt_mode_0_rmmean_3modes)/4 + (uwnd_mode_2_rmmean_3modes - hgt_mode_2_rmmean_3modes)/2/sqrt(2);

%observed anomaly Kelvin and Rossby waves
K_3modes = (uwnd_mode_0_rmseason_3modes - hgt_mode_0_rmseason_3modes)/2;
R_3modes = -(uwnd_mode_0_rmseason_3modes + hgt_mode_0_rmseason_3modes)/4 + (uwnd_mode_2_rmseason_3modes - hgt_mode_2_rmseason_3modes)/2/sqrt(2);

%total wind observed
u_total_obs = ( (K_total - R_total) * psi_0 + 1/sqrt(2) * R_total * psi_2 ) * 50; % unit 50m/s

%wind anomaly observed
u_phy_obs = ( (K_3modes - R_3modes) * psi_0 + 1/sqrt(2) * R_3modes * psi_2 ) * 50; % unit 50m/s

%average total wind observed over WP
u_Wtotal = mean(u_total_obs(Left3_end2:Middle3_end2,:),1);

%average wind anomaly observed over WP
u_W = mean(u_phy_obs(Left3_end2:Middle3_end2,:),1); 

%decompose wind into mean and positive and negative anomalies
uwnd_mean = u_Wtotal-u_W;
uwnd_plus = u_W; uwnd_plus(uwnd_plus<0) = 0; uwnd_plus = uwnd_plus+ uwnd_mean;
uwnd_minus = u_W; uwnd_minus(uwnd_minus>0) = 0; uwnd_minus = uwnd_minus+ uwnd_mean;

%Calculate total wind burst
Total_WB = -uwnd_mean + uwnd_minus + uwnd_plus; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Plot the PDFs of wind %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate the PDFs of intraseasonal wind from observation and model simulation
[fi_ISW,xx_ISW] = ksdensity(u_W); fi_ISW = smooth(fi_ISW);
[fi2_ISW,xx2_ISW] = ksdensity(intra_wind); 

%calculate the PDF of the modeled intraseasonal wind for each set of 40 years to construct the confidence interval 
PDF_intrawind_model_UQ = zeros(L_realization,100);
for i = 1:L_realization
    temp = ksdensity(intra_wind((i-1)*480*50+1: i*480*50),xx2_ISW);
    PDF_intrawind_model_UQ(i,:) = temp;     
end

%calculate the PDFs of tau from observation and model simulation
[fi,xx_tau] = ksdensity(Total_WB); fi = smooth(fi);
[fi2,xx2_tau] = ksdensity(tau); 

%calculate the PDF of tau for each set of 40 years to construct the confidence interval 
PDF_tau_model_UQ = zeros(L_realization,100);
for i = 1:L_realization
    temp = ksdensity(tau((i-1)*480*50+1: i*480*50),xx2_tau);
    PDF_tau_model_UQ(i,:) = temp;     
end

figure
hold on
%plot the PDFs of intraseasonal wind from observation and model simulation
plot(xx_ISW,fi_ISW,'r','LineWidth',1.5)
plot(xx2_ISW,fi2_ISW,'b','LineWidth',1.5)
%upper and lower bounds of the CI
upper_lim = mean(PDF_intrawind_model_UQ) + 2 * std(PDF_intrawind_model_UQ);
lower_lim = mean(PDF_intrawind_model_UQ) - 2 * std(PDF_intrawind_model_UQ); lower_lim(lower_lim<=0) = 0;
patch([xx2_ISW, xx2_ISW(end:-1:1)],[lower_lim,upper_lim(end:-1:1)],'b','facealpha',0.2,'linestyle','none');
box on
title('PDFs of intraseaonal wind')
legend('Observation','Model')
set(gca,'Fontsize',14)

figure
hold on
%plot the PDFs of tau from observation and model simulation
plot(xx_tau,fi,'r','LineWidth',1.5)
plot(xx2_tau,fi2,'b','LineWidth',1.5)
%upper and lower bounds of the CI
upper_lim = mean(PDF_tau_model_UQ) + 2 * std(PDF_tau_model_UQ);
lower_lim = mean(PDF_tau_model_UQ) - 2 * std(PDF_tau_model_UQ); lower_lim(lower_lim<=0) = 0;
patch([xx2_tau, xx2_tau(end:-1:1)],[lower_lim,upper_lim(end:-1:1)],'b','facealpha',0.2,'linestyle','none');
box on
title('PDFs of \tau')
legend('Observation','Model')
set(gca,'Fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Comparison of the wind time series %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_obs = 1980+1/12: 1/12: 2019;
LL = length(t_obs); % let the length of the time series from the model be the same as observations 
range_model = 1+12*362: LL+12*362;
range_tau = range_model(1)*k_dt:range_model(end)*k_dt;
t_model = range_model/12; % time window of plotting variables except the wind 
tau_model = range_tau/k_dt/12;
tau_temp=tau;
u_interannual_temp=u_interannual*50;

%plot the time series fo SSTA, interannual wind, and interseasonal wind
figure
subplot(3,1,1)
plot(t_model,T_E_3R(range_model),'r','linewidth',2);title('(a) T_E')
box on;set(gca,'fontsize',14)
xlim([t_model(1),t_model(end)])
subplot(3,1,2)
plot(tau_model,u_interannual_temp(range_tau),'k','linewidth',2);title('(b) Interannual wind')
box on;set(gca,'fontsize',14)
xlim([tau_model(1),tau_model(end)])
subplot(3,1,3)
plot(tau_model,(tau_temp(range_tau)-u_interannual_temp(range_tau)),'g','linewidth',2);title('(c) Intraseasonal wind')
box on;set(gca,'fontsize',14)
xlim([tau_model(1),tau_model(end)])
