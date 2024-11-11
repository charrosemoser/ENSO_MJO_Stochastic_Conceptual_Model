%Plotting spatiotemporal results of MJO and SST from observational data 
load Obs_fine_2020 %SST, wind, and thermocline data
load uwnd_new_data.mat % wind data
load hgt_new_data.mat % geopotential height data
%subset of times we use from observational data. Since the two data sets
%start at different times we use different indices for each so that they
%line up.
time2 = (1982-1982)*365+1:(1987-1982)*365+1; %use these indices for wind, SST, and thermocline
time3 = 1095+time2; %use these indicies for MJO

ttt = Y; % time variable same indicies as SST thermocline and wind

% Finding the indices of the Pacific ocean domain in observational data sets
Left1_end = (120+2.5) / 2.5 * 10; % Pacific ocean left boundary 120E for SST
Right1_end = 280 / 2.5 * 10; % Pacific ocean right boundary 80W for SST

Right3_end = 280 / 2.5; % Pacific ocean right boundary 80W for wind bursts
Left3_end2 = (140+2.5) / 2.5; % Western Pacific ocean left boundary 140E for wind bursts
Middle3_end2 = 180 / 2.5; % Western Pacific ocean right boundary 180W for wind bursts
Left3_1 = (25+2.5)/2.5; % Left boundary of the Indian Ocean for MJO 25E

psi_0 = sqrt(2) * pi^(-1/4); % meridional basis psi_0 at equator
psi_2 = -(4*pi)^(-1/4); % meridional basis psi_2 at equator

%observed total Kelvin and Rossby waves 
K_total = (uwnd_mode_0_rmmean_3modes - hgt_mode_0_rmmean_3modes)/2;
R_total = -(uwnd_mode_0_rmmean_3modes + hgt_mode_0_rmmean_3modes)/4 + (uwnd_mode_2_rmmean_3modes - hgt_mode_2_rmmean_3modes)/2/sqrt(2);

%observed anomaly of Kelvin and Rossby waves 
K_3modes = (uwnd_mode_0_rmseason_3modes - hgt_mode_0_rmseason_3modes)/2;
R_3modes = -(uwnd_mode_0_rmseason_3modes + hgt_mode_0_rmseason_3modes)/4 + (uwnd_mode_2_rmseason_3modes - hgt_mode_2_rmseason_3modes)/2/sqrt(2);

%observed wind total
u_total_obs = ( (K_total - R_total) * psi_0 + 1/sqrt(2) * R_total * psi_2 ) * 50; % unit 50m/s

%observed wind anomaly
u_phy_obs = ( (K_3modes - R_3modes) * psi_0 + 1/sqrt(2) * R_3modes * psi_2 ) * 50; % unit 50m/s

%average observed wind over Western Pacific total
u_Wtotal = mean(u_total_obs(Left3_end2:Middle3_end2,:),1); 

%average observed wind over Western Pacific anomaly
u_W = mean(u_phy_obs(Left3_end2:Middle3_end2,:),1);

%decompose wind into mean and positive and negative anomalies
uwnd_mean = u_Wtotal-u_W;
uwnd_plus = u_W; uwnd_plus(uwnd_plus<0) = 0; uwnd_plus = uwnd_plus+ uwnd_mean;
uwnd_minus = u_W; uwnd_minus(uwnd_minus>0) = 0; uwnd_minus = uwnd_minus+ uwnd_mean;

year_l = [1982,1992,2002,2008];%start for each set of 10 years

%plot the 4 ten year periods
figure
colormap jet
for jj = 1:4
    %defining the 10 year period we are working with
    time2 = (year_l(jj)-1982)*365+1:(year_l(jj)+10-1982)*365;%for SST hW and wind
    time3 = 1095+time2; %for MJO
  
    %plot MJO over indian and pacific ocean
    subplot(1,8,(jj-1)*2+2)
    [xx,yy] = meshgrid(X_uwnd(Left3_1:Right3_end),ttt(time2));
    contourf(xx,yy,MJO_obs(Left3_1:Right3_end,time3)','LineStyle','none')
    set(gca,'FontSize',14)
    title('MJO')
    cb = colorbar('southoutside');
    hold on
    plot([120 120], ylim, 'r-', 'LineWidth', 2)
    hold off
    clim([-0.25,0.25])
    set(gca,'xtick',50:100:250);
    xlabel('Longitude')
    box on
    
    %plot SST and wind
    subplot(1,8,(jj-1)*2+1)
    hold on
    [xx,yy] = meshgrid(X_sst(Left1_end:Right1_end),ttt(time2));
    contourf(xx,yy,sst_a_fine(time2,Left1_end:Right1_end),'LineStyle','none')%plot SST
    plot(180 + 8*uwnd_mean(time3),ttt(time2),'k','LineWidth',1)%plot interannual wind
    plot(180 + 8*uwnd_plus(time3),ttt(time2),'r','LineWidth',1)%plot positive interseasonal variation
    plot(180 + 8*uwnd_minus(time3),ttt(time2),'b','LineWidth',1)%plot negative interseasonal variation
    xlim([X_sst(Left1_end),X_sst(Right1_end)]);
    set(gca,'FontSize',14)
    title('SST')
    cb = colorbar('southoutside');
    clim([-4,4])
    set(gca,'xlim',[120 280]); % x-range is the longitude
    set(gca,'xtick',120:60:280); 
    xlabel('Longitude')
    box on
end
