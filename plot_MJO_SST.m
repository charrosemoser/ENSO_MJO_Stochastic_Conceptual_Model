%Plotting spatiotemporal patterns for MJO and ENSO from model simulation
load lon_sst.mat %array of longitudes over the pacific

x = [0:Na-1] * dx;%7 grid points for MJO in km

%Time snapshot used in figure 4.2
year_add = 10;
year_l=[112,122,132,142,152];

%Time snapshot used in figure S4 e
% year_add = 20;
% year_l=[490,510];

%Time snapshot used in figure S6
% year_add = 10;
% year_l = [60,70,80,90]+0;

%isolating the positive and negative interseasonal wind variation
temp = tau - u_interannual*50;
tau_plus = temp; tau_plus(temp<0) = 0; tau_plus = tau_plus+ u_interannual*50;
tau_minus = temp; tau_minus(temp>0) = 0; tau_minus = tau_minus+ u_interannual*50;

figure
for jj = 1:4
    window = 12;%yearly window
    total_loop = (year_add - 1) * 12/window;%number of windows (one year here) in one hovmoller diagram
    range_model = 1+12*year_l(jj):12*(year_l(jj)+year_add);%range for each panel
    
    range_tau = range_model(1)*k_dt:range_model(end)*k_dt; % the wind tau has a different temporal resolution and needs to be handled separately
    t_model = range_model/12; % time window for variables except tau
    tau_model = range_tau/k_dt/12; % time wind for tau
    Hov = zeros(length(lon_sst),length(range_model)); % the matrix is for the spatiotemporal pattern, namely the Hovmoller diagram
    
    xx = [T_C_3R(range_model),T_E_3R(range_model)]; % range of the longitude across equatorial Pacific
    % computing the Hovmoller diagram from the bivariate regression
    for i = 1:coarse_x_n
        Hov(i,:) = xx * SST_reg(i,:)'; 
    end

    % plot the Hovmoller diagram of SST
    subplot(1,8,(jj-1)*2+1)
    [xx,yy] = meshgrid(t_model,lon_sst);
    contourf(yy,xx,Hov,30,'linestyle','none');
    hold on
    temp_tau = range_tau(1:k_dt/5:end); % plot tau on top of SST Hovmoller diagram; zero-wind is at the dateline; negative value means easterly and positive value means westerly
    plot(180+u_interannual(temp_tau)*50*10,tau_model(1:k_dt/5:end),'k','linewidth',0.5);
    plot(180+tau_plus(temp_tau)*10,tau_model(1:k_dt/5:end),'r','linewidth',0.5);
    plot(180+tau_minus(temp_tau)*10,tau_model(1:k_dt/5:end),'b','linewidth',0.5);

    %plot colored boxes for each type of ENSO event
    for k = 1:total_loop
        if mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) > 1.0 && mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))>mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))
            plot([120,120],[t_model(1-4+k*window),t_model(1+1+window*k)],'r','linewidth',10)
        elseif mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) > 0.5 && mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))>mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))
            plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'m','linewidth',10)
        elseif mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) > 0.5 && mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))>mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window))
            plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'color',[255 97 0]/255,'linewidth',10)
        elseif mean(T_E_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) < -0.5 || mean(T_C_3R(range_model(1)-1+k*window:range_model(1)+1+k*window)) < -0.5
            plot([120,120],[t_model(1-4+window*k),t_model(1+1+window*k)],'b','linewidth',10)
        end
    end


    set(gca,'xlim',[120 280]); % x-range is the longitude
    set(gca,'xtick',120:60:280);
    


    title('SST')
    xlabel('Longitude');
    set(gca,'fontsize',14,'linewidth',1);
    caxis([-4,4])
    if jj == 1
        ylabel('model year');
    end
    colorbar('southoutside')
    
    %plot MJO panel for the same range of time
    temp = MJO(:,[range_tau(1):60:range_tau(end)])';
    [xx,yy] = meshgrid(x*dim_x*360/40000,[range_tau(1):60:range_tau(end)]*dt/6); 
    subplot(1,8,(jj-1)*2+2)
    hold on
    contourf(xx,yy,temp,40, 'linestyle','none')
    title('MJO')
    set(gca,'fontsize',14)
    xlabel('Longitude');
    colormap jet
    colorbar('southoutside')
    plot([120,120],ylim,'r','LineWidth',2) 
    xlim([30,280])
    caxis([-0.2,0.2])
end