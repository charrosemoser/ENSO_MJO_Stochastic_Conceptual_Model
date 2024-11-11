%Plotting the MJO variable power spectra first from the model and the
%observational data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Model Power Spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the spectrum of different variables in the MJO skeleton model
% The four variables are u, theta, Q and A (all are anomalies)

N2 = 300000;
if mod(N2,2) == 0
    zero_mode2 = (N2)/2+1;
else
    zero_mode2 = (N2+1)/2;
end

figure
for i = 1:4
    % choose the variable to plot the spectrum
    if i == 1
        variable = u_phy(:,100001:400000);
        variable = fftshift(fft(variable));
    elseif i == 2
        variable = theta_phy(:,100001:400000);
        variable = fftshift(fft(variable));
    elseif i == 3
        variable = Q(:,100001:400000);
    elseif i == 4
        variable = Hbar*A(:,100001:400000);
    end
    % define a matrix for computing the spectrum
    Spectrum_variable = zeros(size(variable));
    % along each longitude direction, compute the fft to get the spectrum
    for j = 1:Na    
        Spectrum_variable(j,:) = fftshift(fft(variable(j,:)'));
    end
    % the actual spectural is computed by taking the L2 norm, which is the
    % amplitude of the above fft coefficients
    Spectrum_variable = sqrt(Spectrum_variable .* conj(Spectrum_variable))/length(variable(1,:));
    % smooth out the solution to make it look slightly nicer
    for j = 1:Na
        Spectrum_variable(j,:) = smooth(Spectrum_variable(j,:));
    end

    % the range of plotting the result: the spatial modes
    % here, we keep modes +- 1, 2 and 3
    xscale = -3:3;

    % the range of plotting the result: the temporal direction
    temporal_freq_fastest = length(variable(1,:))/2 - 1;  % the largest frequency, meaning the fastest oscillation in time
    yscale = [1:temporal_freq_fastest]/(N2*dt*dim_t); % frequency in time from lowest to highest
    [xx_scale,yy_scale] = meshgrid(xscale,yscale); % generating the meshgrid
    % plot the result
    subplot(2,2,i)
    hold on
    contourf(xx_scale, yy_scale, log(Spectrum_variable(zero_mode+xscale, zero_mode2 + [1:temporal_freq_fastest])'),40, 'linestyle','none');
    % together with the theoretic values from the linear solution 
    plot(0:5, Eig_Store(1:2,1:6)/2/pi/dim_t,'ko','linewidth',2)
    plot(-[0:5], -Eig_Store(3:4,1:6)/2/pi/dim_t,'ko','linewidth',2)
    xlim([-3,3])
    ylim([0,0.1])
    % band of the 30 to 90 days
    plot([-3,3],[1/30,1/30],'--k','linewidth',2)
    plot([-3,3],[1/90,1/90],'--k','linewidth',2)
    plot([0,0],[0,0.1],'--k','linewidth',2)
    if i == 1
        title('Zonal velocity')
    elseif i == 2
        title('Temperature')
    elseif i == 3
        title('Moisture')
    elseif i == 4
        title('Convective activiy')
    end
    colorbar
end
colormap jet
sgtitle({'Model Power Spectra'},'fontsize',14);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Plot observed MJO power spectra %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load hgt_new_data.mat % geopontential height data
load uwnd_new_data %zonal velocity data
load olr_new_data %outgoing longwave radiation - used as a surrogate for convective activity

%observed Kelvin and Rossby waves
K_obs = (uwnd_mode_0_rmseason - hgt_mode_0_rmseason)/2;
R_obs = -(uwnd_mode_0_rmseason+ hgt_mode_0_rmseason)/4 + (uwnd_mode_2_rmseason - hgt_mode_2_rmseason)/2/sqrt(2);
psi_0 = sqrt(2) * pi^(-1/4); % meridional basis psi_0 at equator
psi_2 = -(4*pi)^(-1/4); % meridional basis psi_2 at equator

%observed convective activity and zonal velocity
A_obs = olr_mode_0_rmseason;
u_phy_obs = ( (K_obs - R_obs) * psi_0 + 1/sqrt(2) * R_obs * psi_2 ) * 50; % unit 50m/s
A_phy_obs = A_obs* psi_0 * 15;

%using only first 2000 timesteps for u and A
u_phy_obs = u_phy_obs(:,1:2000);
A_phy_obs = A_phy_obs(:,1:2000);

%locating the zero mode in time
jj = size(u_phy_obs,2);
if mod(jj,2) == 0
    zero_mode2 = (jj)/2+1;
else
    zero_mode2 = (jj+1)/2;
end

%locating the zero mode in space
jj = size(u_phy_obs,1);
if mod(jj,2) == 0
    zero_mode1 = jj/2+1;
else
    zero_mode1 = (jj+1)/2;
end

figure
for i = 1:2 
    % choose the variable to plot the spectrum
    if i == 1
        variable = u_phy_obs;
        variable = fftshift(fft(variable));%fourier transform over space
    elseif i == 2
        variable = A_phy_obs;
        variable = fftshift(fft(variable));%fourier transform over space
    end
    Spectrum_variable = zeros(size(variable));
    for j = 1:size(u_phy_obs,1)    
        Spectrum_variable(j,:) = fftshift(fft(variable(j,:)'));%fourier transform over time
    end
    % the actual spectural is computed by taking the L2 norm, which is the
    % amplitude of the above fft coefficients
    Spectrum_variable = sqrt(Spectrum_variable .* conj(Spectrum_variable))/length(variable(1,:));
    % smooth out the solution to make it look slightly nicer
    for j = 1:size(u_phy_obs,1)
        Spectrum_variable(j,:) = smooth(Spectrum_variable(j,:));%smooth the spectrum
    end
    
    % the range of plotting the result: the spatial modes
    % here, we keep modes +- 1, 2 and 3
    time_modes = round(length(variable(1,:))/10);%time modes
    xscale = -3:3; %wave number range
    yscale = (1:time_modes)/(length(variable(1,:)));%frequency range
    [xx_scale,yy_scale] = meshgrid(xscale,yscale);
    %plot the spectrum for the given variable
    subplot(1,2,i)
    hold on
    contourf(xx_scale, yy_scale, log(Spectrum_variable(zero_mode1+xscale, zero_mode2 + (1:time_modes))'),40, 'linestyle','none');
    %plot dispersion curve calculated in Eigen_Solver
    plot(0:5, Eig_Store(1:2,1:6)/2/pi/dim_t,'ko','linewidth',2)
    plot(-[0:5], -Eig_Store(3:4,1:6)/2/pi/dim_t,'ko','linewidth',2)
    xlim([-3,3])
    ylim([0,0.1])
    %plot the lines where frequency correspond to 30 and 90 days
    plot([-3,3],[1/30,1/30],'--k','linewidth',2)
    plot([-3,3],[1/90,1/90],'--k','linewidth',2)
    %plot a line marking wave number 0
    plot([0,0],[0,0.1],'--k','linewidth',2)
    box on
    set(gca,'linewidth',1)
    if i == 1
        title('Zonal velocity')
    elseif i == 2
        title('Convective activity')
    end
    colorbar
end
colormap jet
sgtitle({'Observed Power Spectra'},'fontsize',14);

