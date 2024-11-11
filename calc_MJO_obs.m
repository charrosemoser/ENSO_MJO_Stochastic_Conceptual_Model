%Calculating the MJO signal from Observational data

load hgt_new_data %geopolential height data
load uwnd_new_data %zonal wind data
load Q_new_data %moisture data
load olr_new_data %outgoing longwave radiation data

% define kelvin and rossby waves using zonal wind and geopotential height
K_rmseason = (uwnd_mode_0_rmseason - hgt_mode_0_rmseason)/2; 
R_rmseason = -(uwnd_mode_0_rmseason + hgt_mode_0_rmseason)/4 + (uwnd_mode_2_rmseason - hgt_mode_2_rmseason)/sqrt(2)/2;

%we use outgoing longwave radiation as a surrogate for convective activity A
A_rmseason = olr_mode_0_rmseason;

%define moisture
Q_rmseason = Q_mode_0_rmseason;

% Take the data, in which the annual mean, first three harmonic and 120-day
% moving averaged are removed, and only zonal wave numbers 1, 2, 3 (and -1, -2, -3) are kept.
K_3modes = (uwnd_mode_0_rmseason_3modes - hgt_mode_0_rmseason_3modes)/2;
R_3modes = -(uwnd_mode_0_rmseason_3modes + hgt_mode_0_rmseason_3modes)/4 + (uwnd_mode_2_rmseason_3modes - hgt_mode_2_rmseason_3modes)/sqrt(2)/2;
A_3modes = olr_mode_0_rmseason_3modes;
Q_3modes = Q_mode_0_rmseason_3modes;


% lt is the total time length of the data from 1979-2011
lt = size(K_3modes,2);
t = 1:lt;
t = transpose(t);

% The eigenvectors from skeleton model, corresponding to the MJO. 
% Each variable has three entries, representing that of mode 1-3, respectively.
eg_K = eg_K_MJO;
eg_R = eg_R_MJO;
eg_Q = eg_Q_MJO;
eg_A = eg_A_MJO;

% Frequencies corresponding to the MJO (mode 1-3)
omega = [0.0219, 0.0239, 0.0238];

%Fourier transforms across space
K2 = fftshift(fft(K_3modes, [], 1), 1);
R2 = fftshift(fft(R_3modes, [], 1), 1);
Q2 = fftshift(fft(Q_3modes, [], 1), 1);
A2 = fftshift(fft(A_3modes, [], 1), 1);

%identify which spatial index has the zero mode
zero_mode1 = (144)/2+1;

%identify which time index has the zero mode
if mod(lt,2) == 0
    zero_mode2 = (lt)/2+1;
else
    zero_mode2 = (lt+1)/2;
end

%initialize variables for fourier transform of MJO and MJO
MJO1_F_obs = zeros(144, 1); MJO1_temp = zeros(144,lt); %modes 1,2,3
MJO2_F_obs = zeros(144, 1); MJO2_temp = zeros(144,lt); %model -1,-2,-3
MJO1_obs = zeros(144, lt); %modes 1,2,3
MJO2_obs = zeros(144, lt); %model -1,-2,-3

%loop through time to compute MJO components for each time step
for i = 1:lt
    %Fourier transform of MJO of modes 1,2,3
    MJO1_F_obs(zero_mode1+1) =  K2(zero_mode1+1,i) * eg_K(1)' + R2(zero_mode1+1,i) * eg_R(1)' + Q2(zero_mode1+1,i) * eg_Q(1)' + A2(zero_mode1+1,i) * eg_A(1)';
    MJO1_F_obs(zero_mode1+2) =  K2(zero_mode1+2,i) * eg_K(2)' + R2(zero_mode1+2,i) * eg_R(2)' + Q2(zero_mode1+2,i) * eg_Q(2)' + A2(zero_mode1+2,i) * eg_A(2)';
    MJO1_F_obs(zero_mode1+3) =  K2(zero_mode1+3,i) * eg_K(3)' + R2(zero_mode1+3,i) * eg_R(3)' + Q2(zero_mode1+3,i) * eg_Q(3)' + A2(zero_mode1+3,i) * eg_A(3)';
    %MJO modes 1,2,3 in physical space
    MJO1_temp(:,i) = ifft(ifftshift(MJO1_F_obs)); 
    
    %Fourier transform of MJO of modes -1,-2,-3
    MJO2_F_obs(zero_mode1-1) =  K2(zero_mode1-1,i) * eg_K(1) + R2(zero_mode1-1,i) * eg_R(1) + Q2(zero_mode1-1,i) * eg_Q(1) + A2(zero_mode1-1,i) * eg_A(1);
    MJO2_F_obs(zero_mode1-2) =  K2(zero_mode1-2,i) * eg_K(2) + R2(zero_mode1-2,i) * eg_R(2) + Q2(zero_mode1-2,i) * eg_Q(2) + A2(zero_mode1-2,i) * eg_A(2);
    MJO2_F_obs(zero_mode1-3) =  K2(zero_mode1-3,i) * eg_K(3) + R2(zero_mode1-3,i) * eg_R(3) + Q2(zero_mode1-3,i) * eg_Q(3) + A2(zero_mode1-3,i) * eg_A(3);

    %MJO modes -1,-2,-3 in physical space
    MJO2_temp(:,i) = ifft(ifftshift(MJO2_F_obs));

end

%temporal filtering 
bd1 = -90; bd2 = -30; %periods we are concerned with in the modes 1-3

%loop through the spatial grid and filter
for i = 1:144
    MJO1_temporal = fftshift(fft(MJO1_temp(i,:)));%fourier transform across time
    MJO1_temporal2 = MJO1_temporal * 0;
    MJO1_temporal2(zero_mode2+round(lt/bd2): zero_mode2+round(lt/bd1)) = MJO1_temporal(zero_mode2+round(lt/bd2): zero_mode2+round(lt/bd1)); %filter out any signal that does not have a period between -30 and -90 days
    MJO1_obs(i,:) = real(ifft(ifftshift(MJO1_temporal2)));% inverse fourier transform 
end

bd1 = 30; bd2 = 90;%periods we are interested in for modes -1,-2,-3

%loop through the spatial grid and filter
for i = 1:144
    MJO2_temporal = fftshift(fft(MJO2_temp(i,:)));%fourier transform accross time
    MJO2_temporal2 = MJO2_temporal * 0;
    MJO2_temporal2(zero_mode2+round(lt/bd2): zero_mode2+round(lt/bd1)) = MJO2_temporal(zero_mode2+round(lt/bd2): zero_mode2+round(lt/bd1));%filter out any signal that does not have a period between 30 and 90 days
    MJO2_obs(i,:) = real(ifft(ifftshift(MJO2_temporal2)));%inverse fourier transform
end
%note we take the real part of MJO1 and MJO2 because they are complex
%conjugates and the complex part should cancel. If any complex part is left
%it is due to numerical error, so we just take the real part.
MJO_obs = MJO1_obs+MJO2_obs;
