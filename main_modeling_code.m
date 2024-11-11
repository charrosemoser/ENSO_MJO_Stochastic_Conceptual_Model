% Conceptual model for ENSO is based on Chen, Fang & Yu (2022) npj paper.
% Stochastic Skeleton Model for the MJO (Modified from Thual et al JAS
% paper). The difference is that in this code, the noise is not generated
% by a Markov jump process. Instead, a stochastic process with
% multiplicative noise is adopted. Mathematically, we can show the
% equivalence in the limit when dt goes to zero. See Chen & Majda 2015 MWR
% paper.

%Find the eigen vectors and values from the MJO skeleton model
Eigen_Solver 
%Find the regression coefficients for the CP and EP
SST_bivar

rng(91) % fix the random number seed to reproduce the results 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% One time unit in the model is two months %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0.001; % time step; dt = 0.001 is roughly 60*0.001*24 = 1.44 hours.
T = 12000; % total time length, i.e., 12000*2 = 24000 months = 2000 years.
N = round(T/dt); % total number of numerical integration steps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  The setup of the decadal variability I  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The decadal variability is represented by a linear SDE with
% multiplicative noise.
% Prescribe a PDF for the decadal variability I. It is a uniform distribution
% within the window I in [0,4].
x_gap = 0.1; % discritization 
xx = 0:x_gap:4; % range of the random variable
p = ones(size(xx)) * 1/4; % PDF of the random variable
m = trapz(xx,xx.*p); % mean value of the random variable

lambda = 2/60;% the damping in the linear SDE with multiplicative noise
% The inverse of labmda, i.e., 1/lambda is roughly the decorrelation time
% Here, it is 60/2=30 time units. Since each time unit is 2 months. The
% memory time of the decadal variable is about 5 years.

% Determine the multiplicative noise
n = length(xx);
Phi = zeros(1,n);
sgm = zeros(1,n);
for i = 2:n
    Phi(i) = trapz(xx(1:i), (xx(1:i)-m) .* p(1:i));
    sgm(i) = real(sqrt(2/p(i) * (-lambda * Phi(i))));
end

% Simulate the decadal variability using the linear model with
% multiplicative noise
I = zeros(1,N);
for i = 2:N
    temp = round(I(i-1)/x_gap) + 1;
    % Two boundary conditions just for the numerical purpose since a finite
    % dt is used
    if temp < 0
        temp = 1;
    elseif temp > n
        temp = n;
    end
    sgm_x = sgm(temp);
    % Linear model with multiplicative noise; numerical integration
    I(i) = I(i-1) + (-lambda * (I(i-1) - m) * dt) + sgm_x * randn *sqrt(dt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  Setup of the ENSO model  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters of the deterministic part (see the npj paper for details)
c = 1; 
gamma = 0.45; 
r = 0.15;  
alpha_2 = 0.075; 
alpha_1 = 0.0225;  
b_0 = 2.5;
mu = 0.5; 

Cu = 0.009; % the constant factor to guarantee the mean states of the variables are zero
d_tau = 2; % damping in the wind process; 1/d_tau = 0.5 = 1 month is the decorrelation time


% State variable
% WP: Western Pacific; CP: Central Pacific; EP: Eastern Pacific; 
u_3R = zeros(1,N); % ocean current velocity in CP
h_W_3R = zeros(1,N); % thermocline depth in WP
T_C_3R = zeros(1,N); % SST in CP
T_E_3R = zeros(1,N); % SST in EP
tau = zeros(1,N); % atmospheric wind



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Setup of the MJO model  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Spatial grid setup
dn = 1; % a factor that can refine the grids
Na = 7*dn; % atmosphere grid points; 
total_x = 40000; % total length of the equator (km)
dim_x = 15000; % one dimension unit of x axis (km) 
L = total_x/dim_x; % total units in the equator, which is 8/3
dx = L/Na; % distance between every two grid points
x = [0:Na-1] * dx; % grid points in the x axis

%Dimensional units
dim_t = 60; % one dimension unit of time  
dim_u = 50; % dimensional unit of velocity
dim_Q = 15; % dimensional unit of Q
dim_theta = 15; % dimensional unit of theta
dim_A = 15; % dimensional unit of A

%Parameters
Gamma = 4980/17; % Gamma parameter in skeleton model
Qbar = 0.9; % Q_bar in skeleton model
d_k = 30/17; % atmosphere damping (this is optional); if d_k = 1, that means the damping time is 2 months.
psi_0 = sqrt(2) * pi^(-1/4); % meridional basis psi_0 at equator
psi_2 = -(4*pi)^(-1/4); % meridional basis psi_2 at equator
Hbar = 660/17; % H bar parameter
lambda_A = 2; % one month

% warm pool background 
s_q = 66*(1+0.6*cos(2*pi*(1/L*x-1/4))')/17;  
s_theta = 66*(1+.6*cos(2*pi*(1/L*x-1/4))')/17;  

%Fourier wave numbers depending on if Na is even or odd
if mod(Na,2) == 0
    k = 2 * pi / (Na * dx) * (-Na/2:Na/2-1); 
else
    k = 2 * pi / (Na * dx) * (-(Na-1)/2:(Na-1)/2); 
end

% zero-th Fourier mode, also depending on if Na is even or odd
if mod(Na,2) == 0
    zero_mode = Na/2+1;
else
    zero_mode = (Na+1)/2;
end
 
 
% initial values; cannot all be zero since otherwise the multiplicative
% noise will not be triggered
% All variables with _phy means the variable in the physical space. Those
% without _phy means the variable in the Fourier space
% Z is the auxiliary variable to compute Q
%Eq is the latent heat
K = zeros(Na,N);
R = zeros(Na,N);
K_phy = zeros(Na,N);
R_phy = zeros(Na,N);
Q = zeros(Na,N); Q(Na,1) = 0.0;
Q_phy = zeros(Na,N); Q_phy(:,1) = 0.1;
A = zeros(Na,N); A(Na,1) = 0.001;
A_phy = zeros(Na,N); A_phy(:,1) = 0.001;
Z = zeros(Na,N);
Z_phy = zeros(Na,N);
u_phy = zeros(Na,N);
A_bar = zeros(Na,N);
u_interannual = zeros(1,N);
Eq_all = zeros(1,N);

% The variables with _old are used as the current value in the numerical
% integration to solve the value after dt, which will be defined as
% variables with _new
K_old = K(:,1);
R_old = R(:,1);
Q_old = Q(:,1);
A_old = A(:,1);
Z_old = Z(:,1);
K_phy_old = K_phy(:,1);
R_phy_old = R_phy(:,1);
Q_phy_old = Q_phy(:,1);
A_phy_old = A_phy(:,1);
Z_phy_old = Z_phy(:,1);
noise_Q_all = zeros(1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Numerical integration of the model  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load lon_sst %longitudes corresponding to the pacific ocean
for i = 2: N
    if mod(i-1,100000) == 0
        disp(['Steps computed = ', num2str(i-1), ', which is ', num2str((i-1) * dt / 6 ), ' years', '    Steps in total = ', num2str(N), ', which is ', num2str((N) * dt / 6 ), ' years']);
    end
    sigma = I(i-1)*0.12;
    
    
    % Quadratic damping in the CP area to help the PDF in CP
    c1 = 25 * (T_C_3R(i-1) + 0.1).^2 + 0.6;   
    c1 = c1 * (1 + 0.6 * sin(i*dt*2*pi/6) );
    c1 = c1  * 0.78;
 
    % Wind noise coefficients
    sigma_tau = 0.9* (tanh(4.5 * T_C_3R(i-1)) + 1 ) * (1 + 0.3 * cos(i * dt * 2 * pi / 6) );
    beta_E = 2 - I(i-1)/5;
    beta_E = beta_E * 4/3 * sqrt(0.6); 
    beta_u = -0.2 * beta_E;
    beta_h = -0.4 * beta_E;
    beta_C =  0.8 * beta_E;
    
    % Linear damping in the EP area
    c2 = 0.84 * (1 + 0.4 * sin(i * dt * 2 * pi / 6 + 2 * pi / 6) + 0.4 * sin(i * dt * 4 * pi / 6 + 2 * pi / 6) );

    % governing ocean equations 
    % Note that the decadal variable I has been pre-determined. The values
    % are directly used here. The evolution of I was in the above module. 
    u_3R(i) = u_3R(i-1) + ( - r * u_3R(i-1) - alpha_1 * b_0 * mu / 2 * (T_C_3R(i-1) + T_E_3R(i-1)) ) * dt + beta_u * tau(i-1) * dt;
    h_W_3R(i) = h_W_3R(i-1) + ( - r * h_W_3R(i-1) - alpha_2 * b_0 * mu / 2 * (T_C_3R(i-1) + T_E_3R(i-1)) ) * dt  + beta_h * tau(i-1) * dt;
    T_C_3R(i) = T_C_3R(i-1) + ( (gamma * b_0 * mu / 2 - c1) * T_C_3R(i-1) + gamma * b_0 * mu / 2 * T_E_3R(i-1) + gamma * h_W_3R(i-1) + sigma * u_3R(i-1) + Cu*2 ) * dt + beta_C * tau(i-1) * dt;
    T_E_3R(i) = T_E_3R(i-1) + ( gamma * h_W_3R(i-1) + (1 * gamma * b_0 * mu / 2 +gamma * b_0 * mu - c2) * T_E_3R(i-1) + (1 * gamma * b_0 * mu / 2 - gamma * b_0 * mu) * T_C_3R(i-1)) * dt + beta_E * tau(i-1) * dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Governing equations for the MJO: starting here
    alpha_q = 0.9;
    Eq = alpha_q * (T_C_3R(i-1));
    Eq_all(i) = Eq;
    A_bar_mean = 1 / Hbar / (1 - Qbar) * (Eq + s_q - Qbar * s_theta); % mean state of A;
    A_bar_mean(A_bar_mean <= 0) = 1e-5; % positivity; 
    A_bar(:,i) = A_bar_mean;
    % Note: the positivity is automatically satisfied when MJO appears without ENSO. 
    % However, when MJO is coupled with ENSO, ENSO will modulate A_bar. Then this positivity condition needs to be imposed.
    % Note 2: Usually, A_total is positive (which has to be satisfied), as will be seen below. 
    % A_bar being positive is a more strict requirement, which is still physically explainable.

    % update K in Fourier space    
    force = -1/2 * Hbar * A_phy_old;
    c_k = 300/17;
    K_old = fftshift(fft(K_phy_old));
    force = fftshift(fft(force)); 
    
    K_new = K_old .* exp( - d_k * dt - 1i * k' * c_k * dt) + ( force ./ ( d_k + 1i * k' * c_k) ) .* ( 1 - exp( - d_k * dt - 1i * k' * c_k * dt ) );
    if d_k == 0
        K_new(zero_mode) = K_old(zero_mode) + force(zero_mode) * dt;
    else
        K_new(zero_mode) = K_old(zero_mode) * exp( - d_k * dt ) + force(zero_mode) / d_k * ( 1 - exp( - d_k * dt ) );
    end
    K_phy_new = ifft(ifftshift(K_new));
    K_phy_new = real(K_phy_new); 
    
    
    % update R in Fourier space
    force = -1/3 * Hbar * A_phy_old;
    c_k = -100/17;
    R_old = fftshift(fft(R_phy_old));
    force = fftshift(fft(force)); 
    
    R_new = R_old .* exp( - d_k * dt - 1i * k' * c_k * dt) + ( force ./ ( d_k + 1i * k' * c_k) ) .* ( 1 - exp( - d_k * dt - 1i * k' * c_k * dt ) );
    if d_k == 0
        R_new(zero_mode) = R_old(zero_mode) + force(zero_mode) * dt;
    else
        R_new(zero_mode) = R_old(zero_mode) * exp( - d_k * dt ) + force(zero_mode) / d_k * ( 1 - exp( - d_k * dt ) );
    end
    R_phy_new = ifft(ifftshift(R_new));
    R_phy_new = real(R_phy_new); 
    
    % update Z in physical space, where Z is the transformation variable of Q
    noise_Q = max(Eq*3.0,0.01);%*6
    noise_Q = 1*(1-exp(-noise_Q/3));%*2
    noise_Q_all(i) = noise_Q;
    Z_phy_new = Z_phy_old + ( - d_k * Z_phy_old - Hbar * (1 - Qbar) * A_phy_old ) * dt + noise_Q * randn(Na,1)* sqrt(dt);
    % use Z to compute Q
    Q_phy_new = Z_phy_new + Qbar * ( K_phy_new + R_phy_new );
    Q_new = fftshift(fft(Q_phy_new));
    
    % update A in physical space
    A_total = A_phy_old + A_bar(:,i); % total convective activity A
    A_total(A_total <= 0 ) = 1e-5; % positivity
    % the govering equation of A_temp, which is A' (the anomly of A) is
    % given by a linear stochastic equation with multiplicative noise
    A_temp = A_phy_old + (- lambda_A * A_phy_old + Gamma * A_total .* Q_phy_old) * dt + sqrt(0.04*Gamma .* abs(Q_phy_old) .* A_total) .* randn(Na,1) * sqrt(dt);
     
    A_temp(A_temp <= -A_bar(:,i)) = - A_bar(A_temp <= -A_bar(:,i)) + 1e-5; % positivity
    A_phy_new = A_temp; 
    A_new = fftshift(fft(A_phy_new));    
    % put _new variables to _old variables for the next integration step
    % and save these variables
    K_old = K_new; K_phy_old = K_phy_new; 
    R_old = R_new; R_phy_old = R_phy_new; 
    Q_old = Q_new; Q_phy_old = Q_phy_new; 
    A_old = A_new; A_phy_old = A_phy_new; 
    Z_phy_old = Z_phy_new; 
    
    K(:,i) = K_new;
    R(:,i) = R_new;
    Q(:,i) = Q_new;
    A(:,i) = A_new;
    K_phy(:,i) = K_phy_new;
    R_phy(:,i) = R_phy_new;
    Q_phy(:,i) = Q_phy_new;
    A_phy(:,i) = A_phy_new;        

    u_phy(:,i) = ( (K_phy(:,i) - R_phy(:,i)) * psi_0 + 1/sqrt(2) * R_phy(:,i) * psi_2 );
    bd1 = round(140/360*Na);
    bd2 = round(180/360*Na);
    % Governing equation for the MJO: ending here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Wind reconstruction
    u_interannual(i) = u_interannual(i-1) - lambda_A/6 * u_interannual(i-1) * dt + 0.03*Eq * dt + 0.0264*randn*sqrt(dt);
    tau(i) = u_interannual(i) + mean(u_phy(bd1:bd2,i));
end

%averaging the interseasonal wind over the pacific
intra_wind = mean(u_phy(bd1:bd2,:))*dim_u;

%A_bar in fourier space
A_bar_Fourier = real(fftshift(fft(A_bar)));

% Reconstructing the physical variables using meridional bases
u_phy = ( (K_phy - R_phy) * psi_0 + 1/sqrt(2) * R_phy * psi_2 ) * dim_u; % unit 50m/s
theta_phy = ((- K_phy - R_phy) * psi_0 - 1/sqrt(2) * R_phy* psi_2) * dim_theta; % unit 15K 
Q_phy = Q_phy * dim_Q * psi_0; % unit 15K
A_phy = A_phy * dim_A * psi_0; % unit 10K^-1
A_phy_total = A_phy + A_bar * dim_A * psi_0; % unit 10K^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  Postprocessing: Constructing the MJO  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(N,2) == 0
    zero_mode2 = (N)/2+1;
else
    zero_mode2 = (N+1)/2;
end

MJO1_Fourier = zeros(Na, N);
MJO2_Fourier = zeros(Na, N);
MJO1 = zeros(Na, N);
MJO2 = zeros(Na, N);

% MJO is one of the corresponding eigenvectors of the original physical
% system. Therefore, its reconstruction is given by a linear combination of
% K, R, Q and A (all are anomalies) determined by the eigenvector

% modes 1, 2, 3
MJO1_Fourier(zero_mode+1,:) =  K(zero_mode+1,:) * eg_K_MJO(1)' + R(zero_mode+1,:) * eg_R_MJO(1)' + Q(zero_mode+1,:) * eg_Q_MJO(1)' + A(zero_mode+1,:) * eg_A_MJO(1)';
MJO1_Fourier(zero_mode+2,:) =  K(zero_mode+2,:) * eg_K_MJO(2)' + R(zero_mode+2,:) * eg_R_MJO(2)' + Q(zero_mode+2,:) * eg_Q_MJO(2)' + A(zero_mode+2,:) * eg_A_MJO(2)';
MJO1_Fourier(zero_mode+3,:) =  K(zero_mode+3,:) * eg_K_MJO(3)' + R(zero_mode+3,:) * eg_R_MJO(3)' + Q(zero_mode+3,:) * eg_Q_MJO(3)' + A(zero_mode+3,:) * eg_A_MJO(3)';
MJO1_temp = zeros(Na,N);
for i = 1:N
    MJO1_temp(:,i) = ifft(ifftshift(MJO1_Fourier(:,i))); 
end
% modes -1, -2, -3 are the complex conjugates of modes 1, 2, 3
MJO2_Fourier(zero_mode-1,:) =  K(zero_mode-1,:) * eg_K_MJO(1) + R(zero_mode-1,:) * eg_R_MJO(1) + Q(zero_mode-1,:) * eg_Q_MJO(1) + A(zero_mode-1,:) * eg_A_MJO(1);
MJO2_Fourier(zero_mode-2,:) =  K(zero_mode-2,:) * eg_K_MJO(2) + R(zero_mode-2,:) * eg_R_MJO(2) + Q(zero_mode-2,:) * eg_Q_MJO(2) + A(zero_mode-2,:) * eg_A_MJO(2);
MJO2_Fourier(zero_mode-3,:) =  K(zero_mode-3,:) * eg_K_MJO(3) + R(zero_mode-3,:) * eg_R_MJO(3) + Q(zero_mode-3,:) * eg_Q_MJO(3) + A(zero_mode-3,:) * eg_A_MJO(3);
MJO2_temp = zeros(Na,N);
for i = 1:N
    MJO2_temp(:,i) = ifft(ifftshift(MJO2_Fourier(:,i)));
end
% MJO1_temp and MJO2_temp are complex conjugates, which are the spatial
% projections to the leading three modes

% To get the MJO signal, a second projection to the temporal band between
% 30 and 90 days is needed.
% Remember: MJO is the eastward propogating waves.

% Recall the plane wave solution exp(i*k*x - i*omega*t)
% Spatial modes 1, 2 and 3 should correspond to minus 30 to minus 90 days
% due to the minus sign in front of i * omega * t. Be careful for selecting
% the temporal band here. This can equivalently be regarded as considering
% the positive time band 30 to 90 days but put the minus sign into the
% Fourier wavenumber k. That means the modes -1, -2 and -3 should be 
% filtered within the band of 30 to 90 days  
bd11 = -90/(dim_t*dt); bd22 = -30/(dim_t*dt);
for i = 1:Na
    MJO1_temporal = fftshift(fft(MJO1_temp(i,:)));
    MJO1_temporal2 = MJO1_temporal * 0;
    MJO1_temporal2(zero_mode2+round(N/bd22): zero_mode2+round(N/bd11)) = MJO1_temporal(zero_mode2+round(N/bd22): zero_mode2+round(N/bd11));
    MJO1(i,:) = ifft(ifftshift(MJO1_temporal2));
end

bd11 = 30/(dim_t*dt); bd22 = 90/(dim_t*dt);
for i = 1:Na
    MJO2_temporal = fftshift(fft(MJO2_temp(i,:)));
    MJO2_temporal2 = MJO2_temporal * 0;
    MJO2_temporal2(zero_mode2+round(N/bd22): zero_mode2+round(N/bd11)) = MJO2_temporal(zero_mode2+round(N/bd22): zero_mode2+round(N/bd11));
    MJO2(i,:) = ifft(ifftshift(MJO2_temporal2));
end
% The total MJO is the sum of MJO1 and MJO2 which are complex conjugates.
% The sum should be a real-valued field. Yet, due to the roundoff numerical
% error, we put a real(.) function here to elimite those numerical error 
% in the imaginary part, which should be of order 1e-16.
MJO = real(MJO1+MJO2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Postprocessing: Calculating the monthly averaged data  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data points in each month
k_dt = floor(0.5/dt); % Note that one time unit is two months. Therefore, 0.5 time unit is one month.
k_dt_daily = floor(0.5 /30/ dt); % For daily data
% Note that we still keep the temporal high-resolution data of the wind,
% since it is intraseasonal
T_C_og = T_C_3R;
T_E_og = T_E_3R;
h_W_og = h_W_3R;
u_og = u_3R;
I_og = I;

T_C_3R = movmean(T_C_3R,k_dt)'; T_C_3R = T_C_3R(1:k_dt:end);
T_E_3R = movmean(T_E_3R,k_dt)'; T_E_3R = T_E_3R(1:k_dt:end);
h_W_3R = movmean(h_W_3R,k_dt)'; h_W_3R = h_W_3R(1:k_dt:end);
u_3R = movmean(u_3R,k_dt)'; u_3R = u_3R(1:k_dt:end);
I = movmean(I,k_dt)'; I = I(1:k_dt:end);

T_C_3R_daily = movmean(T_C_og,k_dt_daily)'; T_C_3R_daily = T_C_3R_daily(1:k_dt_daily:end);
T_E_3R_daily = movmean(T_E_og,k_dt_daily)'; T_E_3R_daily = T_E_3R_daily(1:k_dt_daily:end);
h_W_3R_daily = movmean(h_W_og,k_dt_daily)'; h_W_3R_daily = h_W_3R_daily(1:k_dt_daily:end);
u_3R_daily = movmean(u_og,k_dt_daily)'; u_3R_daily = u_3R_daily(1:k_dt_daily:end);
I_daily = movmean(I_og,k_dt_daily)'; I_daily = I_daily(1:k_dt_daily:end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Postprocessing: Multiplying the units of state varaibles  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%monthly
T_C_3R = T_C_3R * 7.5;
T_E_3R = T_E_3R * 7.5;
h_W_3R = h_W_3R * 150;
u_3R = u_3R * 1.5;

%daily
T_C_3R_daily = T_C_3R_daily * 7.5;
T_E_3R_daily = T_E_3R_daily * 7.5;
h_W_3R_daily = h_W_3R_daily * 150;
u_3R_daily = u_3R_daily * 1.5;

%wind
tau = tau * dim_u;

% Compute the variance in each month
T_E_3R_seasonal = reshape(T_E_3R,12,[]);
T_E_3R_seasonal = var(T_E_3R_seasonal');
T_C_3R_seasonal = reshape(T_C_3R,12,[]);
T_C_3R_seasonal = var(T_C_3R_seasonal');

% Displaying the Nino 3 and Nino 4 SST statistics from the model
disp('Stats of TC model:')
disp([mean(T_C_3R), var(T_C_3R),skewness(T_C_3R),kurtosis(T_C_3R)]);
disp('Stats of TE model:')
disp([mean(T_E_3R), var(T_E_3R),skewness(T_E_3R),kurtosis(T_E_3R)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Observations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nino 3 observational data
% Load the data
data = importdata('nino3.txt');
nino3 = reshape(data(:,2:end)',[],1);
nino3 = nino3( (1950- 1870) *12 + 1 : end - 12);
% Compute the variance in each month
nino3_seasonal = reshape(nino3,12,[]);
nino3_seasonal = var(nino3_seasonal');

% Nino 4 observational data
% Load the data
data = importdata('nino4.txt');
nino4 = reshape(data(:,2:end)',[],1);
nino4 = nino4( (1950- 1870) *12 + 1 : end - 12);
% Compute the variance in each month
nino4_seasonal = reshape(nino4,12,[]);
nino4_seasonal = var(nino4_seasonal');

% Displaying the Nino 3 and Nino 4 SST statistics from the observations
disp('Stats of Nino4 obs:')
disp([mean(nino4), var(nino4),skewness(nino4),kurtosis(nino4)]);
disp('Stats of Nino3 obs:')
disp([mean(nino3), var(nino3),skewness(nino3),kurtosis(nino3)]);
 

% Displaying the model and observational results for h_W and u as well
disp('Stats of hW model:')
disp([mean(h_W_3R), var(h_W_3R),skewness(h_W_3R),kurtosis(h_W_3R)]);
disp('Stats of u model:')
disp([mean(u_3R), var(u_3R),skewness(u_3R),kurtosis(u_3R)]);

load obs_data
disp('Stats of hW obs:')
disp([mean(h_W_obs), var(h_W_obs),skewness(h_W_obs),kurtosis(h_W_obs)]);
disp('Stats of u obs:')
disp([mean(u_obs), var(u_obs),skewness(u_obs),kurtosis(u_obs)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Comparison the correlations in model and observations  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corr1 = corrcoef(T_E_obs,T_C_obs);
disp(['Corr TE TC obs:   ', num2str(corr1(1,2))]);
corr2 = corrcoef(T_E_3R,T_C_3R);
disp(['Corr TE TC model:   ', num2str(corr2(1,2))]);
corr3 = corrcoef(h_W_obs(1:468),u_obs);
disp(['Corr hW u obs:   ', num2str(corr3(1,2))]);
corr4 = corrcoef(h_W_3R,u_3R);
disp(['Corr hW u model:   ', num2str(corr4(1,2))]);
corr5 = corrcoef(T_E_obs(end-480+1:end),h_W_obs);
disp(['Corr TE hW obs:   ', num2str(corr5(1,2))]);
corr6 = corrcoef(T_E_3R,h_W_3R);
disp(['Corr TE hW model:   ', num2str(corr6(1,2))]);
corr7 =  corrcoef(T_C_obs(end-480+1:end-12),u_obs);
disp(['Corr TC u obs:   ', num2str(corr7(1,2))]);
corr8 = corrcoef(T_C_3R,u_3R);
disp(['Corr TC u model:   ', num2str(corr8(1,2))]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Showing the figures and statistics  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ENSO statistics
ENSO_statistics_results_code




