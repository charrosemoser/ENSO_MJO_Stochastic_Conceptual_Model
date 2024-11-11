%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Model MJO STD by Longitude %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%taking daily MJO results from the model
MJO_new = MJO(:,1:k_dt_daily:end);

%computing the spatiotemporal SST results from the model using the
%bivariate regression coefficients computed in SST_bivar
Hov = zeros(length(lon_sst), length(T_E_3R_daily));
xx = [T_C_3R_daily, T_E_3R_daily];
for i = 1:coarse_x_n
    Hov(i,:) = xx * SST_reg(i,:)'; 
end

window = 30; % look at SST every month
total_loop=floor(length(T_E_3R_daily)/30)-1; % number of months in timeseries

%initializing arrays for counting events
EN_obs=zeros(total_loop,1);
EPEN_obs=zeros(total_loop,1);
CPEN_obs=zeros(total_loop,1);
LN_obs=zeros(total_loop,1);

%take a sample for each month and evalutate to see if an event occurred
for k = 1:total_loop
    if mean(T_E_3R_daily(-0+k*window:+2+k*window)) > 0.5 || mean(T_C_3R_daily(-0+k*window:+2+k*window)) > 0.5
        EN_obs(k)=1;
        if mean(T_E_3R_daily(-0+k*window:+2+k*window))>mean(T_C_3R_daily(-0+k*window:+2+k*window))
            EPEN_obs(k)=1;
        else
            CPEN_obs(k)=1;
        end
    elseif mean(T_E_3R_daily(-0+k*window:+2+k*window)) < -0.5 || mean(T_C_3R_daily(-0+k*window:+2+k*window)) < -0.5
        LN_obs(k)=1;
    end
end

%total model events 
Total_obs = EN_obs | LN_obs;

%model SST spatiotemporal patterns
SST_new = Hov;

% Extract MJO segments for each ENSO type
window_months = 1;
ep_mjo = extract_mjo_segments(MJO_new, EPEN_obs, window_months);
cp_mjo = extract_mjo_segments(MJO_new, CPEN_obs, window_months);
ln_mjo = extract_mjo_segments(MJO_new, LN_obs, window_months);
total_mjo = extract_mjo_segments(MJO_new, Total_obs, window_months);

ep_sst = extract_mjo_segments(SST_new, EPEN_obs, window_months);
cp_sst = extract_mjo_segments(SST_new, CPEN_obs, window_months);
ln_sst = extract_mjo_segments(SST_new, LN_obs, window_months);

% Calculate variance for different longitudes
ep_variance = std(ep_mjo');
cp_variance = std(cp_mjo');
ln_variance = std(ln_mjo');
total_variance = std(total_mjo');



% Plot results
figure;
plot(0:360/Na:360-360/Na, ep_variance, 'r', 'LineWidth', 2);
hold on;
plot(0:360/Na:360-360/Na, cp_variance, 'b', 'LineWidth', 2);
xlabel('Longitude Index');
ylabel('MJO STD');
title('MJO STD for Different ENSO Types from Model Simulation');
legend('EP El Nino', 'CP El Nino');
grid on;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Observed MJO STD by Longitude %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Obs_fine_2020.mat
time2 = 1:13929;% time indices for SST
time3 = 1095+time2; % time indices for MJO

% Pacific ocean domain
Left1_end = (120+2.5) / 2.5 * 10; % Pacific ocean left boundary 120E for SST
Right1_end = 280 / 2.5 * 10; % Pacific ocean right boundary 80W for SST
Middle1_end = 200 / 2.5 * 10; % Middle of Pacific ocean 160W for SST

Left1_C = 160/2.5*10; % CP left boundary index 160W for TC
Right1_C = 210/2.5*10; % CP right boundary index 30E


% Time series
T_E_3R_obs = mean(sst_a_fine(:,Middle1_end:Right1_end),2);% SST averaged over the eastern Pacific
T_C_3R_obs = mean(sst_a_fine(:,Left1_C:Right1_C),2);% SST averaged over the eastern Pacific

window = 30; %take monthly samples to determine ENSO event occurrance
total_loop=floor(length(T_E_3R_obs)/30);%number of months in the observational data set

%initialize observational arrays for counting events
EN_obs=zeros(total_loop,1);
EPEN_obs=zeros(total_loop,1);
CPEN_obs=zeros(total_loop,1);
LN_obs=zeros(total_loop,1);

%count if event occured in a given month
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

%total number of oberved events
Total_obs = EN_obs| LN_obs;

MJO_new = MJO_obs(:,time3);
SST_new = sst_a_fine(time2,Left1_end:Right1_end);
SST_new = SST_new';

% Extract MJO segments for each ENSO type
window_months = 1;
ep_mjo = extract_mjo_segments(MJO_new, EPEN_obs, window_months);
cp_mjo = extract_mjo_segments(MJO_new, CPEN_obs, window_months);
ln_mjo = extract_mjo_segments(MJO_new, LN_obs, -1);% Narrow narrow the range of La Nina selection to avoid overlaps
total_mjo = extract_mjo_segments(MJO_new, Total_obs, window_months);

ep_sst = extract_mjo_segments(SST_new, EPEN_obs, window_months);
cp_sst = extract_mjo_segments(SST_new, CPEN_obs, window_months);
ln_sst = extract_mjo_segments(SST_new, LN_obs, -1);

% Calculate variance for different longitudes
ep_variance = std(ep_mjo');
cp_variance = std(cp_mjo');
ln_variance = std(ln_mjo');
total_variance = std(total_mjo');


% Plot results
figure
plot(0:2.5:360-2.5, ep_variance, 'r-', 'LineWidth', 2);
hold on;
plot(0:2.5:360-2.5, cp_variance, 'b-', 'LineWidth', 2);
xlabel('Longitude Index');
ylabel('MJO STD');
title('MJO STD for Different ENSO Types from Observation');
legend('EP El Nino', 'CP El Nino');
grid on;


%function that extracts the segments of MJO during ENSO events
function mjo_segments = extract_mjo_segments(MJO, event_index, window_months)
    event_starts = find(diff([0; event_index; 0]) == 1);
    event_ends = find(diff([0; event_index; 0]) == -1) - 1;
    
    mjo_segments = [];
    for i = 1:length(event_starts)
        start_month = max(1, event_starts(i) - window_months);
        end_month = min(length(event_index), event_ends(i) + window_months);
        
        start_day = (start_month-1) * 30 + 1;
        end_day = min(end_month * 30, size(MJO, 2));
        
        mjo_segments = [mjo_segments, MJO(:, start_day:end_day)];
    end
end
