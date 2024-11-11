% Figures and Statistics from the results by running the code main_modeling_code.m
% Note you must run SST_bivar.m to obtain the bivariate regression
% coefficients. 
% The code will generate and show the following results:
% Part I: Statistics
% Part II: Time Series
% Part III: Spatiotemporal Patterns
% Part IV: Counting the event number
% Part V. Bivariate plot of event strength and location

L_realization = round(T/40/12/0.5);%number of 40 year sets in model simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~ Part I: Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SST Spectrums %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spectrums are computed using standard FFT
figure('color','white')
set(gcf,'unit','centimeters','position',[10 5 18 10])
% Spectrum of T_E
subplot(2,4,1)
save_spec = zeros(257, L_realization);
% Find the power spectrum for every set of 40 years
for i = 1:L_realization
    ts = 1*12:1*12:480*12;
    Fs = 1*12;
    index1 = T_E_3R((i-1)*480+1: i*480);
    L = length(ts);
    NFFT = 2^nextpow2(L); 
    Y_y1 = fft(index1,NFFT)/L;
    f_y1 = Fs/2*linspace(0,1,NFFT/2+1);
    tpp_y = 2*abs(Y_y1(1:NFFT/2+1));
    save_spec(:,i) = tpp_y;
end
%mean and std across across the sets of 40 years to plot the confidence
%interval
mean_model= mean(save_spec');
std_model = std(save_spec');
hold on
plot(1./f_y1(end:-1:1) ,mean_model(end:-1:1),'b','linewidth',2) 
low = mean_model-std_model;low(low<1e-10)=1e-10;
x_axis = [f_y1(end:-1:1), fliplr(f_y1(end:-1:1))]; x_axis(x_axis==0)=0.00001;x_axis = 1./x_axis; 
y_axis = [low(end:-1:1) fliplr(mean_model(end:-1:1)+std_model(end:-1:1))];  
facealpha = 0.2;
hh1 = fill(x_axis, y_axis,'b','facealpha',facealpha);
hh1.EdgeColor = 'none';
set(gca,'xscale','log')
xlim([0.8,10]);
set(gca,'fontsize',12)
box on
xlabel('Year')
set(gca,'xTick',[1:6,8,10])
set(gca,'xTicklabel',  [1:6,8,10]);
hold on
%compute the observed power spectrum for EP SSTA
ts = 1*12:1*12:480*12; 
Fs = 1*12; 
L = length(ts);
NFFT = 2^nextpow2(L); 
Y_y1 = fft(nino3,NFFT)/L;
f_y1 = Fs/2*linspace(0,1,NFFT/2+1);
tpp_y = 2*abs(Y_y1(1:NFFT/2+1));
hh2=plot(1./f_y1(end:-1:1),tpp_y(end:-1:1),'r','linewidth',2);
temp_y = 2*abs(Y_y1(1:NFFT/2+1));
tp_y = temp_y; 
title('(a) Spectrum of T_E','fontsize',10)
xlabel('Year')
ylabel('power x frequency')
set(gca,'xtick',[1,2,3,4,5,6,8,10])
set(gca,'xticklabel',[1,2,3,4,5,6,8,10])

% Spectrum of T_C
subplot(2,4,5)
% compute the power spectrum for each set of 40 years
save_spec = zeros(257, L_realization);
for i = 1:L_realization
    ts = 1*12:1*12:480*12;
    Fs = 1*12;
    index1 = T_C_3R((i-1)*480+1: i*480);
    L = length(ts);
    NFFT = 2^nextpow2(L); 
    Y_y1 = fft(index1,NFFT)/L;
    f_y1 = Fs/2*linspace(0,1,NFFT/2+1);
    tpp_y = 2*abs(Y_y1(1:NFFT/2+1));
    save_spec(:,i) = tpp_y;
end
%mean and std across the sets of 40 years to get confidence interval
mean_model= mean(save_spec');
std_model = std(save_spec');
hold on
plot(1./f_y1(end:-1:1) ,mean_model(end:-1:1),'b','linewidth',2) 
low = mean_model-std_model;low(low<1e-10)=1e-10;
x_axis = [f_y1(end:-1:1), fliplr(f_y1(end:-1:1))]; x_axis(x_axis==0)=0.00001;x_axis = 1./x_axis; 
y_axis = [low(end:-1:1) fliplr(mean_model(end:-1:1)+std_model(end:-1:1))];  
facealpha = 0.2;
hh1 = fill(x_axis, y_axis,'b','facealpha',facealpha);
hh1.EdgeColor = 'none';
set(gca,'xscale','log')
xlim([0.8,10]);
set(gca,'fontsize',12)
box on
xlabel('Year')
set(gca,'xTick',[1:6,8,10])
set(gca,'xTicklabel',  [1:6,8,10]);
hold on
%Compute observed power spectrum for CP SSTA
ts = 1*12:1*12:480*12; 
Fs = 1*12; 
L = length(ts);
NFFT = 2^nextpow2(L); 
Y_y1 = fft(nino4,NFFT)/L;
f_y1 = Fs/2*linspace(0,1,NFFT/2+1);
tpp_y = 2*abs(Y_y1(1:NFFT/2+1));
hh2=plot(1./f_y1(end:-1:1),tpp_y(end:-1:1),'r','linewidth',2);
temp_y = 2*abs(Y_y1(1:NFFT/2+1));
tp_y = temp_y; 
title('(a) Spectrum of T_C','fontsize',10)
xlabel('Year')
ylabel('power x frequency')
set(gca,'xtick',[1,2,3,4,5,6,8,10])
set(gca,'xticklabel',[1,2,3,4,5,6,8,10])

 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SST PDFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PDFs of TE and TC from model and observation
[PDF_T_E_obs, xx_T_E] = ksdensity(nino3);
[PDF_T_E_model, xx_T_E] = ksdensity(T_E_3R,xx_T_E);
[PDF_T_C_obs, xx_T_C] = ksdensity(nino4);
[PDF_T_C_model, xx_T_C] = ksdensity(T_C_3R,xx_T_C);

%calculating the PDF for each set of 40 years in the simulation to create
%confidence interval
PDF_T_E_model_UQ = zeros(L_realization,100);
PDF_T_C_model_UQ = zeros(L_realization,100);
for i = 1:L_realization
    temp = ksdensity(T_E_3R((i-1)*480+1: i*480),xx_T_E);
    PDF_T_E_model_UQ(i,:) = temp;   
    temp = ksdensity(T_C_3R((i-1)*480+1: i*480),xx_T_C);
    PDF_T_C_model_UQ(i,:) = temp;   
end
% plot the Nino 3 PDFs
subplot(2,4,2)
hold on
plot(xx_T_E,PDF_T_E_obs,'r','linewidth',2)
plot(xx_T_E,PDF_T_E_model,'b','linewidth',2)
upper_lim = mean(PDF_T_E_model_UQ) + 2 * std(PDF_T_E_model_UQ);
lower_lim = mean(PDF_T_E_model_UQ) - 2 * std(PDF_T_E_model_UQ); lower_lim(lower_lim<=0) = 0;
patch([xx_T_E, xx_T_E(end:-1:1)],[lower_lim,upper_lim(end:-1:1)],'b','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',9)
title('(c) PDF of T_E','fontsize',10)
set(gca,'xtick',-2:2:4);
set(gca,'ytick',0:0.1:0.5);
% plot the Nino 4 PDFs
subplot(2,4,6)
hold on
plot(xx_T_C,PDF_T_C_obs,'r','linewidth',2)
plot(xx_T_C,PDF_T_C_model,'b','linewidth',2)
upper_lim = mean(PDF_T_C_model_UQ) + 2 * std(PDF_T_C_model_UQ);
lower_lim = mean(PDF_T_C_model_UQ) - 2 * std(PDF_T_C_model_UQ); lower_lim(lower_lim<=0) = 0;
patch([xx_T_C, xx_T_C(end:-1:1)],[lower_lim,upper_lim(end:-1:1)],'b','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',9)
title('(d) PDF of T_C','fontsize',10)
xlabel('¡ãC')
set(gca,'xtick',-2:2:4);
set(gca,'ytick',0:0.2:0.6);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PDFs of h_W and u %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load obs_data % observational data for h_W and u
%calculate PDFs for h_W and u from model and observations
[PDF_h_W_obs, xx_h_W] = ksdensity(h_W_obs);
[PDF_h_W_model, xx_h_W] = ksdensity(h_W_3R,xx_h_W);
[PDF_u_obs, xx_u] = ksdensity(u_obs);
[PDF_u_model, xx_u] = ksdensity(u_3R,xx_u);

%calculate the PDFs of each set of 40 years of model simulation to get
%confidence interval
PDF_h_W_model_UQ = zeros(L_realization,100);
PDF_u_model_UQ = zeros(L_realization,100);
for i = 1:L_realization
    temp = ksdensity(h_W_3R((i-1)*480+1: i*480),xx_h_W);
    PDF_h_W_model_UQ(i,:) = temp;   
    temp = ksdensity(u_3R((i-1)*480+1: i*480),xx_u);
    PDF_u_model_UQ(i,:) = temp;   
end

% plot the PDFs of h_W
subplot(2,4,4)
hold on
plot(xx_h_W,PDF_h_W_obs,'r','linewidth',2)
plot(xx_h_W,PDF_h_W_model,'b','linewidth',2)
upper_lim = mean(PDF_h_W_model_UQ) + 2 * std(PDF_h_W_model_UQ);
lower_lim = mean(PDF_h_W_model_UQ) - 2 * std(PDF_h_W_model_UQ); lower_lim(lower_lim<=0) = 0;
patch([xx_h_W, xx_h_W(end:-1:1)],[lower_lim,upper_lim(end:-1:1)],'b','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',9)
title('(g) PDF of h_W','fontsize',10)
xlabel('m')
set(gca,'xtick',-40:40:40);
set(gca,'ytick',0:0.01:0.03);

% plot the PDFs of u
subplot(2,4,8)
hold on
plot(xx_u,PDF_u_obs,'r','linewidth',2)
plot(xx_u,PDF_u_model,'b','linewidth',2)
upper_lim = mean(PDF_u_model_UQ) + 2 * std(PDF_u_model_UQ);
lower_lim = mean(PDF_u_model_UQ) - 2 * std(PDF_u_model_UQ); lower_lim(lower_lim<=0) = 0;
patch([xx_u, xx_u(end:-1:1)],[lower_lim,upper_lim(end:-1:1)],'b','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',9)
title('(h) PDF of u','fontsize',10)
xlabel('m/s')
set(gca,'xtick',-0.4:0.4:0.4);
set(gca,'ytick',0:3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% SST seasonal variations %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find the seasonal variation from each set of 40 years to compute 
% confidence interval
T_E_3R_seasonal_UQ = zeros(L_realization,12);
T_C_3R_seasonal_UQ = zeros(L_realization,12);
for i = 1:L_realization
    temp = reshape(T_E_3R((i-1)*480+1: i*480),12,[]);
    T_E_3R_seasonal_UQ(i,:) = var(temp');   
    temp = reshape(T_C_3R((i-1)*480+1: i*480),12,[]);
    T_C_3R_seasonal_UQ(i,:) = var(temp');    
end
% plot the seasonal variation of Nino 3 SSTs
subplot(2,4,3)
hold on
plot(1:12,nino3_seasonal,'r','linewidth',2)
plot(1:12,T_E_3R_seasonal,'b','linewidth',2)
upper_lim = mean(T_E_3R_seasonal_UQ) + 2 * std(T_E_3R_seasonal_UQ);
lower_lim = mean(T_E_3R_seasonal_UQ) - 2 * std(T_E_3R_seasonal_UQ); lower_lim(lower_lim<=0) = 0;
patch([1:12, 12:-1:1],[lower_lim,upper_lim(end:-1:1)],'b','facealpha',0.2,'linestyle','none')
box on
set(gca,'fontsize',9)
title('(e) Variance of T_E','fontsize',10)
xlabel('Calendar month')
set(gca,'xtick',1:2:12,'xticklabel',{'J','M','M','J','S','N'})
set(gca,'xlim',[1,12],'ylim',[0 2]);
set(gca,'ytick',0:0.5:2);
% plot the seasonal variation of Nino 4 SSTs
subplot(2,4,7)
hold on
plot(1:12,nino4_seasonal,'r','linewidth',2)
plot(1:12,T_C_3R_seasonal,'b','linewidth',2)
upper_lim = mean(T_C_3R_seasonal_UQ) + 2 * std(T_C_3R_seasonal_UQ);
lower_lim = mean(T_C_3R_seasonal_UQ) - 2 * std(T_C_3R_seasonal_UQ); lower_lim(lower_lim<=0) = 0;
patch([1:12, 12:-1:1],[lower_lim,upper_lim(end:-1:1)],'b','facealpha',0.2,'linestyle','none')
box on
set(gca,'fontsize',9)
title('(f) Variance of T_C','fontsize',10)
xlabel('Calendar month')
set(gca,'xlim',[1,12],'ylim',[0 1]);
set(gca,'xtick',1:2:12,'xticklabel',{'J','M','M','J','S','N'})
set(gca,'ytick',0:0.2:1);

 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% ACFs of SSTs, h_W and u %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find the autocorrelation function from observation and model simulation
ACF_T_E_model = autocorr(T_E_3R,60);
ACF_T_C_model = autocorr(T_C_3R,60);
ACF_T_E_obs = autocorr(nino3,60);
ACF_T_C_obs = autocorr(nino4,60);
ACF_u_model = autocorr(u_3R,60);
ACF_h_W_model = autocorr(h_W_3R,60);
ACF_u_obs = autocorr(u_obs,60);
ACF_h_W_obs = autocorr(h_W_obs,60);

%calculate the autocorrelation function for each set of 40 years from the
%model simulation to construct the confidence interval
ACF_T_E_model_UQ = zeros(L_realization,61);
ACF_T_C_model_UQ = zeros(L_realization,61);
ACF_h_W_model_UQ = zeros(L_realization,61);
ACF_u_model_UQ = zeros(L_realization,61);
for i = 1:L_realization
    temp = autocorr(T_E_3R((i-1)*480+1: i*480),60);
    ACF_T_E_model_UQ(i,:) = temp;   
    temp = autocorr(T_C_3R((i-1)*480+1: i*480),60);
    ACF_T_C_model_UQ(i,:) = temp;
    temp = autocorr(h_W_3R((i-1)*480+1: i*480),60);
    ACF_h_W_model_UQ(i,:) = temp;   
    temp = autocorr(u_3R((i-1)*480+1: i*480),60);
    ACF_u_model_UQ(i,:) = temp;
end

figure
%plot ACFs of T_E
subplot(2,2,1)
hold on
plot([0:60]/12, ACF_T_E_obs,'r','linewidth',2)
plot([0:60]/12, ACF_T_E_model,'b','linewidth',2)
upper_lim = mean(ACF_T_E_model_UQ) + 2 * std(ACF_T_E_model_UQ);
lower_lim = mean(ACF_T_E_model_UQ) - 2 * std(ACF_T_E_model_UQ); 
patch([[0:60]/12, [60:-1:0]/12],[lower_lim,upper_lim(end:-1:1)],'b','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',9)
legend('Nino 3','T_E')
title('ACF','fontsize',9)
%plot ACFs of T_C
subplot(2,2,3)
hold on
plot([0:60]/12, ACF_T_C_obs,'r','linewidth',2)
plot([0:60]/12, ACF_T_C_model,'b','linewidth',2)
upper_lim = mean(ACF_T_C_model_UQ) + 2 * std(ACF_T_C_model_UQ);
lower_lim = mean(ACF_T_C_model_UQ) - 2 * std(ACF_T_C_model_UQ); 
patch([[0:60]/12, [60:-1:0]/12],[lower_lim,upper_lim(end:-1:1)],'b','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',9)
legend('Nino 4','T_C')
%plot ACFs of h_W
subplot(2,2,2)
hold on
plot(0:60,ACF_h_W_obs,'r','linewidth',2)
plot(0:60,ACF_h_W_model,'b','linewidth',2)
upper_lim = mean(ACF_h_W_model_UQ) + 2 * std(ACF_h_W_model_UQ);
lower_lim = mean(ACF_h_W_model_UQ) - 2 * std(ACF_h_W_model_UQ); 
patch([[0:60], [60:-1:0]],[lower_lim,upper_lim(end:-1:1)],'b','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',9)
legend('hW Obs','hW model')
title('ACF h_W')
%plot ACFs of u
subplot(2,2,4)
hold on
plot(0:60,ACF_u_obs,'r','linewidth',2)
plot(0:60,ACF_u_model,'b','linewidth',2)
upper_lim = mean(ACF_u_model_UQ) + 2 * std(ACF_u_model_UQ);
lower_lim = mean(ACF_u_model_UQ) - 2 * std(ACF_u_model_UQ); 
patch([[0:60], [60:-1:0]],[lower_lim,upper_lim(end:-1:1)],'b','facealpha',0.2,'linestyle','none');
box on
set(gca,'fontsize',9)
legend('u Obs','u model')
title('ACF u')


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~ Part II: Time Series ~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Comparison of the model and observational time series %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h_W_obs_new = h_W_obs(1:end-12);
u_obs_new = u_obs;
T_E_obs_new = T_E_obs(30*12+1:end-12);
T_C_obs_new = T_C_obs(30*12+1:end-12);
t_obs = 1980+1/12: 1/12: 2019;
LL = length(t_obs); % let the length of the time series form the model be the same as observations 
%range_model = 1+12*120: LL+12*120;% time snippet for Figure 1
range_model = 1+12*490: LL+12*490; % time snippet for Figure S4
range_tau = range_model(1)*k_dt:range_model(end)*k_dt;
t_model = range_model/12; % time window of plotting variables except the wind 
tau_model = range_tau/k_dt/12; % the time window is more refined for the wind as we did not apply the monthly average

% When plotting the time series, two different y-axes may be used
figure('color','white')
set(gcf,'unit','centimeters','position',[10 5 18 15])
% Nino 3 and 4 SSTs, observations
subplot(5,1,1)
hold on
yyaxis left
plot(t_obs,T_E_obs_new,'r','linewidth',2)
set(gca,'ycolor','r')
xlim([t_obs(1),t_obs(end)])
ylim([-3.5,3.5])
ylabel('^oC')
yyaxis right
plot(t_obs,T_C_obs_new,'g','linewidth',2)
set(gca,'ycolor','g')
ylim([-3.5,3.5])
h = legend('T_E','T_C','orientation','horizontal','location','southwest');
set(h,'box','off');
xlim([t_obs(1),t_obs(end)])
set(gca,'fontsize',9)
ylabel('^oC')
box on
text(1973.5,0,'(a) Obs','fontsize',9)
grid on
grid(gca,'minor')
title('Comparison of the observational time series and model simulations','fontsize',12)
% h_W and u, observations
subplot(5,1,2)
hold on
yyaxis left
plot(t_obs,h_W_obs_new,'b','linewidth',2)
set(gca,'ycolor','b')
xlim([t_obs(1),t_obs(end)])
ylim([-50,50])
ylabel('m')
yyaxis right
plot(t_obs,u_obs_new,'k','linewidth',2)
set(gca,'ycolor','k')
ylim([-0.5,0.5])
ylabel('m/s')
h=legend('h_W','u','orientation','horizontal','location','northwest');
set(h,'box','off');
xlim([t_obs(1),t_obs(end)])
set(gca,'fontsize',9)
box on
text(1973.5,0,'(b) Obs','fontsize',9)
grid on
grid(gca,'minor')
% Nino 3 and 4 SSTs, model
subplot(5,1,3)
hold on
yyaxis left
plot(t_model,T_E_3R(range_model),'r','linewidth',2)
set(gca,'ycolor','r')
xlim([t_model(1),t_model(end)])
ylim([-3.5,3.5])
ylabel('^oC')
yyaxis right
plot(t_model,T_C_3R(range_model),'g','linewidth',2)
set(gca,'ycolor','g')
ylim([-3.5,3.5])
ylabel('^oC')
h=legend('T_E','T_C','orientation','horizontal','location','southwest');
set(h,'box','off');
xlim([t_model(1),t_model(end)])
set(gca,'fontsize',9)
box on
text(t_model(1)-6.5,0,'(c) Model','fontsize',9)
grid on
grid(gca,'minor')
% h_W and u, model
subplot(5,1,4)
hold on
yyaxis left
plot(t_model,h_W_3R(range_model),'b','linewidth',2)
set(gca,'ycolor','b')
xlim([t_model(1),t_model(end)])
ylim([-50,50])
ylabel('m')
yyaxis right
plot(t_model,u_3R(range_model),'k','linewidth',2)
set(gca,'ycolor','k')
ylim([-0.5,0.5])
ylabel('m/s')
h=legend('h_W','u','orientation','horizontal','location','northwest');
set(h,'box','off');
xlim([t_model(1),t_model(end)])
set(gca,'fontsize',9)
box on
text(t_model(1)-6.5,0,'(d) Model','fontsize',9)
grid on
grid(gca,'minor')
% I and tau from model 
subplot(5,1,5)
hold on
yyaxis left
plot(tau_model,tau(range_tau),'c','linewidth',2)
xlim([t_model(1),t_model(end)])
ylim([-25,25])
set(gca,'ycolor','c')
ylabel('m/s')
yyaxis right
plot(t_model,I(range_model),'m','linewidth',2)
xlim([t_model(1),t_model(end)])
set(gca,'ycolor','m')
ylim([0,4])
h = legend('\tau','I','orientation','horizontal','location','southwest');
set(h,'box','off');
box on
set(gca,'fontsize',9)
text(t_model(1)-6.5,2,'(e) Model','fontsize',9)
grid on
grid(gca,'minor')
xlabel('t')


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~ Part III: Spatiotemporal Patterns ~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Spatiotemporal reconstruction using the Nino 3 and 4 indices %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find the spatio-temporal reconstruction from the model simulation using
%the bivariate regression coefficients computed in SST_bivar.m
Hov = zeros(length(lon_sst), length(T_E_3R));
xx = [T_C_3R, T_E_3R];
for i = 1:coarse_x_n
    Hov(i,:) = xx * SST_reg(i,:)'; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Plotting the SST spatiotemporal patterns %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Length of period for each subpanel in showing the spatiotemporal pattern
k_dt = 0.5 / dt;
year_add = 20;
t_temp = 1980+1/12: 1/12: 1980+year_add;  
LL = length(t_temp);

window = 12;
total_loop = (year_add - 1) * 12/window;
figure('color','white')
set(gcf,'unit','centimeters','position',[10 5 18 15])
colormap(jet)
% Plot 5 different periods
for j = 1:5

    if j == 1
        range_model = 1+12*60:LL+12*60;
    elseif j == 2
        range_model = 1+12*80:LL+12*80;
    elseif j == 3
        range_model = 1+12*100:LL+12*100;
    elseif j == 4
        range_model = 1+12*120:LL+12*120;
    elseif j == 5
        range_model = 1+12*140:LL+12*140;
    end
    range_tau = range_model(1)*k_dt:range_model(end)*k_dt; % the wind tau has a different temporal resolution and needs to be handled separately
    t_model = range_model/12; % time window for variables except tau
    tau_model = range_tau/k_dt/12; % time wind for tau
    Hov = zeros(coarse_x_n,length(range_model)); % the matrix is for the spatiotemporal pattern, namely the Hovmoller diagram
    
    xx = [T_C_3R(range_model),T_E_3R(range_model)]; % range of the longitude across equatorial Pacific
    % computing the Hovmoller diagram from the bivariate regression
    for i = 1:coarse_x_n
        Hov(i,:) = xx * SST_reg(i,:)';
    end
    
    % plot the Hovmoller diagram of SST
    subplot(1,5,j)
    [xx,yy] = meshgrid(t_model,X_SST_coarse);
    contourf(yy,xx,Hov,30,'linestyle','none')
    hold on
    plot([180 180],[t_model(1) t_model(end)],'m--','linewidth',2);
    temp_tau = range_tau(1:k_dt/5:end); % plot tau on top of SST Hovmoller diagram; zero-wind is at the dateline; negative value means easterly and positive value means westerly
    plot(180+tau(temp_tau)*20,tau_model(1:k_dt/5:end),'k','linewidth',0.5);
    % Use different colored labels to point out different type of ENSO
    % events, including Strong EP El Nino, moderate EP El Nino, CP El Nino, and La Nina
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
    xlabel('longitude');
    set(gca,'fontsize',9,'linewidth',1);
    caxis([-3,3])
    if j == 1
        ylabel('model year');
    end
end
sgtitle('Hovmoller diagrams of the standard run','fontsize',12);
colorbar('eastoutside','position',[.92,0.12,0.01,0.8],'yTick',-4:1:4,'fontsize',9,'fontname','arial');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~ Part IV: Counting the event number ~~~~~~~~~~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Constructing the SST spatiotemporal patterns for a simulation %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%produce spatio temporal results from the bivariate regression with the
%model results plugged in 
range_model = 1:length(T_E_3R);
Hov = zeros(length(lon_sst), length(T_E_3R));
xx = [T_C_3R, T_E_3R];
for i = 1:coarse_x_n
    Hov(i,:) = xx * SST_reg(i,:)'; 
end

total_loop = length(T_E_3R)/12-1;
sst_max_long = zeros(total_loop,2);
sst_min_long_LN = zeros(total_loop,2);

EN = zeros(total_loop,1); % total number of El Nino events 
EPEN = zeros(total_loop,1); % total number of EP El Nino events
CPEN = zeros(total_loop,1); % total number of CP El Nino events
LN = zeros(total_loop,1); % total number of La Nina events
EEN = zeros(total_loop,1);% total number of extreme El Nino events. Nino 3 SST > 2.5 degree Celsius.

%populating a one in each array if the respective event occured
for k = 1:total_loop
    if max(T_E_3R(-8+k*window:+3+k*window)) >= 2.5
        EEN(k) = 1;
    end
    if mean(T_E_3R(-0+k*window:+2+k*window)) > 0.5 || mean(T_C_3R(-0+k*window:+2+k*window)) > 0.5
        EN(k) = 1;
        Hov_cal = mean(Hov(:,-0+k*window:+2+k*window),2);
        sst_max_long(k,1) = max(Hov_cal); % value
        sst_max_long(k,2) = find(Hov_cal == max(Hov_cal)); % longitude
        if mean(T_E_3R(-0+k*window:+2+k*window)) > mean(T_C_3R(-0+k*window:+2+k*window))
            EPEN(k) = 1;
        else
            CPEN(k) = 1;
        end
    elseif mean(T_E_3R(-0+k*window:+2+k*window)) < -0.5 || mean(T_C_3R(-0+k*window:+2+k*window)) < -0.5
        LN(k) = 1;  
        Hov_cal = mean(Hov(:,-0+k*window:+2+k*window),2);
        sst_min_long_LN(k,1) = min(Hov_cal); % value
        sst_min_long_LN(k,2) = find(Hov_cal == min(Hov_cal)); % longitude
    end
end
MYEN = EN(1:end-1) .* EN(2:end); % total number of multi-year El Nino events
MYLN = LN(1:end-1) .* LN(2:end); % total number of multi-year La Nina events
for i = 1:length(MYEN)-1
    if MYEN(i) == 1 && MYEN(i+1) == 1
        MYEN(i+1) = 0;
    end
    if MYLN(i) == 1 && MYLN(i+1) == 1
        MYLN(i+1) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Observed Peak SSTA Spatiotemporal Results %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obs_total_loop = floor(length(Y_coarse)/12)-1;
sst_max_long_obs = zeros(obs_total_loop,2);

%populating a one in each array if the respective event occured
for k = 1:obs_total_loop
    if mean(CP_obs(-0+k*window:+2+k*window)) > 0.5 || mean(EP_obs(-0+k*window:+2+k*window)) > 0.5
        Hov_cal = mean(SST_coarse(-0+k*window:+2+k*window,:),1);
        sst_max_long_obs(k,1) = max(Hov_cal); % value
        sst_max_long_obs(k,2) = find(Hov_cal == max(Hov_cal)); % longitude

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Displaying the number of each event %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

disp('Total number of each time of events every 70 years')
disp('El Nino,   EP,    CP,    Extreme,    Multi-year El Nino,    La Nina,    Multi-year La Nina')
disp('Model result')
disp([sum(EN)/(T/6/70), sum(EPEN)/(T/6/70), sum(CPEN)/(T/6/70), sum(EEN)/(T/6/70), sum(MYEN)/(T/6/70), sum(LN)/(T/6/70), sum(MYLN)/(T/6/70)])
disp('Observations')
disp([24.00001, 14.00001, 10.00001, 4.00001, 5.00001, 24.00001, 8.00001])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Bar Chart of the number of each event %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find standard deviation of each event
%find number of sets of 70 years in the model simulation
n_tseg = floor(length(EN)/70);

%initializing the variables 
EEN_UQ = zeros(n_tseg,1);
EN_UQ = zeros(n_tseg,1);
EPEN_UQ = zeros(n_tseg,1);
CPEN_UQ = zeros(n_tseg,1);
LN_UQ = zeros(n_tseg,1);
MYEN_UQ = zeros(n_tseg,1);
MYLN_UQ = zeros(n_tseg,1);

%populating each element with the number of events in that set of 70 years
for i = 1:n_tseg
    EEN_UQ(i) = sum(EEN(70*(i-1)+1:70*i));
    EN_UQ(i) = sum(EN(70*(i-1)+1:70*i));
    EPEN_UQ(i) = sum(EPEN(70*(i-1)+1:70*i));
    CPEN_UQ(i) = sum(CPEN(70*(i-1)+1:70*i));
    LN_UQ(i) = sum(LN(70*(i-1)+1:70*i));
    MYEN_UQ(i) = sum(MYEN(70*(i-1)+1:70*i));
    MYLN_UQ(i) = sum(MYLN(70*(i-1)+1:70*i));
end


%calculating the lower and upper limits of the error bars by the std of the
%respective event count
EEN_upper_lim = mean(EEN_UQ) + 2 * std(EEN_UQ);
EEN_lower_lim = mean(EEN_UQ) - 2 * std(EEN_UQ);EEN_lower_lim(EEN_lower_lim<=0) = 0;%positivity

EN_upper_lim = mean(EN_UQ) + 2 * std(EN_UQ);
EN_lower_lim = mean(EN_UQ) - 2 * std(EN_UQ); EN_lower_lim(EN_lower_lim<=0) = 0;

EPEN_upper_lim = mean(EPEN_UQ) + 2 * std(EPEN_UQ);
EPEN_lower_lim = mean(EPEN_UQ) - 2 * std(EPEN_UQ); EPEN_lower_lim(EPEN_lower_lim<=0) = 0;


CPEN_upper_lim = mean(CPEN_UQ) + 2 * std(CPEN_UQ);
CPEN_lower_lim = mean(CPEN_UQ) - 2 * std(CPEN_UQ); CPEN_lower_lim(CPEN_lower_lim<=0) = 0;

LN_upper_lim = mean(LN_UQ) + 2 * std(LN_UQ);
LN_lower_lim = mean(LN_UQ) - 2 * std(LN_UQ); LN_lower_lim(LN_lower_lim<=0) = 0;

MYEN_upper_lim = mean(MYEN_UQ) + 2 * std(MYEN_UQ);
MYEN_lower_lim = mean(MYEN_UQ) - 2 * std(MYEN_UQ); MYEN_lower_lim(MYEN_lower_lim<=0) = 0;


MYLN_upper_lim = mean(MYLN_UQ) + 2 * std(MYLN_UQ);
MYLN_lower_lim = mean(MYLN_UQ) - 2 * std(MYLN_UQ); MYLN_lower_lim(MYLN_lower_lim<=0) = 0;


% creating the array to plot and the upper and lower limits for the error
% bars
counts = [sum(EN)/(T/6/70), 24; sum(EPEN)/(T/6/70), 14; sum(CPEN)/(T/6/70), 10; sum(EEN)/(T/6/70), 4; sum(MYEN)/(T/6/70), 5; sum(LN)/(T/6/70), 24; sum(MYLN)/(T/6/70), 8];
upper_lims = [EN_upper_lim, NaN; EPEN_upper_lim, NaN; CPEN_upper_lim, NaN; EEN_upper_lim, NaN; MYEN_upper_lim, NaN; LN_upper_lim, NaN; MYLN_upper_lim, NaN];
lower_lims = [EN_lower_lim, NaN; EPEN_lower_lim, NaN; CPEN_lower_lim, NaN; EEN_lower_lim, NaN; MYEN_lower_lim, NaN; LN_lower_lim, NaN; MYLN_lower_lim, NaN];

figure
bar(counts, 'grouped')
title('ENSO Event Frequencies Over a 70 year Period')
xlabel('Type of Event')
ylabel('Number of events per 70 year period')
set(gca, 'XTickLabel', {'El Nino' 'EP El Nino' 'CP El Nino' 'Extreme El Nino' 'Multi-Year El Nino' 'La Nina' 'Multi-Year La Nina'});
hold on

% Add error bars only to the model bars
numGroups = size(counts, 1);
numBars = size(counts, 2);

% Calculate the width of each group of bars
groupwidth = min(0.8, numBars/(numBars + 1.5));

for i = 1:numBars
    % Calculate the center of each bar
    b = (1:numGroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numBars); 
    er = errorbar(b, counts(:,i), counts(:,i) - lower_lims(:,i), upper_lims(:,i) - counts(:,i), 'k', 'linestyle', 'none');
end

hold off


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%~~~~~~~~~ Part V. Bivariate plot of event strength and location ~~~~~~~~%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
load obs_peak %observed SST peaks

% locations of the El Nino event peak each year (if it's an El Nino year)
%model
sst_max_long_num = zeros(length(lon_sst),1);
for i = 1:length(lon_sst)
    sst_max_long_num(i) = sum(sst_max_long(:,2) == i);
end

% Strength of the event
%model
sst_max_num = zeros(50,1);
for i = 1:50
    for k = 1:total_loop
        if sst_max_long(k,1) > (i-1) * 0.1 && sst_max_long(k,1) <= i * 0.1 % find the value of the maximum
            sst_max_num(i) = sst_max_num(i)+1;
        end
    end
end

% locations of the La Nina event peak each year (if it's a La Nina year)
sst_max_long_num_LN = zeros(length(lon_sst),1);
for i = 1:length(lon_sst)
    sst_max_long_num_LN(i) = sum(sst_min_long_LN(:,2) == i);
end
% Strength of the event
sst_min_num_LN = zeros(50,1);
for i = 1:50
    for k = 1:total_loop
        if sst_min_long_LN(k,1) < (i-1) * (-0.1) && sst_min_long_LN(k,1) >= i*(-0.1) % find the value of the minimum
            sst_min_num_LN(i) = sst_min_num_LN(i) + 1;
        end
    end
end

% El Nino events
figure('color','white')
set(gcf,'unit','centimeters','position',[20 5 18 15])
sgtitle({'Bivariate distribution of DJF El Nino SSTA peaks';'2000yr Ctrl, averaged 5S-5N'},'fontsize',14);

% distribution of the event peak 
subplot(3,3,1:2)
plot(lon_sst,sst_max_long_num,'k','linewidth',2);
hold on
plot([160 160],[0 100],'k--','linewidth',0.5);
plot([210 210],[0 100],'k--','linewidth',0.5);
plot([270 270],[0 100],'k--','linewidth',0.5);
set(gca,'xlim',[120 280],'ylim',[0 100]);
set(gca,'xtick',120:20:280,'xticklabel',[]);
set(gca,'ytick',0:20:100);
ylabel('event count');
set(gca,'fontsize',14);

% bivariate plot
subplot(3,3,[4:5 7:8])
hold on
scatter(lon_sst(sst_max_long(find(sst_max_long(:,2)~=0),2)),sst_max_long(find(sst_max_long(:,2)~=0),1),'b*');
scatter(peak_long,peak_temp,'r', 'filled');
plot([160 160],[0 6],'k--','linewidth',0.5);
plot([210 210],[0 6],'k--','linewidth',0.5);
plot([270 270],[0 6],'k--','linewidth',0.5);
set(gca,'xlim',[120 280],'ylim',[0 5]);
set(gca,'xtick',120:20:280);
set(gca,'ytick',0:1:5);
ylabel('peak SSTA (degC)');
xlabel('longitude of peak SSTA');
box on
set(gca,'fontsize',14);

% distribution of the strength
subplot(3,3,[6 9])
plot(sst_max_num,0:0.1:4.9,'k','linewidth',2);
set(gca,'xlim',[0 100],'ylim',[0 5]);
set(gca,'xtick',0:20:100);
set(gca,'ytick',0:1:5);
xlabel('event count');
set(gca,'fontsize',14);
title([num2str(sum(EN)) ' total warm events'],'fontsize',10);

% La Nina events
figure('color','white')
set(gcf,'unit','centimeters','position',[20 5 18 15])
sgtitle({'Bivariate distribution of DJF La Nina SSTA peaks';'2000yr Ctrl, averaged 5S-5N'},'fontsize',14);

% distribution of the event peak 
subplot(3,3,1:2)
plot(lon_sst,sst_max_long_num_LN,'k','linewidth',2);
hold on
plot([160 160],[0 200],'k--','linewidth',0.5);
plot([210 210],[0 200],'k--','linewidth',0.5);
plot([270 270],[0 200],'k--','linewidth',0.5);
set(gca,'xlim',[120 280],'ylim',[0 200]);
set(gca,'xtick',120:20:280,'xticklabel',[]);
set(gca,'ytick',0:50:200);
ylabel('event count');
set(gca,'fontsize',14);

% bivariate plot
subplot(3,3,[4:5 7:8])
scatter(lon_sst(sst_min_long_LN(find(sst_min_long_LN(:,2)~=0),2)),sst_min_long_LN(find(sst_min_long_LN(:,2)~=0),1),'k*');
hold on
plot([160 160],[-5 0],'k--','linewidth',0.5);
plot([210 210],[-5 0],'k--','linewidth',0.5);
plot([270 270],[-5 0],'k--','linewidth',0.5);
set(gca,'xlim',[120 280],'ylim',[-5 0]);
set(gca,'xtick',120:20:280);
set(gca,'ytick',-5:1:0);
ylabel('peak SSTA (degC)');
xlabel('longitude of peak SSTA');
box on
set(gca,'fontsize',14);

% distribution of the strength
subplot(3,3,[6 9])
plot(sst_min_num_LN,0:-0.1:-4.9,'k','linewidth',2);
set(gca,'xlim',[0 100],'ylim',[-5 0]);
set(gca,'xtick',0:20:100);
set(gca,'ytick',-5:1:0);
xlabel('event count');
set(gca,'fontsize',14);
title([num2str(sum(LN)) ' total cold events'],'fontsize',10);



