%Here we compute the bivariate regression coefficients from the
%observational data. We then plot the coefficients and show the accuraccy
%of the bivariate regression by reconstructing observations and finding the
%residual.

load SST_h_data %contains data for SST we need to do the bivariate regression
load lon_sst %longitudes used when plotting SSTA 
%The data is collected daily. We only take monthly data for T_E and T_C so
%we make the data match. 
[fine_y_n,fine_x_n]= size(sst_a_fine);
k_dt2 = round(365/12); %monthly time steps
Y_coarse = Y(1:k_dt2:end); %monthly time
[dummy,coarse_y_n] = size(Y_coarse);
SST_Y_coarse = zeros(coarse_y_n,fine_x_n);%The array to be populated with the coarser time data and fine location data
for i = 1:fine_x_n
    SST_Y_coarse_temp = movmean(sst_a_fine(:, i),k_dt2);%smooth the fine data across time for given location i
    SST_Y_coarse(:, i) = SST_Y_coarse_temp(1:k_dt2:end); %collect monthly smoothed data from each location
end

% the data is for the entire equitorial band.
%We are only interested in the pacific (longitude 120 to 288 or column 481
%to column 1153). This results in 672 entries which is much more than what
%we need, so we instead only take a subset of these to improve computational efficiency. We will
%collect every 4 values. 
coarse_x_n = 170;%number of locations in the coarse
SST_coarse = zeros(coarse_y_n,coarse_x_n);%final coarse array to update the SST_Y_coarse to also have coarse location data

for i = 1:coarse_y_n
    SST_coarse_temp = movmean(SST_Y_coarse(i, :),4);
    SST_coarse(i,:) = SST_Y_coarse(i,477:4:1153); %collects every 4 locations for each time series 
end

%to keep track of the longitude
X_SST_coarse = X_sst(477:4:1153); %collects every 4 columns 

%do the bivariate regression
SST_reg = zeros(coarse_x_n, 2); %array to be populated with regression coefficients
matrix = [mean(SST_coarse(:,42:92),2), mean(SST_coarse(:,92:152),2)];
for i = 1: coarse_x_n
    b = SST_coarse(:,i);
    coeffs = matrix\b;
    SST_reg(i,:) = transpose(coeffs);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Plotting the Bivariate Regression Coefficients %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
plot(SST_reg(:,1), 'b', LineWidth=2)
plot(SST_reg(:,2), 'r', LineWidth=2)
title('Regression Coefficients', 'FontSize',12)
legend('a(x)','b(x)')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Testing Bivariate Reconstruction on Observations %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CP_obs = mean(SST_coarse(:,42:92),2);%observed SSTA in the CP
EP_obs = mean(SST_coarse(:,92:152),2);%observed SSTA in the EP

Hov = zeros(coarse_x_n, coarse_y_n);
xx = [CP_obs, EP_obs];
% computing the Hovmoller diagram from the bivariate regression on the
% observed CP and EP SSTA
for i = 1:coarse_x_n
    Hov(i,:) = xx * SST_reg(i,:)';
end

%computing the residual 
residual = Hov'- SST_coarse;

%plot the observational data
[xx1,yy1] = meshgrid(X_SST_coarse,Y_coarse);
figure
colormap jet
hold on
subplot(1,3,1)
contourf(xx1,yy1,SST_coarse,40,'linestyle','none')
xlim([120,280])
caxis([-3,3])
title('(b) SSTA Obs')
xlabel('x (longitude)')
ylabel('year')

%plot the bivariate regression reconstructed observational data 
subplot(1,3,2)
contourf(xx1,yy1,Hov',40,'linestyle','none')
xlim([120,280])
caxis([-3,3])
title('(c) Reconstructed')
xlabel('x (longitude)')
yticklabels([])

%plot the residual
subplot(1,3,3)
contourf(xx1,yy1,residual,40,'linestyle','none')
xlim([120,280])
caxis([-3,3])
title('(d) Residual')
xlabel('x (longitude)')
yticklabels([])
