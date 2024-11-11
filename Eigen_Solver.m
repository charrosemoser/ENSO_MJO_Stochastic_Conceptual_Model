% Computing the eigenvalues and eigenvectors of the linearized skeleton model

% The original MJO Skeleton model had a time unit be 34 days.
% Rescale the time unit to make one unit be 2 months, consistent with the
% ENSO model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Defining the parameters %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Skeleton model, parameters can be found in Thual et al JAS paper
% Whenever a parameter appears on the right-hand side of the governing
Q_tilde = 0.9; % Q_bar parameter in the paper. I try to avoid using bar for parameters since bar indicates the mean state
Gamma = 4980/17; % Capital Gamma parameter
epsilon = 0.0; % coefficient of the additional damping, which should be small
d_bar = 1; % additional damping if needed
H_bar = 660/17; % H_bar parameter
A_bar = .1; % the mean state of variable A
L = 8/3; % non-dimensional form of the total length of the equatorial band
dim_t = 60; % unit of time, which is 2 months
Eig_Store = zeros(4,11); % the four eigenvalues for the first 11 wavenumbers
Ev1 = zeros(4,11); % amplitude of the first component of the eigenvalue 
Ev2 = zeros(4,11); % amplitude of the second component of the eigenvalue 
Ev3 = zeros(4,11); % amplitude of the third component of the eigenvalue 
Ev4 = zeros(4,11); % amplitude of the fourth component of the eigenvalue 
eg_K_MJO = zeros(1,3); % the Kelvin wave component of the eigenvector of the MJO
eg_R_MJO = zeros(1,3); % the Rossby wave component of the eigenvector of the MJO
eg_Q_MJO = zeros(1,3); % the moisture component of the eigenvector of the MJO
eg_A_MJO = zeros(1,3); % the convective activity component of the eigenvector of the MJO
eg_K_Rossby = zeros(1,3); % the Kelvin wave component of the eigenvector of the moisture Rossby wave 
eg_R_Rossby = zeros(1,3); % the Rossby wave component of the eigenvector of the moisture Rossby wave 
eg_Q_Rossby = zeros(1,3); % the moisture component of the eigenvector of the moisture Rossby wave 
eg_A_Rossby = zeros(1,3); % the convective activity component of the eigenvector of the moisture Rossby wave 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Solving the eigenvalue problem %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for kk = 0:10
    if kk == 0
        kk = 0.01; % to prevent 0 appearing in the denominator
    end
    % matrix of the linear equation (the part generating oscillation
    % patterns)
    M_Omega = [30/17, 0, 0, 0;
        0, -10/17, 0, 0;        
        Q_tilde * 30/17, -Q_tilde*10/17, 0, 0;
        0, 0, 0, 0];
    % matrix of the linear equation (the part mostly for the damping)
    M_D = [-epsilon * d_bar * 30/17, 0, 0, -H_bar/2;
        0, -epsilon * d_bar * 30/17, 0, -H_bar/3;        
        0, 0, -epsilon * d_bar * 30/17, (-1+Q_tilde/6) * H_bar;
        0, 0, Gamma * A_bar, 0];
    % overall matrix
    M = M_Omega * 1i * 20 *kk/L*pi - M_D;
    % solving the eigenvalue problem
    [VecEig, Lambda] = eig(M);
    if kk == 0.01
        kk = 0;
    end
    % ranking the eigenvalues: from large positive to large negative values
    % MJO and Kelvin wave go eastward
    % Kelvin wave is dry and has a larger speed. So MJO ranks 2.
    % Moisture Rossby wave and (dry) Rossby wave go westward
    % Dry Rossby wave has a larger speed. So Moisture Rossby wave ranks 3.
    [sorting, ranking] = sort(diag(imag(Lambda)),'descend');
    Eig_Store(:,kk+1) = sorting;
    Ev1(:,kk+1) = sqrt((real(VecEig(:,ranking(1)))).^2 + (imag(VecEig(:,ranking(1)))).^2);
    Ev2(:,kk+1) = sqrt((real(VecEig(:,ranking(2)))).^2 + (imag(VecEig(:,ranking(2)))).^2);
    Ev3(:,kk+1) = sqrt((real(VecEig(:,ranking(3)))).^2 + (imag(VecEig(:,ranking(3)))).^2);
    Ev4(:,kk+1) = sqrt((real(VecEig(:,ranking(4)))).^2 + (imag(VecEig(:,ranking(4)))).^2);
    if kk >=1 && kk<=3
        eg_K_MJO(kk) = VecEig(1,ranking(2));
        eg_R_MJO(kk) = VecEig(2,ranking(2));
        eg_Q_MJO(kk) = VecEig(3,ranking(2));
        eg_A_MJO(kk) = VecEig(4,ranking(2));
        eg_K_Rossby(kk) = VecEig(1,ranking(3));
        eg_R_Rossby(kk) = VecEig(2,ranking(3));
        eg_Q_Rossby(kk) = VecEig(3,ranking(3));
        eg_A_Rossby(kk) = VecEig(4,ranking(3));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Display the eigen modes %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xx,yy] = meshgrid(0:10,1:4);
figure
KK = [0.01,1:10];
KK1 = [1:10];
for i = 1:4
    % frequency (CPD: Cycle per day)
    subplot(4,3,3*i-1)
    plot(KK, Eig_Store(i,:)/2/pi/dim_t ,'bo-','linewidth',2);    
    set(gca,'fontsize',12)
    if i == 1
        title('Frequency (cpd)','fontsize',12)
    end
    if i == 4
        xlabel('Wavenumber','fontsize',12);
    end
    % phase speed (m/s)
    subplot(4,3,3*i-2);
    plot(KK1, Eig_Store(i,2:end)/pi/dim_t*6250/27./KK1,'bo-','linewidth',2);
    set(gca,'fontsize',12)
    if i == 1
        ylabel('Dry Kelvin','fontsize',12);
    elseif i == 2
        ylabel('MJO','fontsize',12);
    elseif i == 3
        ylabel('Moist Rossby','fontsize',12);
    elseif i == 4
        ylabel('Dry Rossby','fontsize',12);
    end
    if i == 1
        title('Phase speed (m/s)','fontsize',12)        
    end
    if i == 4
        xlabel('Wavenumber','fontsize',12);
    end
end
% For each characteristic variable (e.g., eigenmode), the corresponding
% component of the original variables in the system. This is the eigenvector.
subplot(4,3,3);
amp = sqrt(real(Ev1).^2 + imag(Ev1).^2);
temp = amp(3,:);
amp(3,:) = amp(4,:);
amp(4,:) = temp;
contourf(xx,yy,amp,20,'linestyle','none');
set(gca,'fontsize',12)
set(gca,'YTickLabel',{'K','R','A','Q'})
title('Eigenvectors','fontsize',12)
subplot(4,3,6);
amp = sqrt(real(Ev2).^2 + imag(Ev2).^2);
temp = amp(3,:);
amp(3,:) = amp(4,:);
amp(4,:) = temp;
contourf(xx,yy,amp,20,'linestyle','none');
set(gca,'fontsize',12)
set(gca,'YTickLabel',{'K','R','A','Q'})
subplot(4,3,9);
amp = sqrt(real(Ev3).^2 + imag(Ev3).^2);
temp = amp(3,:);
amp(3,:) = amp(4,:);
amp(4,:) = temp;
contourf(xx,yy,amp,20,'linestyle','none');
set(gca,'fontsize',12)
set(gca,'YTickLabel',{'K','R','A','Q'})
subplot(4,3,12);
amp = sqrt(real(Ev4).^2 + imag(Ev4).^2);
temp = amp(3,:);
amp(3,:) = amp(4,:);
amp(4,:) = temp;
contourf(xx,yy,amp,20,'linestyle','none');
set(gca,'fontsize',12)
set(gca,'YTickLabel',{'K','R','A','Q'})
colormap jet
