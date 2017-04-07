clc;
clear;
close all;

% load('aircraft_vars.mat');

%% Calculate 2-D Oblique Shocks

beta = linspace(0, 90, 1001);
% theta = linspace(0, 20, 1001);
gam = 1.4;
M0 = 1.6; % req.cr_M0(1);
M2_des = 0.999;

%% First Shock

theta_1 = atand(2.*cotd(beta).*((M0.*sind(beta)).^2 - 1)./((M0^2.*(gam + cosd(2.*beta))) + 2));
Pt1 = (0.5*(gam + 1).*(M0.*sind(beta)).^2./(1 + 0.5*(gam - 1).*(M0 * sind(beta)).^2)).^(gam./(gam - 1)) .* ((2*gam/(gam+1)).*(M0.*sind(beta)).^2 - (gam - 1)./(gam + 1)).^(-1./(gam - 1));

beta1 = beta(theta_1 >= 0);
Pt1 = Pt1(theta_1 >= 0);
theta1 = theta_1(theta_1 >=0);
M1 = 1./(sind(beta1 - theta1)) .* sqrt((1 + 0.5*(gam - 1).*(M0.*sind(beta1)).^2)./(gam.*(M0.*sind(beta1)).^2 - 0.5*(gam - 1)));

%% Second shock

Pt_tot = [];
Pt2 = [];
M2 = [];
theta2 = [];
beta2 = [];
ind = [];

for i = 1:length(M1)
    % Calculate all possible ramp angles
    theta_i(i,:) =  atand(2.*cotd(beta).*((M1(i).*sind(beta)).^2 - 1)./((M1(i)^2.*(gam + cosd(2.*beta))) + 2));
    Pti(i,:) = (0.5*(gam + 1).*(M1(i).*sind(beta)).^2./(1 + 0.5*(gam - 1).*(M1(i) * sind(beta)).^2)).^(gam./(gam - 1)) .* ((2*gam/(gam+1)).*(M1(i).*sind(beta)).^2 - (gam - 1)./(gam + 1)).^(-1./(gam - 1));
   
    % filter out values that didn't work
    beta_i{i} = beta(theta_i(i,:) >= 0);
    Pt_i{i} = Pti(i, theta_i(i,:) >= 0);
    theta_n{i} = theta_i(i,theta_i(i,:) >=0);
    M_i{i} = 1./(sind(beta_i{i} - theta_n{i})) .* sqrt((1 + 0.5*(gam - 1).*(M1(i).*sind(beta_i{i})).^2)./(gam.*(M1(i).*sind(beta_i{i})).^2 - 0.5*(gam - 1)));

    % Find subsonic Mach Number
    if any(M_i{i} < M2_des)
        M2(end+1) = max(M_i{i}(M_i{i} < M2_des));
        Pt2(end+1) = (Pt_i{i}(M_i{i} == M2(end)));
        Pt_tot(end+1) = Pt1(i) .* Pt2(end);
        theta2(end+1) = theta_n{i}((M_i{i} == M2(end)));
        beta2(end+1) = beta_i{i}(M_i{i} == M2(end));
        ind(end+1) = i;
%     else
%         M2(i) = min(M_i{i});
%         Pt2(i) = Pt_i{i}(M_i{i} == M2(i));
%         Pt_tot(i) = Pt1(i) .* Pt2(i);
    end
    
end

%% Post-Process Results

figure();
plot(theta1(ind), Pt_tot);
title('Total Pressure Ratio Across Both Shocks');
xlabel('\theta_1');
ylabel('Pt_2/Pt_0');

figure();
plot(theta1(ind), theta2);
title('Ramp 1 vs Ramp 2 Angles');
xlabel('\theta_1');
ylabel('\theta_2');

theta1_temp = theta1(ind);
beta1_temp = beta1(ind);

theta1_des = theta1_temp(Pt_tot == max(Pt_tot));
theta2_des = theta2(Pt_tot == max(Pt_tot));
beta1_des = beta1_temp(Pt_tot == max(Pt_tot));
beta2_des = beta2(Pt_tot == max(Pt_tot));

fprintf('\tFirst Shock\n\n');
fprintf('Ramp Angle 1: %0.5f\n', theta1_des);
fprintf('Shock Angle 1: %0.5f\n', beta1_des);

fprintf('\n\tSecond Shock\n\n');
fprintf('Ramp Angle 2: %0.5f\n', theta2_des);
fprintf('Shock Angle 2: %0.5f\n', beta2_des);

fprintf('Max Pressure Recovery Ratio: %0.5f\n', max(Pt_tot));


%% Area
mdotreq = 266.56; %Taken from PERF
P0 =2.480729*144; %lbf/ft^2
R = 53.34; %lbf/lbmR
T = 389.97; % Rankine
gamma = 1.4; % Cp/Cv
gc = 32.174; % Fucking US Constants
M0 =1.6; %Flight Mach Number
area = mdotreq/(P0/(R*T)*sqrt(gamma*gc*R*T)*M0); %Inlet Area

fprintf('Inlet Area Required = %0.2f ft^2 \n', area);
