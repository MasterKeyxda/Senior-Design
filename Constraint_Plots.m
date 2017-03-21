%% Constraint Plots
% 
% clc; 
% clear all;
% close all; 

% Outputs of this script are wing area (S) and engine thrust (T). 
% Five curves are generated based on aircraft stall speed, maximum speed,
% maximum rate of climb, take-off run, ceiling, and turn requirements.
% The intersection of the curves yields the design point and the
% corresponding wing loading (W/S) and thrust loading (T/W)

% allows running of script independently of main script
if ~exist('Wt', 'var')
   clc;
   clear;
   load('aircraft_vars.mat'); 
   close all;
end

constraints.wingLoading = 50:0.5:250; % wing loading parameter (lb/ft^2)

%% Stall Speed Curve (4.3.2 Sadraey)

constraints.Vstall = 125 * 1.688; % stall speed, knots -> ft/s
% constraints.rho_sl = 0.002378; % density (slug/ft^3)
constraints.CLMax = 2.0; % Roskam Part I, p.91; CLmax range = 1.2 - 1.8

% Generate stall speed curve
constraints.Vstall_Curve = 0.5*atm.rho_sl*(constraints.Vstall^2)*constraints.CLMax * ones(1, length(constraints.wingLoading)); % for plotting purposes

%% Maximum Speed Curve (4.3.3 Sadraey)

% Cruise and Max Speeds
constraints.Vmax = 1.1*Wt.fuel.V_max_cr; % max speed (ft/s); assume max speed is 10% greater than
% cruise speed (Sadraey Eqn 4.49)

% Induced drag term, K
e = 0.9; % Oswald efficiency factor 
AR = 3; % aspect ratio 
K = 1 / (pi*e*AR); 

% Zero Lift Drag Coefficient
CD_0 = 0.0200; % from CD_0 Estimate.xlsx or see Table 4.12 in Sadraey

% Sadraey Eqn 4.47 (T/W) 
% maxSpeedCurve = (rhoSeaLevel .* (Vmax.^2) .* CD_0 .* (1 ./ (2 .* wingLoading))) + (((2.*K) ./ (rhoSeaLevel .* sigma .* (Vmax.^2))) .* (wingLoading)); 
constraints.maxSpeedCurve = (atm.rho_sl.*(constraints.Vmax.^2).*CD_0)./(2.*constraints.wingLoading) + (2.*K.*constraints.wingLoading)./(atm.rho_sl.*atm.sig_rho.*(constraints.Vmax.^2));

%% Take-off Run (S_TO) Curve (Section 4.3.4 Sadraey)

% takeoffRun = 6900; % NASA takeoff field length of less than 7000 ft; based
% on FAR 25, aircraft must clear imaginary 35 ft obstacle

frictionCoeff = 0.04; % friction coefficient dry concrete / asphault (Sadraey Table 4.15)
g = 32.17; % gravity (ft/s^2)
% LG gear and Flaps drag coefficients (Equation 4.69b Sadraey)
Cd_landingGear = 0.08; % LG drag coefficient varies 0.06 - 0.012
Cd_flap = 0.005; % flap drag coefficient varies from 0.003 to 0.008

% Zero lift drag coefficient at takeoff (Eqn 4.69a)
Cd_0L_TO = CD_0 + Cd_landingGear + Cd_flap;

% Cruise lift coeff and flap coeff at take-off
CL_Cruise = 0.05; % typical value according to Sadraey
CL_flap_TO = 0.6; % varies from 0.3 to 0.8

% Take-off Lift Coefficient (Eqn 4.69c)
CL_TO = CL_Cruise + CL_flap_TO; 

% Drag coefficient at take-off (Eqn 4.68)
Cd_TO = Cd_0L_TO + (K * (CL_TO^2));

% Aircraft Lift at take-off (Eqn 4.67b)
CLR = constraints.CLMax / 1.21; 

% Eqn 4.67a
CDG = Cd_TO - (frictionCoeff * CL_TO);
rho_cr = atm.sig_rho * atm.rho_sl;
takeOffRunCurve = (frictionCoeff - (frictionCoeff + (CDG ./ CLR)) .*...
    (exp(0.6 .* rho_cr .* g .* req.takeoffRun .* (1 ./ constraints.wingLoading)))) ...
    ./ (1 - (exp(0.6 .*rho_cr .*g .* req.takeoffRun .* (1 ./ constraints.wingLoading)))); % Missing CDg term

%% Rate of Climb (ROC) Curve (4.3.5 Sadraey)

LD_Max = 7; % max lift to drag ratio
ROC = 10000 / 60;  % rate of climb (ft/min) -> ft/s
% Rate of Climb Curve
 ROCCurve = (ROC ./ sqrt((2 ./ (rho_cr .* sqrt(CD_0 ./ K))) .* constraints.wingLoading))...
     + (1 ./ LD_Max);
 
%% Ceiling Curve(4.3.6 Sadraey)
%Keyur Edit - None.
% Curve based on cruise ceiling (Eq. 4.95, Sadraey)
% altService = 42; % kft
% [~,~,sigmaService,aSoundRatio] = AltTable(altService, 'h'); % height input of kft
rhoService = atm.sig_rho * atm.rho_sl; 
ROC_cr = 300; % rate of climb at service ceiling (ft/min) (FAR numbers)
ROC_cr = ROC_cr / 60; % convert ROC to ft/s
ceilingCurve = (ROC_cr ./ (atm.sig_rho .* sqrt((2 ./ (rhoService .* sqrt(CD_0 ./ K)))...
    .* constraints.wingLoading))) + (1 ./ (atm.sig_rho .* LD_Max));

%% Plot Curves

figure();
hold on; 
title('Wing Area (S) and Engine Thrust (T) Constraint Plot')
xlabel('Wing Loading, W/S (lbf / ft^2)')
ylabel('Thrust-to-Weight Ratio, T/W (lbf/lbf)')
% Stall Speed Curve 
stallVerticalAxis = linspace(0,2, length(constraints.wingLoading));
plot(constraints.Vstall_Curve, stallVerticalAxis, '--c')
% Max Speed Curve
plot(constraints.wingLoading, constraints.maxSpeedCurve, 'r')
% Take off Run Curve
plot(constraints.wingLoading, takeOffRunCurve, 'b')
% ROC Curve
plot(constraints.wingLoading, ROCCurve, 'm')
% Ceiling Curve
plot(constraints.wingLoading, ceilingCurve, 'g')
         
legend('V_{stall}','V_{max}', 'Take-off', 'Rate of Climb', 'Ceiling')

%% Wing Area and Engine Thrust

% WTO = 102000; % gross takeoff weight lbf
% [x,y] = ginput(1);
% Design point manually identified from contraint plot 
% designPoint = [x, y]; % [Wing Loading (W/S), Thrust-to-Weight Ratio (T/W)]
test1 = [constraints.Vstall_Curve(1), spline(constraints.wingLoading, ceilingCurve, constraints.Vstall_Curve(1))];
% test2 -> where vstall and vmax curve intersect
test2 = [constraints.Vstall_Curve(1), spline(constraints.wingLoading, constraints.maxSpeedCurve, constraints.Vstall_Curve(1))];

% test3 -> intersection between Ceiling & Vmax curve
ceilCurve_avg = mean(ceilingCurve);
test3 = [spline(constraints.maxSpeedCurve, constraints.wingLoading, ceilCurve_avg), ceilCurve_avg];

if test3(1) <= test2(1)
    designPoint = test3;
else
    designPoint = test2;
end
fprintf('TW Ratio: %0.5f\n', designPoint(2));
fprintf('Wing Loading: %0.5f\n', designPoint(1));
fprintf('The wing area, S is: %0.2f ft^2 \n', Wt.WTO / designPoint(1))
fprintf('The thrust, T is: %0.2f lbf \n', ceilCurve_avg * Wt.WTO)

S_w = Wt.WTO/designPoint(1);
constraints.req_Thr = designPoint(2) * Wt.WTO;