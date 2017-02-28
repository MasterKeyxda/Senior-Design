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
   close all;
   load('aircraft_vars.mat'); 
end

constraints.wingLoading = 50:0.5:250; % wing loading parameter (lb/ft^2)

%% Stall Speed Curve (4.3.2 Sadraey)

constraints.Vstall = 125 * 1.688; % stall speed, knots -> ft/s
% constraints.rho_sl = 0.002378; % density (slug/ft^3)
constraints.CLMax = 2.2; % Roskam Part I, p.91; CLmax range = 1.2 - 1.8

% Generate stall speed curve
constraints.Vstall_Curve = 0.5*atm.rho_sl*(constraints.Vstall^2)*constraints.CLMax * ones(1, length(constraints.wingLoading)); % for plotting purposes
%% Maximum Speed Curve (4.3.3 Sadraey)

% Cruise and Max Speeds
% Vcruise = M*a; % cruise speed (ft/s)
constraints.Vmax = 1.1*Wt.fuel.V_max_cr; % max speed (ft/s); assume max speed is 10% greater than
% cruise speed (Sadraey Eqn 4.49)

% Induced drag term, K
e = 0.9; % Oswald efficiency factor 
AR = 3; % aspect ratio 
K = 1 / (pi*e*AR); 

% Zero Lift Drag Coefficient
CD_0 = 0.020; % from CD_0 Estimate.xlsx or see Table 4.12 in Sadraey

% Sadraey Eqn 4.47 (T/W) 
% maxSpeedCurve = (rhoSeaLevel .* (Vmax.^2) .* CD_0 .* (1 ./ (2 .* wingLoading))) + (((2.*K) ./ (rhoSeaLevel .* sigma .* (Vmax.^2))) .* (wingLoading)); 
constraints.maxSpeedCurve = (atm.rho_sl.*(constraints.Vmax.^2).*CD_0)./(2.*constraints.wingLoading) + (2.*K.*constraints.wingLoading)./(atm.rho_sl.*atm.sig_rho.*(constraints.Vmax.^2));
%% Take-off Run (S_TO) Curve (Section 4.3.4 Sadraey)

% takeoffRun = 6900; % NASA takeoff field length of less than 7000 ft; based
% on FAR 25, aircraft must clear imaginary 35 ft obstacle

frictionCoeff = 0.04; % friction coefficient dry concrete / asphault (Sadraey Table 4.15)
g = 32.17; % gravity (ft/s^2)
% LG gear and Flaps drag coefficients (Equation 4.69b Sadraey)
dragCoeff_landingGear = 0.08; % LG drag coefficient varies 0.06 - 0.012
dragCoeff_flap = 0.005; % flap drag coefficient varies from 0.003 to 0.008

% Zero lift drag coefficient at takeoff (Eqn 4.69a)
dragCoeff_zeroLiftTakeOff = CD_0 + dragCoeff_landingGear + dragCoeff_flap;

% Cruise lift coeff and flap coeff at take-off
liftCoeff_Cruise = 0.05; % typical value according to Sadraey
liftCoeff_flapTakeOff = 0.6; % varies from 0.3 to 0.8

% Take-off Lift Coefficient (Eqn 4.69c)
liftCoeff_TakeOff = liftCoeff_Cruise + liftCoeff_flapTakeOff; 

% Drag coefficient at take-off (Eqn 4.68)
dragCoeff_TakeOff = dragCoeff_zeroLiftTakeOff + (K * (liftCoeff_TakeOff^2));

% Aircraft Lift at take-off (Eqn 4.67b)
CLR = constraints.CLMax / 1.21; 

% Eqn 4.67a
CDG = dragCoeff_TakeOff - (frictionCoeff * liftCoeff_TakeOff);
rho_cr = atm.sig_rho * atm.rho_sl;
takeOffRunCurve = (frictionCoeff - (frictionCoeff + (CDG ./ CLR)) .*...
    (exp(0.6 .* rho_cr .* g .* req.takeoffRun .* (1 ./ constraints.wingLoading)))) ...
    ./ (1 - (exp(0.6 .*rho_cr .*g .* req.takeoffRun .* (1 ./ constraints.wingLoading)))); % Missing CDg term

%% Rate of Climb (ROC) Curve (4.3.5 Sadraey)

liftDragRatioMax = 7; % max lift to drag ratio
ROC = 10000;  % rate of climb (ft/min)
ROC = ROC / 60; % convert ROC to ft/s
% Rate of Climb Curve
 ROCCurve = (ROC ./ sqrt((2 ./ (rho_cr .* sqrt(CD_0 ./ K))) .* constraints.wingLoading))...
     + (1 ./ liftDragRatioMax);
%% Ceiling Curve(4.3.6 Sadraey)
%Keyur Edit - None.
% Curve based on cruise ceiling
altService = 42; % kft
[~,~,sigmaService,aSoundRatio] = AltTable(altService, 'h'); % height input of kft
rhoService = sigmaService * atm.rho_sl; 
ROC_Service = 300; % rate of climb at service ceiling (ft/min) % FAR numbers
ROC_Service = ROC_Service / 60; % convert ROC to ft/s
ceilingCurve = (ROC_Service ./ (sigmaService .* sqrt((2 ./ (rhoService .* sqrt(CD_0 ./ K)))...
    .* constraints.wingLoading))) + (1 ./ (sigmaService .* liftDragRatioMax));
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
ceilCurve_avg = mean(ceilingCurve);
designPoint = spline(constraints.maxSpeedCurve, constraints.wingLoading, ceilCurve_avg);

fprintf('The wing area, S is: %0.2f ft^2 \n', WTO / designPoint)
fprintf('The thrust, T is: %0.2f lbf \n', ceilCurve_avg * WTO)

S_w = WTO/designPoint;
contstraints.req_Thr = ceilCurve_avg * WTO;