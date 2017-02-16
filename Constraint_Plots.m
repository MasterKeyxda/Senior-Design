%% Constraint Plots

clc; 
clear all;
close all; 

% Outputs of this script are wing area (S) and engine thrust (T). 
% Five curves are generated based on aircraft stall speed, maximum speed,
% maximum rate of climb, take-off run, ceiling, and turn requirements.
% The intersection of the curves yields the design point and the
% corresponding wing loading (W/S) and thrust loading (T/W)

wingLoading = 50:0.5:250; % wing loading parameter (lb/ft^2)

%% Stall Speed Curve (4.3.2 Sadraey)

Vstall = 125; % stall speed (knots)
Vstall = Vstall * 1.688; % convert stall speed to ft/s
rhoSeaLevel = 0.002378; % density (slug/ft^3)
CLMax = 2.2; % Roskam Part I, p.91; CLmax range = 1.2 - 1.8

% Generate stall speed curve
stallSpeedCurve = 0.5*rhoSeaLevel*(Vstall^2)*CLMax;
stallSpeedCurve = stallSpeedCurve * ones(1, length(wingLoading)); % for plotting purposes
%% Maximum Speed Curve (4.3.3 Sadraey)

% Adjust density and speed of sound based on ceiling height

% Air Properties
altCruise = 41; % cruise ceiling altitude

% Call Function AltTable to retrieve density ratio and speed of sound ratio
[~,~,sigma,aSoundRatio] = AltTable(altCruise,'h'); % speed of sound ratio
rhoCruise = sigma * rhoSeaLevel; % density (slug/ft^3) at cruise ceiling
aSound = aSoundRatio * 1116.5; % speed of sound (ft/s) at cruise ceiling

% Cruise and Max Speeds
M = 1.6; % cruise Mach number
Vcruise = M*aSound; % cruise speed (ft/s)
Vmax = 1.1*Vcruise; % max speed (ft/s); assume max speed is 20% greater than
% cruise speed (Sadraey Eqn 4.49)

% Induced drag term, K
e = 0.9; % Oswald efficiency factor 
AR = 3; % aspect ratio 
K = 1 / (pi*e*AR); 

% Zero Lift Drag Coefficient
CD_0 = 0.020; % from CD_0 Estimate.xlsx or see Table 4.12 in Sadraey

% Sadraey Eqn 4.47 (T/W) 
% maxSpeedCurve = (rhoSeaLevel .* (Vmax.^2) .* CD_0 .* (1 ./ (2 .* wingLoading))) + (((2.*K) ./ (rhoSeaLevel .* sigma .* (Vmax.^2))) .* (wingLoading)); 
maxSpeedCurve = (rhoSeaLevel.*(Vmax.^2).*CD_0)./(2.*wingLoading) + (2.*K.*wingLoading)./(rhoSeaLevel.*sigma.*(Vmax.^2));
%% Take-off Run (S_TO) Curve (Section 4.3.4 Sadraey)

takeoffRun = 6900; % NASA takeoff field length of less than 7000 ft; based
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
CLR = CLMax / 1.21; 

% Eqn 4.67a
CDG = dragCoeff_TakeOff - (frictionCoeff * liftCoeff_TakeOff);

takeOffRunCurve = (frictionCoeff - (frictionCoeff + (CDG ./ CLR)) .*...
    (exp(0.6 .* rhoCruise .* g .* takeoffRun .* (1 ./ wingLoading)))) ...
    ./ (1 - (exp(0.6 .*rhoCruise .*g .* takeoffRun .* (1 ./ wingLoading)))); % Missing CDg term

%% Rate of Climb (ROC) Curve (4.3.5 Sadraey)

liftDragRatioMax = 7; % max lift to drag ratio
ROC = 10000;  % rate of climb (ft/min)
ROC = ROC / 60; % convert ROC to ft/s
% Rate of Climb Curve
 ROCCurve = (ROC ./ sqrt((2 ./ (rhoCruise .* sqrt(CD_0 ./ K))) .* wingLoading))...
     + (1 ./ liftDragRatioMax);
%% Ceiling Curve(4.3.6 Sadraey)
%Keyur Edit - None.
% Curve based on service ceiling
altService = 42; % kft
[~,~,sigmaService,aSoundRatio] = AltTable(altService, 'h'); % height input of kft
rhoService = sigmaService * rhoSeaLevel; 
ROC_Service = 100; % rate of climb at service ceiling (ft/min)
ROC_Service = ROC_Service / 60; % convert ROC to ft/s
ceilingCurve = (ROC_Service ./ (sigmaService .* sqrt((2 ./ (rhoService .* sqrt(CD_0 ./ K)))...
    .* wingLoading))) + (1 ./ (sigmaService .* liftDragRatioMax));
%% Plot Curves

figure(1)
hold on; 
title('Wing Area (S) and Engine Thrust (T) Constraint Plot')
xlabel('Wing Loading, W/S (lbf / ft^2)')
ylabel('Thrust-to-Weight Ratio, T/W (lbf/lbf)')
% Stall Speed Curve 
stallVerticalAxis = linspace(0,2, length(wingLoading));
plot(stallSpeedCurve, stallVerticalAxis, '--c')
% Max Speed Curve
plot(wingLoading, maxSpeedCurve, 'r')
% Take off Run Curve
plot(wingLoading, takeOffRunCurve, 'b')
% ROC Curve
plot(wingLoading, ROCCurve, 'm')
% Ceiling Curve
plot(wingLoading, ceilingCurve, 'g')
         
legend('V_{stall}','V_{max}', 'Take-off', 'Rate of Climb', 'Ceiling')

%% Wing Area and Engine Thrust

WTO = 102000; % gross takeoff weight lbf
[x,y] = ginput(1);
% Design point manually identified from contraint plot 
designPoint = [x, y]; % [Wing Loading (W/S), Thrust-to-Weight Ratio (T/W)]

fprintf('The wing area, S is: %0.2f ft^2 \n', WTO / designPoint(1))
fprintf('The thrust, T is: %0.2f lbf \n', designPoint(2) * WTO)