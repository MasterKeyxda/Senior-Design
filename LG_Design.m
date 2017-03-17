%% Landing Gear Preliminary Design 
% Parameters to Determine:
% wheel track, wheel Base, height, sizing, placement, load distribution,
% tires size, pressure, strut diameter, number of struts, strut length

% Date: 3/16/17
% Configuration: Tricycle
% Landing Gear Type: Retractable

% References:
% Roskam Airplane Design Part 2
% AIAA Journal of Aircraft: Undercarriage Lateral Tip-over Criteria
% Aircraft Design: A Conceptual Approach

%% Inputs
% Dimensions in ft
% CG Locations
% Xcg locations referenced from nose of landing gear
Xcg = 92; % Location of Xcg 
XcgAft = Xcg + 4; % Most aft Xcg location 
XcgFwd = Xcg - 4; % Most forward Xcg location

% Zcg locations 
ZcgFuse = 4.83; % zcg from bottom of fuselage
Zcg = 5 + ZcgFuse; % referenced from ground

% Landing Gear Geometry
LG.geom.xNose = 30; % distance from nose of airplane to nose LG
LG.geom.B = 98; % LG wheel base (dist between nose gear and main gear)
LG.geom.T = 9; % LG wheel track (dist between main gear wheels)
% XLE_w = 80; % length from nose to leading edge of wing

%% Lateral Tip-Over Criteria 

% phi --> Angle (deg) from nose wheel to main wheel
LG.angle.phi = atand(LG.geom.T / (2 * LG.geom.B));
LG.geom.OE = XcgFwd * sind(LG.angle.phi);
LG.angle.psi = atand(Zcg / LG.geom.OE); % lateral tipover angle (deg)

% Recommended that psi <= 55 deg for most fwd cg location
% However, it is acceptable for psi to range from 40 to 70 deg
fprintf('The lateral tip-over angle is %0.2f degrees. \n', LG.angle.psi)
if LG.angle.psi <= 55
    fprintf('Recommended lateral tip-over criterion met. \n')
elseif (LG.angle.psi >= 40 && LG.angle.psi <= 70)
    fprintf('Acceptable lateral tip-over criterion met. \n')
else
    fprintf('The lateral tip-over criterion is NOT satisfied. \n')
end
fprintf('\n')

%% Longitudinal Tip-Over Criteria

LG.angle.zeta = atand((LG.geom.B - XcgAft) / Zcg);

% Recommended that zeta is approx 15 deg (+-1 deg) for most aft cg location
% However, it is acceptable for zeta to range from 10 to 20 deg
fprintf('The longitudinal tip-over angle is %0.2f degrees. \n', LG.angle.zeta)
if (LG.angle.zeta >= 14) && (LG.angle.zeta <=16);
    fprintf('Recommended longitudinal tip-over criterion met. \n')
elseif (LG.angle.zeta >= 10) && (LG.angle.zeta <=20); 
    fprintf('Acceptable longitudinal tip-over criterion met. \n')
else
    fprintf('The longitudinal tip-over criterion is NOT satisfied. \n')
end
fprintf('\n')

%% Lateral Ground Clearance Criteria

% Front View: distance from main gear wheel assembly to wing tip
LG.geom.xLat = (WING.geom.span / 2) - (LG.geom.T / 2);
LG.geom.hClear = 6; % distance from ground to wing
LG.angle.beta = atand(LG.geom.hClear / LG.geom.xLat);

% Beta must be greater than 5 degrees
fprintf('The lateral ground clearance angle, beta is %0.2f degrees. \n', LG.angle.beta)
if (LG.angle.beta >= 5)
    fprintf('The lateral ground clearance is criterion met. \n')
else
    fprintf('The lateral ground clearance is NOT criterion met. \n')
end
fprintf('\n')

%% Longitudinal Ground Clearance Criteria

% Ensure aircraft tail does not strike ground during takeoff
LG.angle.takeOff = 10; % aircraft take-off angle 

% x Distance from main gear to end of tail
LG.geom.xMGAft = L_F - LG.geom.xNose - LG.geom.B; 
LG.geom.hTail = 8; % distance from ground to tail of fuselage

% Longitudinal Ground Clearance Angle
LG.angle.theta = atand(LG.geom.hTail / LG.geom.xMGAft);  
fprintf('The longitudinal ground clearance angle, beta is %0.2f degrees. \n', LG.angle.theta)
if (LG.angle.theta > LG.angle.takeOff + 1) && (LG.angle.theta <= 15)
    fprintf('The longitudinal ground clearance criterion is met. \n')
else 
    fprintf('The longitudinal ground clearance criterion is NOT met. \n')
end
fprintf('\n')

%% Landing Gear Loads

% See Figure 11.5 Aircraft Conceptual Design, Raymer (p. 234)
LG.geom.Na = XcgAft; % dist from nose gear to aft cg
LG.geom.Nf = XcgFwd; % dist from nose gear to fwd cg
LG.geom.Mf = LG.geom.B - LG.geom.Nf; % dist from fwd cg to main gear

% Number of tires
LG.qty.MainTire = 2; % number of tires on each main gear
LG.qty.NoseTire = 2; % number of tires on nose gear

% Number of Struts
LG.qty.Mainstruts = 2; % number of struts for entire main gear assembly
LG.qty.Nosestruts = 1; % number of struts for nose gear assembly

% Max Static Load (Main Gear) --> eqn 11.1
LG.load.MGMax = Wt.WTO * (LG.geom.Na / LG.geom.B); 
LG.load.MGMax  = LG.load.MGMax / LG.qty.Mainstruts; % Divide load by # of struts
% Load on each wheel; mult by 1.07 for FAR 25
LG.load.MGwheel = (LG.load.MGMax / LG.qty.MainTire)  * 1.07; % load on each wheel

% Max Static Load (Nose Gear) --> eqn 11.2
LG.load.NGMaxStatic = Wt.WTO * (LG.geom.Mf / LG.geom.B); 
LG.load.NGMaxStatic = LG.load.NGMaxStatic / LG.qty.Nosestruts; % Divide load by # of struts

% Nose Gear Dynamic Braking Load --> eqn 11.4
% Assumes braking coeff of 0.3 (hard runway); deceleartion of 10 ft/s^2
g = 32.17; % accel of gravity (ft/s^2)
LG.load.NGBrakeDyn = (10 * Zcg * Wt.WTO) / (g * LG.geom.B); 
LG.load.NGBrakeDyn = LG.load.NGBrakeDyn / LG.qty.Nosestruts; % Divide load by # of struts

% Determine whether to size nose gear for dynamic or static load
LG.load.NG = max(LG.load.NGBrakeDyn, LG.load.NGMaxStatic);

% Load on each nose wheel; mult by 1.07 for FAR 25 requirements
LG.load.NGwheel = (LG.load.NG / LG.qty.NoseTire) * 1.07;

%% Tire Sizing

%% Strut Sizing

% Strut Length
LG.strut.WL





