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

if ~exist('Wt', 'var')
   clc;
   clear;
   load('aircraft_vars.mat'); 
   close all;
end
clc;
%% Inputs
% Dimensions in ft
% CG Locations
% See figure 2 of https://arc.aiaa.org/doi/abs/10.2514/1.46564
% Xcg locations referenced from nose of landing gear
Xcg = cgExcel(1) - x_NoseGear;
XcgAft = Xcg; % Most aft Xcg location 
XcgFwd = Xcg - cgRange; % Most forward Xcg location

% Zcg locations 
%LG.geom.hClear = 0.7; % min clearance with lowest object
ZcgFuse = abs(ZcgEmpty); % Zcg from the bottom of the fuselage (ft)
LG.geom.hClear = 5.54; % clearance between fuselage and ground (same as lowest pt of wing)(ft)
Zcg = LG.geom.hClear + ZcgFuse; % referenced from ground

% Landing Gear Geometry
LG.geom.B = 64;
LG.geom.CO = LG.geom.B - XcgAft; % distance between MG wheels and CG;
%LG.geom.B = x_MainGear - x_NoseGear; % LG wheel base (dist between nose gear and main gear)
LG.geom.T = 9; % LG wheel track (dist between main gear wheels)

fprintf('The wheel base is %0.2f ft \n', LG.geom.B)
fprintf('The distance between MG wheels and the CG is %0.2f ft \n', LG.geom.CO)
fprintf('The wheel track is %0.2f ft \n', LG.geom.T)
fprintf('\n')
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
%LG.angle.zeta = atand((1.4) / Zcg);
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
LG.angle.takeOff = 8; % aircraft take-off angle 

% x Distance from main gear to end of tail
LG.geom.xMGAft = L_F - x_NoseGear - LG.geom.B; 
LG.geom.hTail = LG.geom.hClear + (Dmax/2); % distance from ground to tail of fuselage

% Longitudinal Ground Clearance Angle
LG.angle.theta = atand(LG.geom.hTail / LG.geom.xMGAft);  
fprintf('The takeoff angle is %0.2f degrees. \n', LG.angle.takeOff)
fprintf('The longitudinal ground clearance angle, theta is %0.2f degrees. \n', LG.angle.theta)
if (LG.angle.theta >= LG.angle.takeOff + 0.15) && (LG.angle.theta <= 15)
    fprintf('The longitudinal ground clearance criterion is met. \n')
else 
    fprintf('The longitudinal ground clearance criterion is NOT met. \n')
end
fprintf('\n')

%% Landing Gear Loads

% Based on Raymer Ch 11
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
% Assumes braking coeff of 0.3 (hard runway); size for deceleration of 10 ft/s^2
g = 32.17; % accel of gravity (ft/s^2)
LG.load.NGBrakeDyn = (10 * Zcg * Wt.WTO) / (g * LG.geom.B); 
LG.load.NGBrakeDyn = LG.load.NGBrakeDyn / LG.qty.Nosestruts; % Divide load by # of struts

% Nose Gear Load --> Static + Dynamic load (see p.31 Roskam Pt4)
LG.load.NG = LG.load.NGBrakeDyn + LG.load.NGMaxStatic;

% Load on each nose wheel; mult by 1.07 for FAR 25 requirements
LG.load.NGwheel = (LG.load.NG / LG.qty.NoseTire) * 1.07;

% Desirable load percentages: 5-20% Nose LG; 80-95% Main LG
LG.geom.Bm = LG.geom.B - Xcg; % distance between Xcg and MLG (ft)

% From statics 
% Load percentage for cg most aft
LG.load.NGpercAft = ((LG.geom.Bm / LG.geom.B) * 100); % NG load percentage
LG.load.MGpercAft = ((Xcg / LG.geom.B) * 100); % MLG load percentage
fprintf('Most AFT CG: The percentage of static load on the nose gear is %0.2f percent \n', LG.load.NGpercAft)
fprintf('Most AFT CG: The percentage of static load on the main gear is %0.2f percent \n', LG.load.MGpercAft)

% Load percentage for cg most forward
LG.geom.BmFwd = LG.geom.B - XcgFwd; 
LG.load.NGpercFwd = ((LG.geom.BmFwd / LG.geom.B) * 100); % NG load percentage
LG.load.MGpercFwd = ((XcgFwd / LG.geom.B) * 100); % MLG load percentage
fprintf('Most FWD CG: The percentage of static load on the nose gear is %0.2f percent \n', LG.load.NGpercFwd)
fprintf('Most FWD CG: The percentage of static load on the main gear is %0.2f percent \n', LG.load.MGpercFwd)
fprintf('\n')

%% Tire Sizing

% Plan to select from tire book once airplane design is further along
% For now, the tires will be sized using statistical methods
% Raymer, Table 11.1 Tire Constants for Business Twin
% Tire Diameter Constants
LG.const.ADiam = 2.69;
LG.const.BDiam = 0.251; 

% Tire Width Constants
LG.const.AWidth = 1.17; 
LG.const.BWidth = 0.216;

% Main Gear Assembly Tire Sizes (Dimensions in inches)
LG.tire.MainDiam = LG.const.ADiam *(LG.load.MGwheel^LG.const.BDiam);
LG.tire.MainWidth = LG.const.AWidth * (LG.load.MGwheel^LG.const.BWidth);

% Nose Gear Assembly Tire Sizes (Dimensions in inches)
% Raymer - "Nose tires can assumed to be 60-100% the size of main tires"
LG.tire.scale = 0.75; % scale nose tires dimensions by 70% of main tires
LG.tire.NoseDiam = LG.tire.scale * LG.tire.MainDiam;
LG.tire.NoseWidth = LG.tire.scale * LG.tire.MainWidth;

% Tire pressures

% Tire Parameters for later design 
% LG.tire.defl % max allowable tire deflection 

% Print sizing results
fprintf('TIRE SIZING \n')
fprintf('The tire diameter for the main gear assembly is %0.2f ft \n', LG.tire.MainDiam/12)
fprintf('The tire width for the main gear assembly is %0.2f ft \n', LG.tire.MainWidth/12)
fprintf('The tire diameter for the nose gear assembly is %0.2f ft \n', LG.tire.NoseDiam/12)
fprintf('The tire width for the nose gear assembly is %0.2f ft \n', LG.tire.NoseWidth/12)
fprintf('\n')
%% Strut Sizing

% Strut Parameters
% Based on method in Roskam Pt. 4, Section 2.5

% Strut Length --> eqn 2.11, 2.12 Roskam Pt.4
LG.tire.eff = 0.47; % tire energy absorption efficiency (Table 2.17 Roskam)
LG.strut.LdFactor = 1.5; % FAR 25 landing gear load factor (1.5 - 2.0)
LG.strut.Wtouch = 12; % design vertical touchdown rate (ft/s) --> FAR 25.723
LG.strut.eff = 0.75; % Raymer p.242 Table 11.4 Oleo-pneumatic)

% Main Landing Gear Strut
% MISSING TERMS 
LG.strut.WLMain =  LG.qty.Mainstruts * LG.load.MGMax; % MLG landing weight (lb) --> eqn 2.10
% Shock absorber efficiency 
% Shock absorber length (ft) --> eqn 2.12 Roskam Pt.4
%LG.strut.lengthMain = (0.5*(LG.strut.WLMain / g)*(LG.strut.Wtouch^2) / (LG.qty.Mainstruts * LG.load.MGMax * LG.strut.LdFactor)); 
LG.strut.lengthMain = LG.geom.hClear - LG.tire.MainDiam/12; % required fuselage clearance - tire diameter; 
% Strut Diameter --> eqn 2.13 Roskam Pt.4
LG.strut.diamMain = 0.041 + (0.0025 * (LG.load.MGMax)^0.5); % diameter of strut

% Nose Landing Gear Strut
LG.strut.WLNose = LG.qty.Nosestruts * LG.load.NG; % NLG landing weight (lb) --> eqn 2.10
% Nose Strut Diameter
LG.strut.diamNose = 0.041 + (0.0025 * (LG.load.NG)^0.5); % diameter of strut
fprintf('Strut Sizing \n')
fprintf('The main gear strut diameter is %0.2f ft \n',LG.strut.diamMain)
fprintf('The main gear strut length is %0.2f ft \n', LG.strut.lengthMain)
% fprintf('The nose gear strut length is %0.2f ft \n')
fprintf('The nose gear strut diameter is %0.2f ft \n', LG.strut.diamNose)



