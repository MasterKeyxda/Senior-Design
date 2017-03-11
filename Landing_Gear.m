%% Landing Gear Preliminary Design 
% Parameters to Determine:
% Wheel Track, Wheel Base, Height, Sizing, Location, Load Distribution,
% Tires size and pressure

%% META
% 3/5/17
% Determined Tricycle Landing Configuration
% Parameters Needed for Progression:
% 1. Acceleration/ Decelleration for TO/L
% 2. CG Locations
% 3. Fuel Location

%% CG Calculation
% See Sadraey, Figure 11.11 p.598 
x_LE = 72.5; % ft distance from nose to leading edge at root
x_MACpercent = .16; % Assume aircraft CG is located at 16% MAC
x_MAC_cg = x_MACpercent * WING.geom.MAC; % aircraft cg location in terms of MAC
C_r = 23.425; %  root chord (ft)
x_RootLE_MACLE = (C_r - WING.geom.MAC)/2; % Distance from LE Root to LE MAC
x_cg = x_MAC_cg + x_RootLE_MACLE + x_LE; % xcg relative to aircraft nose
x_cgShift = 8.33; % Based on extreme value Roskam Pt 2, p.243, Table 10.3
x_cgRange = [x_cg-x_cgShift,x_cg+x_cgShift]; % [for cg, aft cg]

%% LANDING GEAR LOAD DISTRIBUTION

% Main gear and nose gear placement
LG.geom.x_MainG = C_r + x_LE - 10; % distance from nose of airplane to center of main gear (ft)
LG.geom.Bm = LG.geom.x_MainG - x_cg; % distance between main gear and aircraft CG (ft)
LG.geom.x_Nose = 30; % distance from nose to nose gear (ft)
LG.geom.Bn = x_cg - LG.geom.x_Nose; % distance from cg to nose gear (ft)
LG.geom.B = LG.geom.Bm + LG.geom.Bn; % distance between nose gear and main gear, wheel base (ft)

% Percentage of Static Load
LG.load.Fn = LG.geom.Bm  ./ LG.geom.B * Wt.WTO; % Force of Static Load on Nose Gear (lbf)
LG.load.Fm = LG.geom.Bn ./ LG.geom.B * Wt.WTO; % Force of Static Load on Main Gear (lbf)
PercentFn = LG.load.Fn / Wt.WTO * 100; % percentage of loading on nose gear
PercentFm = LG.load.Fm / Wt.WTO * 100; % percentage of loading on main gear

fprintf('WHEEL BASE and LOADS \n')
fprintf('The wheel base is %0.2f ft \n',LG.geom.B)
fprintf('The load percentage on the nose gear is %0.2f percent \n',PercentFn)
fprintf('The load percentage on the main gear is %0.2f percent \n',PercentFm)
fprintf('\n')


% See Sadraey, Figure 9.18, p.506
% UDPATE ONCE Cg excursion known
LG.geom.height = 6.50; % height from bottom of tire to zcg location
% Nose Gear distances
LG.geom.BnMin = LG.geom.Bn - x_cgShift; % dist from nose gear to fwd cg (ft)
LG.geom.BnMax = LG.geom.Bn + x_cgShift; % dist from nose gear to aft cg (ft)
% Main Gear Distances
LG.geom.BmMin = LG.geom.B - LG.geom.BnMax; % dist from main gear to aft cg (ft)
LG.geom.BmMax = LG.geom.B - LG.geom.BnMin; % dist from main gear to fwd cg (ft)

% Gear Static and Dynamic Loading
g = 32.17; % acceleration of gravity (ft/s^2)
aLanding = 3; % assume deceleration for landing (Sadraey)
aTakeoff = 4; % assumed takeoff accel
% Gear Loads (lbf)
LG.load.FnMax = ((LG.geom.BmMax / LG.geom.B) * Wt.WTO) + ((aLanding * Wt.WTO * LG.geom.height) / (g * LG.geom.B)); 
LG.load.FmMax = ((LG.geom.BnMax / LG.geom.B) * Wt.WTO) + ((aTakeoff * Wt.WTO * LG.geom.height) / (g * LG.geom.B)); 

%% LANDING GEAR HEIGHT
% Sadraey Example 9.1 p. 502
% Requirements:
% Landing Gear Height Provides Clearance During Taxi
% Landing Gear Provides Clearance During Takeoff Rotation
% Landing Gear Reccommended Clearances for Rear Fuselage and Wing (AVG)
% Rear Fuselage: 1.14 ft, Wing: 2.46 ft
LG.geom.x_Upsweep = 120.83; % distance from nose to upsweep rear fuselage (ft)
LG.geom.alpha_TO = 10; % takeoff angle degrees
LG.geom.alpha_c = 10 +.5; % clearance angle
AB = LG.geom.x_Upsweep - LG.geom.x_MainG; % distance from main landing gear to upsweep (ft)
Hf = AB*tand(LG.geom.alpha_c); % distance between bottom tire of LG and fuselage

fprintf('LANDING GEAR HEIGHT \n')
fprintf('The landing gear height is Hf %0.2f ft \n',Hf) 
fprintf('\n')

%% Tip-Over Criteria

% Longitudinal Tip-Over Criterion
% Roskam Pt 2, p.219, Fig 9.1a
% Angle between most aft cg and vertical line extending from center of main
% gear wheel (degrees)
LG.geom.beta = atand(LG.geom.Bm / LG.geom.height); 
LG.crit.LT_Tip = 15; % longitudinal tip-over criteria, min angle (degrees) 
% Check if criterion is satisfied
fprintf('TIP-OVER CRITERIA \n')
fprintf('The  longitudinal tip-back angle is %0.2f deg \n', LG.geom.beta)
if LG.geom.beta > LG.crit.LT_Tip; 
    fprintf('The longitudinal tip-over criterion is satisfied. \n')
else
    fprintf('The longitudinal tip-over criterion is NOT satisfied. \n')
end
fprintf('\n')

% Lateral Tip-Over Criterion
LG.crit.Lat_Tip = 55; % lateral tip-over criteria, max angle (degrees)


%% Ground Clearance Criteria

% Note: Both Longitudinal and Lateral ground clearance requirements must be
% satisfied for tires and struts deflated. 

% Longitudinal Ground Clearance Criterion
LG.geom.tail = 15; % dist from ground to tail at end of fuselage (ft)
LG.geom.x_MGaft = L_F - LG.geom.B - LG.geom.x_Nose; % dist from main gear to end of fuselage (ft)
LG.crit.LT_Ground = LG.geom.alpha_TO; % takeoff angle 
LG.geom.theta = atand(LG.geom.tail / LG.geom.x_MGaft); % tail clearance angle (degrees)
% Check if criterion is satisfied
fprintf('GROUND CLEARANCE CRITERIA \n')
fprintf('The longitudinal ground clearance angle is %0.2f deg \n', LG.geom.theta)
if LG.geom.theta > LG.crit.LT_Ground;
    fprintf('The longitudinal ground clearance criterion is satisfied. \n')
else 
    fprintf('The longitudinal ground clearance criterion is NOT satisfied. \n')
end
fprintf('\n')
% Lateral Ground Clearance Criterion

%% Tire Selection

% Number of tires
LG.geom.noseTire = 2; % nose gear
LG.geom.mainTire = 4; % main gear (entire main gear assembly)

% Load on each tire / wheel
LG.load.noseT = LG.load.FnMax / LG.geom.noseTire;
LG.load.mainT = LG.load.FmMax / LG.geom.mainTire; 

% Tire Sizing
% Table 11.1 Raymer (dimensions in inches)
% Constants for Business Twin
LG.diam.A = 2.69;
LG.diam.B = 0.251; 

% Nose Gear Wheel Diameter
LG.geom.noseDiam = LG.diam.A * (LG.load.noseT^LG.diam.B);
LG.geom.noseDiam = LG.geom.noseDiam / 12; % convert to ft

% Nose Gear Wheel Width

% Main Gear Wheel Diameter

% Main Gear Wheel Width

fprintf('TIRE SELECTION \n')
fprintf('The number of tires on the nose gear assembly is %0.0f \n', LG.geom.noseTire)
fprintf('The number of tires on the main gear assembly is %0.0f \n', LG.geom.mainTire)



% See Table 11.2 in Raymer, Aircraft Design
