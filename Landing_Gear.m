%% Landing Gear Preliminary Design 
% All references correspond to Sadraey
% Parameters to Determine:
% Wheel Track, Wheel Base, Height, Sizing, Location, Load Distribution,
% Tires

%% META
% 3/5/17
% Determined Tricycle Landing Configuration
% Parameters Needed for Progression:
% 1. Acceleration/ Decelleration for TO/L
% 2. CG Locations
% 3. Fuel Location

%% CG Calculation
% See Figure 11.11 pg 598 
alpha_TO = 12; % Takeoff angle degrees
alpha_c = 12.1; % Clearance angle degrees
x_LE = 72.5; % ft distance from nose to leading edge at root
x_MACpercent = .16; %  Aircraft CG is located at 16% MAC
x_MAC_cg = x_MACpercent * WING.geom.MAC; % aircraft cg location in terms of MAC
C_r = 23.425; %  root chord (ft)
x_RootLE_MACLE = (C_r - WING.geom.MAC)/2; % Distance from LE Root to LE MAC
x_cg = x_MAC_cg + x_RootLE_MACLE + x_LE; % CG of entire aircraft relative to nose
%x_cg = [0.8 1 1.2] * x_cg; % [forward takeoff aft] shift in CG during flight

%% Statics for Wheel Base
% First Starting With Placement of main landing gear
x_MainGear = C_r + x_LE - 10; % distance from nose of airplane to center of main gear (ft)
Bm = x_MainGear - x_cg; % distance between main gear and aircraft CG
x_NoseGear = 37.33; %(ft) distance from nose to nose gear
Bn = x_cg - x_NoseGear; % distance from cg to nose gear
B = Bm + Bn; % distance between nose gear and main gear

% Percentage of Static Load
Fn = Bm./B * Wt.WTO; % Force of Static Load on Nose Gear (lbf)
Fm = Bn./B * Wt.WTO; % Force of Static Load on Main Gear (lbf)
PercentFn = Fn/Wt.WTO * 100; % percentage of loading on nose gear
PercentFm = Fm/Wt.WTO * 100; % percentage of loading on main gear

fprintf('The Wheel Base is %0.2f ft \n',B)
fprintf('The load percentage on the nose gear is %0.2f percent \n',PercentFn)
fprintf('The load percentage on the main gear is %0.2f percent \n',PercentFm)

%% Determine Landing Gear Height
% EX 9.1 pg. 502
% Requirements:
% Landing Gear Height Provides Clearance During Taxi
% Landing Gear Provides Clearance During Takeoff Rotation
% Landing Gear Reccommended Clearances for Rear Fuselage and Wing (AVG)
% Rear Fuselage: 1.14 ft, Wing: 2.46 ft
x_upsweep = 120.83; % distance from nose to upsweep rear fuselage
alpha_TO = 10; % takeoff angle degrees
alpha_c = 10 +.3; % clearance angle
AB = x_upsweep - x_MainGear; % distance from aft landing gear to upsweep (ft)
Hf = AB*tand(alpha_c); % distance between landing gear and fuselage
fprintf('The landing gear height is Hf %0.2f ft \n',Hf) 
