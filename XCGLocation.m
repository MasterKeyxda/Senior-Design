%% XCG Location

% All references and moment arms in ft
% Xarm refers to the x dist from the airplane nose to the start of component /
% component MAC. 

% Lifting Components - referenced to desired CG
% Wing
mXcgi.Wing = [Wt.Struc.Wing, (cglocAC+0.25*cRootSub)];%0.60 * L_F; % distance to LE MAC of wing
% Tail
mXcgi.Tail = [Wt.Struc.HT+Wt.Struc.VT, TAIL.Lopt]; %
XLE_w = 80;% length from nose to leading edge of wing
% Fuselag
mXcgi.Fuselage = [Wt.Struc.Fuselage, -XLE_w + cglocAC + 0.40*L_F]; % ranges from 0.40 - 0.50 length of fuselage; Roskam Pt.5, p.114; rear fuselage mounted engines
% Engine
length_engine = 142/12; % ft length of PW TF33-P7
x_engine = 100; % location of inlet of engine from nose
mXcgi.Engine2 = [Wt.Pwr.Engine*(2/3),-XLE_w + cglocAC + x_engine + length_engine*.3]; % engine_cg can vary from .3 -.4 of length of engine from inlet 
mXcgi.Engine1 = [Wt.Pwr.Engine*(1/3),-XLE_w + cglocAC + x_engine + 1.2*length_engine + length_engine*.3]; % single engine
% Nacelle
nacLength = (2.4077*(constraints.req_Thr^0.3876))/12; % ft
nacDMax = 1.0827*(constraints.req_Thr^0.4134)/12;
mXcgi.Nacelle2 = [Wt.Struc.Nacelle*(2/3),-XLE_w + cglocAC + x_engine + nacLength*.38];
mXcgi.Nacelle1 = [Wt.Struc.Nacelle*(1/3),-XLE_w + cglocAC + x_engine + 1.2*length_engine + nacLength*.38];
% Nose and Main Gear
x_NoseGear = 30; % ft nosegear distance from nose
x_MainGear = 85.92; % ft main gear distance from nose
mXcgi.NoseGear = [Wt.Struc.NoseGear,-XLE_w + cglocAC + x_NoseGear];
mXcgi.MainGear = [Wt.Struc.MainGear,-XLE_w + cglocAC + x_MainGear];
% Fuel System
% Fuel tank should have minimal effect on CG loaded and empty
LengthFuelTank = 12; % ft
% Put fuel tank CG on CG
mXcgi.FuelSystem = [Wt.Pwr.FuelSystem,0]; 
% Propulsion System
%mXcgi.=[Wt.Struc,-XLE_w + cglocAC + x_engine + length_engine*.4];

% Also Omit Starting Systems

% Fixed Equipment Weight
% Flight Control System
x_cockpit = 25 + 5.5/3; % center of cockpit from nose
mXcgi.FCsys = [Wt.Feq.FCsysGD,-XLE_w + cglocAC + x_cockpit]; % In Cockpit

% Hydraulic System for Landing Gear
Nose_rat = Wt.Struc.NoseGear/(Wt.Struc.NoseGear+Wt.Struc.MainGear); % weight ratio of nose landing gear
Main_rat = Wt.Struc.MainGear/(Wt.Struc.NoseGear+Wt.Struc.MainGear); % weight ratio of main landing gear
mXcgi.HydraulicNose = [Wt.Feq.Hydraulic*Nose_rat,-XLE_w + cglocAC + x_NoseGear];
mXcgi.HydraulicMain = [Wt.Feq.Hydraulic*Main_rat,-XLE_w + cglocAC + x_MainGear];

% Electrical System
mXcgi.ElectricalSystem = [Wt.Feq.Iae,-XLE_w + cglocAC + x_cockpit]; % Under Cockpit

% Air Conditioning, Pressurazation, and de-icing systems
mXcgi.Api = [Wt.Feq.Api,-XLE_w + cglocAC + x_engine + length_engine*.4]; % located at engine

% Oxygen System
Wt.Oxygen = [Wt.Feq.OxygenGD,-XLE_w + cglocAC + (53.33+25)/2]; % CG placed in center of habitable cabin 

% Auxiliary Power Unit
% Placed at Engines
% mXcgi.Apu2 = [Wt.Feq.Apu*(2/3),-XLE_w + cglocAC + x_engine + length_engine*.4];
mxcgi.Apu1 = [Wt.Feq.Apu,-XLE_w + cglocAC + x_engine + length_engine*.4];

%Aircraft Furnishings
mXcgi.Furn =[Wt.Feq.Furn,(-XLE_w + cglocAC +53.33+25)/2]; % CG placed in center of habitable cabin

% Operational Weight( Food, Drinks, ETC)
mXcgi.Oper = [Wt.Feq.Oper,-XLE_w + cglocAC + 26]; % CG placed in front of cabin

% Paint Weight
mXcgi.Paint = [Wt.Feq.Paint,-XLE_w + cglocAC + 0.40*L_F]; % Fuselage CG paint

% TOTAL AIRCRAFT CG
CG_types = fieldnames(mXcgi);
moment = 0;
wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mXcgi.(CG_types{i})(1) .* mXcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mXcgi.(CG_types{i})(1);
 end
 cg = moment/wt;
 fprintf('The Difference Between Required and Actual CG %0.2f ft\n', cg);

% Xcg.Wing = mXcg_i.Wing + 0.37*WING.geom.MAC; % 35-42% MAC (Sadraey Table 11.2)
% mXcg_i.HT = 138; 
% Xcg.HT = mXcg_i.HT + (0.35*TAIL.ch); % 30-40 % HT MAC (Sadraey Table 11.2)
% mXcg_i.VT = 140; 
% Xcg.VT = mXcg_i.VT + (0.35*TAIL.cv); % Is TAIL.cv the MAC? 30-40 % VT MAC (Sadraey Table 11.2)
%
% XLE_w = 72.5;
% mXcgi.Fuselage = [Wt.Struc.Fuselage, -XLE_w + cglocAC + 0.45*L_F]; % ranges from 0.45 - 0.50 length of fuselage; Roskam Pt.5, p.114; rear fuselage mounted engines
% mXcgi.Engines = [Wt.Pwr.Engine, 
% Nacelle Length
% http://adg.stanford.edu/aa241/AircraftDesign.html Sect 9.2.2)
% nacLength = (2.4077*(constraints.req_Thr^0.3876))/12; % ft
% nacDMax = 1.0827*(constraints.req_Thr^0.4134)/12;
