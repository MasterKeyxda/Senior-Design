%% XCG Location
clear CG_types
clear mXcgi

% All references and moment arms in ft
% Lifting Components - referenced to desired CG

% Wing
mXcgi.Wing = [Wt.Struc.Wing, cglocAC]; 

% Tail
mXcgi.Tail = [Wt.Struc.HT+Wt.Struc.VT, TAIL.Lopt]; 
XLE_w = 80; % length from nose to leading edge of wing

% Fuselage
% Ranges from 0.40 - 0.50 length of fuselage;
mXcgi.Fuselage = [Wt.Struc.Fuselage, -XLE_w - 0.25*cRootSub + cglocAC + 0.40*L_F];

% Engine
length_engine = 142/12; % ft length of PW TF33-P7
x_engine = 101.25; % location of inlet of engine from nose
% Engine_cg can vary from .3 -.4 of length of engine from inlet 
mXcgi.Engine2 = [Wt.Pwr.Engine*(2/3),-XLE_w - 0.25*cRootSub + cglocAC + x_engine + length_engine*.3]; 
% Single engine (most aft)
mXcgi.Engine1 = [Wt.Pwr.Engine*(1/3),-XLE_w - 0.25*cRootSub + cglocAC + x_engine + 1.3*length_engine + length_engine*.3]; 

% Nacelle
% http://adg.stanford.edu/aa241/AircraftDesign.html Sect 9.2.2)
nacLength = (2.4077*(constraints.req_Thr^0.3876))/12; % nacelle length ft
nacDMax = 1.0827*(constraints.req_Thr^0.4134)/12; % nacelle diameter ft
mXcgi.Nacelle2 = [Wt.Struc.Nacelle*(2/3),-XLE_w - 0.25*cRootSub + cglocAC + x_engine + nacLength*.38];
mXcgi.Nacelle1 = [Wt.Struc.Nacelle*(1/3),-XLE_w - 0.25*cRootSub + cglocAC + x_engine + 1.2*length_engine + nacLength*.38];

% Nose and Main Gear
x_NoseGear = 30; % ft nosegear distance from nose
x_MainGear = 85.92; % ft main gear distance from nose
mXcgi.NoseGear = [Wt.Struc.NoseGear,-XLE_w - 0.25*cRootSub + cglocAC + x_NoseGear];
mXcgi.MainGear = [Wt.Struc.MainGear,-XLE_w - 0.25*cRootSub + cglocAC + x_MainGear];

% Fuel System
% Fuel tank should have minimal effect on CG loaded and empty
LengthFuelTank = 12; % ft

% Put fuel tank CG on CG
mXcgi.FuelSystem = [Wt.Pwr.FuelSystem,0]; 

% Fixed Equipment Weight
% Flight Control System
x_cockpit = 25 + 5.5/3; % center of cockpit from nose
mXcgi.FCsys = [Wt.Feq.FCsysGD,-XLE_w - 0.25*cRootSub + cglocAC + x_cockpit]; % In Cockpit

% Hydraulic System for Landing Gear
Nose_rat = Wt.Struc.NoseGear/(Wt.Struc.NoseGear+Wt.Struc.MainGear); % weight ratio of nose landing gear
Main_rat = Wt.Struc.MainGear/(Wt.Struc.NoseGear+Wt.Struc.MainGear); % weight ratio of main landing gear
mXcgi.HydraulicNose = [Wt.Feq.Hydraulic*Nose_rat,-XLE_w  - 0.25*cRootSub + cglocAC + x_NoseGear];
mXcgi.HydraulicMain = [Wt.Feq.Hydraulic*Main_rat,-XLE_w  - 0.25*cRootSub + cglocAC + x_MainGear];

% Instrumentation, Avionics
mXcgi.Instrumentation = [Wt.Feq.Iae,-XLE_w - 0.25*cRootSub + cglocAC + x_cockpit]; % Under Cockpit

% Electrical System 
% mXcgi.ElecSys = [Wt. 

% Air Conditioning, Pressurazation, and de-icing systems
mXcgi.Api = [Wt.Feq.ApiGD,-XLE_w - 0.25*cRootSub + cglocAC + x_engine + length_engine*.4]; % located at engine

% Oxygen System
% CG placed in center of habitable cabin 
mXcgi.Oxygen = [Wt.Feq.OxygenGD,-XLE_w - 0.25*cRootSub + cglocAC + (53.33+25)/2]; 

% Auxiliary Power Unit
% Placed at Engines
mXcgi.Apu1 = [Wt.Feq.Apu,-XLE_w - 0.25*cRootSub + cglocAC + x_engine + length_engine*.4];

% Aircraft Furnishings (Seats, Lavatory, etc.)
% CG placed in center of habitable cabin
mXcgi.Furn =[Wt.Feq.Furn,(-XLE_w - 0.25*cRootSub + cglocAC +53.33+25)/2]; 

% Operational Weight( Food, Drinks, ETC)
mXcgi.Oper = [Wt.Feq.Oper,-XLE_w - 0.25*cRootSub + cglocAC + 26]; % CG placed in front of cabin

% Paint Weight
mXcgi.Paint = [Wt.Feq.Paint,-XLE_w - 0.25*cRootSub + cglocAC + 0.40*L_F]; % Fuselage CG paint

% Empty Weight CG
CG_types = fieldnames(mXcgi);
moment = 0;
wt = 0;

for i = 1:numel(CG_types)
   moment = moment - mXcgi.(CG_types{i})(1) .* mXcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mXcgi.(CG_types{i})(1);
end

wtEmpty = wt; 
cgEmpty = moment/wt;

% If positive - nose heavy
% If negative - tail heavy
fprintf('The empty weight cg is %0.4f \n', cgEmpty)
fprintf('The Difference Between Required and Actual CG %0.4f ft \n', cgEmpty);

%% Empty Weight + Crew CG (OEW CG)
mXcgi.Crew = [Wt.oew.crew,-XLE_w - 0.25*cRootSub + cglocAC + x_cockpit]; % crew in cockpit 
CG_types = fieldnames(mXcgi);
moment = 0;
wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mXcgi.(CG_types{i})(1) .* mXcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mXcgi.(CG_types{i})(1);
 end
wtOEW = wt; 
cgOEW = moment/wt;
fprintf('The operating empty weight cg is %0.4f \n', cgOEW)
%% Empty Weight + Crew + Fuel CG


% How are we distributing the fuel again?

%% Empty Weight + Crew + PAX (No baggage) CG

% Assume passenger distribution places cg at cabin center
mXcgi.pldNoBag = [Wt.pld.n_pass * Wt.pld.apw, -XLE_w - 0.25*cRootSub + cglocAC + (53.33+25)/2]; 
CG_types = fieldnames(mXcgi);
moment = 0;
wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mXcgi.(CG_types{i})(1) .* mXcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mXcgi.(CG_types{i})(1);
 end
wtPAX_NoBag = wt; 
cgPAX_NoBag = moment/wt;
fprintf('The cg for OEW + PAX (No baggage) is %0.4f \n', cgPAX_NoBag)
%% Empty Weight + Crew + PAX (w/ baggage) CG

% Assume passenger distribution places cg at cabin center
mXcgi.pldBag = [Wt.pld.w_tot, -XLE_w - 0.25*cRootSub + cglocAC + (53.33+25)/2]; 
CG_types = fieldnames(mXcgi);
moment = 0;
wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mXcgi.(CG_types{i})(1) .* mXcgi.(CG_types{i})(2); % Convert sign NEGATIVE = FLIP OVER
   wt = wt + mXcgi.(CG_types{i})(1);
 end
wtPAX_Bag = wt; 
cgPAX_Bag = moment/wt;
fprintf('The cg for OEW + PAX (with baggage) is %0.4f \n', cgPAX_Bag)

%% Empty Weight + Crew + PAX (No baggage) + Fuel  CG

%% Empty Weight + Crew + PAX (w/ baggage) + Fuel (MTOW) CG

%% CG Excursion Diagram

wtArray = [wtEmpty; wtOEW; wtPAX_NoBag; wtPAX_Bag];
cgArray = [cgEmpty; cgOEW; cgPAX_NoBag; cgPAX_Bag]; 

figure()
title('Aircraft CG Excursion')
xlabel('CG Location, ft')
ylabel('Weight (lbs)')
hold on;
% Plot CG Diagram
for iPlot = 1:length(wtArray)
    plot(cgArray(iPlot),wtArray(iPlot), 'o')
end
plot(cgArray,wtArray, 'b--') % plot line through points