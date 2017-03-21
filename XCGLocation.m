%% XCG Location
clear CG_types
clear mXcgi

% All references and moment arms in ft
<<<<<<< Updated upstream
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
x_engine2 = 116; % location of inlet double engines from nose
x_engine = 1.5*length_engine+x_engine2; % location of inlet of single engine from nose

% Engine_cg can vary from .3 -.4 of length of engine from inlet 
mXcgi.Engine2 = [Wt.Pwr.Engine*(2/3),-XLE_w - 0.25*cRootSub + cglocAC + x_engine2 + length_engine*.3]; 
% Single engine (most aft)
mXcgi.Engine1 = [Wt.Pwr.Engine*(1/3),-XLE_w - 0.25*cRootSub + cglocAC + x_engine + 0.3*length_engine + length_engine*.3]; 

% Nacelle
% http://adg.stanford.edu/aa241/AircraftDesign.html Sect 9.2.2)
nacLength = (2.4077*(constraints.req_Thr^0.3876))/12; % nacelle length ft
nacDMax = 1.0827*(constraints.req_Thr^0.4134)/12; % nacelle diameter ft
mXcgi.Nacelle2 = [Wt.Struc.Nacelle*(2/3),-XLE_w - 0.25*cRootSub + cglocAC + x_engine2 + nacLength*.38];
mXcgi.Nacelle1 = [Wt.Struc.Nacelle*(1/3),-XLE_w - 0.25*cRootSub + cglocAC + x_engine + 0.2*length_engine + nacLength*.38];

% Nose and Main Gear
x_NoseGear = 32; % ft nosegear distance from nose
x_MainGear = 95.92; % ft main gear distance from nose
mXcgi.NoseGear = [Wt.Struc.NoseGear,-XLE_w - 0.25*cRootSub + cglocAC + x_NoseGear];
mXcgi.MainGear = [Wt.Struc.MainGear,-XLE_w - 0.25*cRootSub + cglocAC + x_MainGear];

% Fuel System
% Fuel tank should have minimal effect on CG loaded and empty
LengthFuelTank = 12; % ft

% Put fuel tank CG on CG
mXcgi.FuelSystem = [Wt.Pwr.FuelSystem,0]; 

% Fixed Equipment Weight
% Flight Control System

x_cockpit = L_N + L_CP/3; % length of quiet spike needle +  center of cockpit from nose
mXcgi.FCsys = [Wt.Feq.FCsysGD,-XLE_w - 0.25*cRootSub + cglocAC + x_cockpit]; % In Cockpit
=======
% Xarm refers to the x dist from the airplane nose to the start of component /
% component MAC. 

% Lifting Components - referenced to desired CG
% Wing
mXcgi.Wing = [Wt.Struc.Wing, (cglocAC+0.25*cRootSub)]; 
% Tail
mXcgi.Tail = [Wt.Struc.HT+Wt.Struc.VT, TAIL.Lopt]; %
XLE_w = 80;% length from nose to leading edge of wing
% Fuselage
mXcgi.Fuselage = [Wt.Struc.Fuselage, -XLE_w + cglocAC + 0.40*L_F]; % ranges from 0.40 - 0.50 length of fuselage;
% Engine
length_engine = 142/12; % ft length of PW TF33-P7
x_engine = 101.25; % location of inlet of engine from nose
mXcgi.Engine2 = [Wt.Pwr.Engine*(2/3),-XLE_w + cglocAC + x_engine + length_engine*.3]; % engine_cg can vary from .3 -.4 of length of engine from inlet 
mXcgi.Engine1 = [Wt.Pwr.Engine*(1/3),-XLE_w + cglocAC + x_engine + 1.3*length_engine + length_engine*.3]; % single engine
% Nacelle
% http://adg.stanford.edu/aa241/AircraftDesign.html Sect 9.2.2)
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
>>>>>>> Stashed changes

% Hydraulic System for Landing Gear
Nose_rat = Wt.Struc.NoseGear/(Wt.Struc.NoseGear+Wt.Struc.MainGear); % weight ratio of nose landing gear
Main_rat = Wt.Struc.MainGear/(Wt.Struc.NoseGear+Wt.Struc.MainGear); % weight ratio of main landing gear
<<<<<<< Updated upstream
mXcgi.HydraulicNose = [Wt.Feq.Hydraulic*Nose_rat,-XLE_w  - 0.25*cRootSub + cglocAC + x_NoseGear];
mXcgi.HydraulicMain = [Wt.Feq.Hydraulic*Main_rat,-XLE_w  - 0.25*cRootSub + cglocAC + x_MainGear];

% Instrumentation, Avionics
mXcgi.Instrumentation = [Wt.Feq.Iae,-XLE_w - 0.25*cRootSub + cglocAC + x_cockpit]; % Under Cockpit

% Electrical System 
mXcgi.ElecSys = [Wt.Feq.ElecSys,-XLE_w - 0.25*cRootSub + cglocAC + x_cockpit]; % Also under cockpit


% Air Conditioning, Pressurazation, and de-icing systems
mXcgi.Api = [Wt.Feq.ApiGD,-XLE_w - 0.25*cRootSub + cglocAC + x_engine + length_engine*.4]; % located at engine

% Oxygen System
% CG placed in center of habitable cabin 
mXcgi.Oxygen = [Wt.Feq.OxygenGD,-XLE_w - 0.25*cRootSub + cglocAC + (L_C/2) + L_N + L_CP]; 

% Auxiliary Power Unit
% Placed at Engines
mXcgi.Apu1 = [Wt.Feq.Apu,-XLE_w - 0.25*cRootSub + cglocAC + x_engine + length_engine*.4];

% Aircraft Furnishings (Seats, Lavatory, etc.)
% CG placed in center of habitable cabin
mXcgi.Furn =[Wt.Feq.Furn,-XLE_w - 0.25*cRootSub + cglocAC + (L_C/2) + L_N + L_CP]; 

% Operational Weight( Food, Drinks, ETC)
mXcgi.Oper = [Wt.Feq.Oper,-XLE_w - 0.25*cRootSub + cglocAC + L_N + L_CP]; % CG placed in front of cabin

% Paint Weight
mXcgi.Paint = [Wt.Feq.Paint,-XLE_w - 0.25*cRootSub + cglocAC + 0.40*L_F]; % Fuselage CG paint
=======
mXcgi.HydraulicNose = [Wt.Feq.Hydraulic*Nose_rat,-XLE_w + cglocAC + x_NoseGear];
mXcgi.HydraulicMain = [Wt.Feq.Hydraulic*Main_rat,-XLE_w + cglocAC + x_MainGear];

% Electrical System
mXcgi.ElectricalSystem = [Wt.Feq.Iae,-XLE_w + cglocAC + x_cockpit]; % Under Cockpit

% Air Conditioning, Pressurazation, and de-icing systems
mXcgi.Api = [Wt.Feq.ApiGD,-XLE_w + cglocAC + x_engine + length_engine*.4]; % located at engine

% Oxygen System
mXcgi.Oxygen = [Wt.Feq.OxygenGD,-XLE_w + cglocAC + (53.33+25)/2]; % CG placed in center of habitable cabin 

% Auxiliary Power Unit
% Placed at Engines
% mXcgi.Apu2 = [Wt.Feq.Apu*(2/3),-XLE_w + cglocAC + x_engine + length_engine*.4];
mXcgi.Apu1 = [Wt.Feq.Apu,-XLE_w + cglocAC + x_engine + length_engine*.4];

%Aircraft Furnishings
mXcgi.Furn =[Wt.Feq.Furn,(-XLE_w + cglocAC +53.33+25)/2]; % CG placed in center of habitable cabin

% Operational Weight( Food, Drinks, ETC)
mXcgi.Oper = [Wt.Feq.Oper,-XLE_w + cglocAC + 26]; % CG placed in front of cabin

% Paint Weight
mXcgi.Paint = [Wt.Feq.Paint,-XLE_w + cglocAC + 0.40*L_F]; % Fuselage CG paint
>>>>>>> Stashed changes

% Empty Weight CG
CG_types = fieldnames(mXcgi);
moment = 0;
wt = 0;
<<<<<<< Updated upstream

for i = 1:numel(CG_types)
   moment = moment - mXcgi.(CG_types{i})(1) .* mXcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mXcgi.(CG_types{i})(1);
end

wtEmpty = wt; 
cgEmpty = moment/wt;

% If positive - nose heavy
% If negative - tail heavy
fprintf('The empty weight Xcg is %0.4f \n', cgEmpty)
fprintf('The Difference Between Required and Actual XCG %0.4f ft \n', cgEmpty);

%% Empty Weight + Crew CG (OEW CG)
mXcgi.Crew = [Wt.oew.crew,-XLE_w - 0.25*cRootSub + cglocAC + x_cockpit]; % crew in cockpit 
=======
 for i = 1:numel(CG_types)
   moment = moment - mXcgi.(CG_types{i})(1) .* mXcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mXcgi.(CG_types{i})(1);
 end
wtEmpty = wt; 
cgEmpty = moment/wt;
fprintf('The empty weight cg is %0.4f \n', cgEmpty)
fprintf('The Difference Between Required and Actual CG %0.4f ft \n', cgEmpty);

%% Empty Weight + Crew CG (OEW CG)
mXcgi.Crew = [Wt.oew.crew,-XLE_w + cglocAC + x_cockpit]; % crew in cockpit 
>>>>>>> Stashed changes
CG_types = fieldnames(mXcgi);
moment = 0;
wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mXcgi.(CG_types{i})(1) .* mXcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mXcgi.(CG_types{i})(1);
 end
wtOEW = wt; 
cgOEW = moment/wt;
<<<<<<< Updated upstream
fprintf('The operating empty weight Xcg is %0.4f \n', cgOEW)
%% Empty Weight + Crew + Fuel CG

wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mXcgi.(CG_types{i})(1) .* mXcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mXcgi.(CG_types{i})(1);
 end
wtFuel = wt + Wt.fuel.w_tot; 
cgFuel = moment/wt;
fprintf('The Xoew + fuel cg is %0.4f \n', cgFuel)
%% Empty Weight + Crew + PAX (No baggage) CG

% Assume passenger distribution places cg at cabin center
mXcgi.pldNoBag = [Wt.pld.n_pass * Wt.pld.apw, -XLE_w - 0.25*cRootSub + cglocAC + (L_C/2) + L_N + L_CP]; 
=======
fprintf('The operating empty weight cg is %0.4f \n', cgOEW)
%% Empty Weight + Crew + Fuel CG


% How are we distributing the fuel again?

%% Empty Weight + Crew + PAX (No baggage) CG

% Assume passenger distribution places cg at cabin center
mXcgi.pldNoBag = [Wt.pld.n_pass * Wt.pld.apw, -XLE_w + cglocAC + (53.33+25)/2]; 
>>>>>>> Stashed changes
CG_types = fieldnames(mXcgi);
moment = 0;
wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mXcgi.(CG_types{i})(1) .* mXcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mXcgi.(CG_types{i})(1);
 end
wtPAX_NoBag = wt; 
cgPAX_NoBag = moment/wt;
<<<<<<< Updated upstream
fprintf('The Xcg for OEW + PAX (No baggage) is %0.4f \n', cgPAX_NoBag)
%% Empty Weight + Crew + PAX (w/ baggage) CG
x_baggage = 4; % distance of baggage from end of cabin
% Assume passenger distribution places cg at cabin center
mXcgi.pldBag = [Wt.pld.lug * Wt.pld.n_pass, -XLE_w - 0.25*cRootSub + cglocAC + (L_C/2) + L_N + L_CP + x_baggage]; 
=======
fprintf('The cg for OEW + PAX (No baggage) is %0.4f \n', cgPAX_NoBag)
%% Empty Weight + Crew + PAX (w/ baggage) CG

% Assume passenger distribution places cg at cabin center
mXcgi.pldBag = [Wt.pld.w_tot, -XLE_w + cglocAC + (53.33+25)/2]; 
>>>>>>> Stashed changes
CG_types = fieldnames(mXcgi);
moment = 0;
wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mXcgi.(CG_types{i})(1) .* mXcgi.(CG_types{i})(2); % Convert sign NEGATIVE = FLIP OVER
   wt = wt + mXcgi.(CG_types{i})(1);
 end
wtPAX_Bag = wt; 
cgPAX_Bag = moment/wt;
<<<<<<< Updated upstream
fprintf('The Xcg for OEW + PAX (with baggage) is %0.4f \n', cgPAX_Bag)

%% Empty Weight + Crew + PAX (w/ baggage) + Fuel (MTOW) CG

 for i = 1:numel(CG_types)
   moment = moment - mXcgi.(CG_types{i})(1) .* mXcgi.(CG_types{i})(2); % Convert sign NEGATIVE = FLIP OVER
   wt = wt + mXcgi.(CG_types{i})(1);
 end
wtMTOW = wt + Wt.fuel.w_tot; 
cgMTOW = moment/wt;
fprintf('The Xcg for MTOW is %0.4f \n', cgMTOW)

%% CG Excursion Diagram

wtArray = [wtEmpty; wtOEW; wtPAX_NoBag; wtPAX_Bag; wtFuel; wtMTOW];
cgArray = [cgEmpty; cgOEW; cgPAX_NoBag; cgPAX_Bag; cgFuel; cgMTOW]; 
cgRange = max(cgArray) - min(cgArray); 
fprintf('The cg range is %0.2f ft \n', cgRange)
figure()
title('Aircraft CG Excursion')
xlabel('CG Location (F.S.), ft')
ylabel('Weight (lbs)')
hold on;

% Plot CG Diagram
plot(cgArray(1), wtArray(1), 'ro', 'MarkerSize',10)
plot(cgArray(2), wtArray(2), 'go', 'MarkerSize',10)
plot(cgArray(3), wtArray(3), 'o', 'MarkerSize',10, 'Color',[0.5 0 1])
plot(cgArray(4), wtArray(4), 'o', 'MarkerSize',10, 'Color',[0 0.4 0.3])
plot(cgArray(5), wtArray(5), 'o', 'MarkerSize',10, 'Color',[1 0.5 0])
plot(cgArray(6), wtArray(6), 'o', 'MarkerSize',10, 'Color',[0 1 1])
plot(cgArray,wtArray, 'b--') % plot line through points
legend('Empty Wt CG','OEW Wt CG','OEW Wt + Fuel CG', ...
    'OEW Wt + Fuel + PAX (No Baggage) CG', ... 
    'OEW Wt + Fuel + PAX (With Baggage) CG',...
    'MTOW CG','Location','Northwest')
legend('boxoff')
=======
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
>>>>>>> Stashed changes
