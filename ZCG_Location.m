%% ZCG Location
clear CG_types
clear mZcgi

% All references and moment arms in ft 
% Zcg calc performed with reference to bottom of fuselage
% Lifting Components - referenced to desired CG
% Wing
hWing = 0.5; % distance from bottom of fuselage to bottom of wing
% Assume wing zcg is at half of wing thickness
mZcgi.Wing = [Wt.Struc.Wing, 0.5*trSub + 0.5 + hWing];

% Horizontal Tail (NOTE: need to account for taper of fuselage -- ASK MATT)
TAIL.thickRatio = 0.05; % Tail thickness ratio
mZcgi.HTail = [Wt.Struc.HT, Dmax + TAIL.bv + 0.5*TAIL.thickRatio*TAIL.ch]; 

% Vertical Tail
% need moment arm
% mZcgi.VTail = [Wt.Struc.VT]

% Fuselage
mZcgi.Fuselage = [Wt.Struc.Fuselage, Dmax/2];

% Engine
D_eng = 54/12; % engine diameter (ft) P&W TF33-P-7
% Side engines placed slightly above halfway point of fuselage
mZcgi.Engine2 = [Wt.Pwr.Engine*(2/3), (Dmax/2) + (0.25*D_eng)]; 

% Aft engine placed to at back of fuselage
mZcgi.Engine1 = [Wt.Pwr.Engine*(1/3), (Dmax/2) + (0.5 * D_eng)]; 

% Nacelle
% http://adg.stanford.edu/aa241/AircraftDesign.html Sect 9.2.2)
nacLength = (2.4077*(constraints.req_Thr^0.3876))/12; % ft
nacDMax = 1.0827*(constraints.req_Thr^0.4134)/12;
mZcgi.Nacelle2 = [Wt.Struc.Nacelle*(2/3),-XLE_w - 0.25*cRootSub + cglocAC + x_engine + nacLength*.38];
mZcgi.Nacelle1 = [Wt.Struc.Nacelle*(1/3),-XLE_w - 0.25*cRootSub + cglocAC + x_engine + 1.2*length_engine + nacLength*.38];
% Nose and Main Gear
x_NoseGear = 30; % ft nosegear distance from nose
x_MainGear = 85.92; % ft main gear distance from nose
mZcgi.NoseGear = [Wt.Struc.NoseGear,-XLE_w - 0.25*cRootSub + cglocAC + x_NoseGear];
mZcgi.MainGear = [Wt.Struc.MainGear,-XLE_w - 0.25*cRootSub + cglocAC + x_MainGear];
% Fuel System
% Fuel tank should have minimal effect on CG loaded and empty
LengthFuelTank = 12; % ft
% Put fuel tank CG on CG
mZcgi.FuelSystem = [Wt.Pwr.FuelSystem,0]; 
% Propulsion System
%mXcgi.=[Wt.Struc,-XLE_w + cglocAC + x_engine + length_engine*.4];

% Also Omit Starting Systems

% Fixed Equipment Weight
% Flight Control System
x_cockpit = 25 + 5.5/3; % center of cockpit from nose
mZcgi.FCsys = [Wt.Feq.FCsysGD,-XLE_w - 0.25*cRootSub + cglocAC + x_cockpit]; % In Cockpit

% Hydraulic System for Landing Gear
Nose_rat = Wt.Struc.NoseGear/(Wt.Struc.NoseGear+Wt.Struc.MainGear); % weight ratio of nose landing gear
Main_rat = Wt.Struc.MainGear/(Wt.Struc.NoseGear+Wt.Struc.MainGear); % weight ratio of main landing gear
mZcgi.HydraulicNose = [Wt.Feq.Hydraulic*Nose_rat,-XLE_w  - 0.25*cRootSub + cglocAC + x_NoseGear];
mZcgi.HydraulicMain = [Wt.Feq.Hydraulic*Main_rat,-XLE_w  - 0.25*cRootSub + cglocAC + x_MainGear];

% Instrumentation, Avionics
mZcgi.Instrumentation = [Wt.Feq.Iae,-XLE_w - 0.25*cRootSub + cglocAC + x_cockpit]; % Under Cockpit

% Electrical System 

% Air Conditioning, Pressurazation, and de-icing systems
mZcgi.Api = [Wt.Feq.ApiGD,-XLE_w - 0.25*cRootSub + cglocAC + x_engine + length_engine*.4]; % located at engine

% Oxygen System
mZcgi.Oxygen = [Wt.Feq.OxygenGD,-XLE_w - 0.25*cRootSub + cglocAC + (53.33+25)/2]; % CG placed in center of habitable cabin 

% Auxiliary Power Unit
% Placed at Engines
% mXcgi.Apu2 = [Wt.Feq.Apu*(2/3),-XLE_w + cglocAC + x_engine + length_engine*.4];
mZcgi.Apu1 = [Wt.Feq.Apu,-XLE_w - 0.25*cRootSub + cglocAC + x_engine + length_engine*.4];

%Aircraft Furnishings
mZcgi.Furn =[Wt.Feq.Furn,(-XLE_w - 0.25*cRootSub + cglocAC +53.33+25)/2]; % CG placed in center of habitable cabin

% Operational Weight( Food, Drinks, ETC)
mZcgi.Oper = [Wt.Feq.Oper,-XLE_w - 0.25*cRootSub + cglocAC + 26]; % CG placed in front of cabin

% Paint Weight
mZcgi.Paint = [Wt.Feq.Paint,-XLE_w - 0.25*cRootSub + cglocAC + 0.40*L_F]; % Fuselage CG paint

% Empty Weight CG
CG_types = fieldnames(mZcgi);
moment = 0;
wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mZcgi.(CG_types{i})(1) .* mZcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mZcgi.(CG_types{i})(1);
 end
wtEmpty = wt; 
cgEmpty = moment/wt;

% If positive - nose heavy
% If negative - tail heavy
fprintf('The empty weight cg is %0.4f \n', cgEmpty)
fprintf('The Difference Between Required and Actual CG %0.4f ft \n', cgEmpty);

%% Empty Weight + Crew CG (OEW CG)
mZcgi.Crew = [Wt.oew.crew,-XLE_w - 0.25*cRootSub + cglocAC + x_cockpit]; % crew in cockpit 
CG_types = fieldnames(mZcgi);
moment = 0;
wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mZcgi.(CG_types{i})(1) .* mZcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mZcgi.(CG_types{i})(1);
 end
wtOEW = wt; 
cgOEW = moment/wt;
fprintf('The operating empty weight cg is %0.4f \n', cgOEW)
%% Empty Weight + Crew + Fuel CG


% How are we distributing the fuel again?

%% Empty Weight + Crew + PAX (No baggage) CG

% Assume passenger distribution places cg at cabin center
mZcgi.pldNoBag = [Wt.pld.n_pass * Wt.pld.apw, -XLE_w - 0.25*cRootSub + cglocAC + (53.33+25)/2]; 
CG_types = fieldnames(mZcgi);
moment = 0;
wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mZcgi.(CG_types{i})(1) .* mZcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mZcgi.(CG_types{i})(1);
 end
wtPAX_NoBag = wt; 
cgPAX_NoBag = moment/wt;
fprintf('The cg for OEW + PAX (No baggage) is %0.4f \n', cgPAX_NoBag)
%% Empty Weight + Crew + PAX (w/ baggage) CG

% Assume passenger distribution places cg at cabin center
mZcgi.pldBag = [Wt.pld.w_tot, -XLE_w - 0.25*cRootSub + cglocAC + (53.33+25)/2]; 
CG_types = fieldnames(mZcgi);
moment = 0;
wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mZcgi.(CG_types{i})(1) .* mZcgi.(CG_types{i})(2); % Convert sign NEGATIVE = FLIP OVER
   wt = wt + mZcgi.(CG_types{i})(1);
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