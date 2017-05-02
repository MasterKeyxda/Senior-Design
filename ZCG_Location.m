%% ZCG Location
clear CG_types
clear mZcgi

% All references and moment arms in ft 
% Zcg calc performed with reference to bottom of fuselage
% Lifting Components - referenced to desired CG
% Wing

% Assume wing Zcg is at half of wing thickness
mZcgi.Wing = [Wt.Struc.Wing, 0.5*trSub];

% Horizontal Tail (NOTE: need to account for taper of fuselage -- ASK MATT)
TAIL.thickRatio = 0.05; % Tail thickness ratio
mZcgi.HTail = [Wt.Struc.HT, Dmax + TAIL.bv + 0.5*TAIL.thickRatio*TAIL.ch]; 

% Vertical Tail
% Assume VT zcg is 1/3 of its span 
mZcgi.VTail = [Wt.Struc.VT, Dmax + 0.33*bv];

% Fuselage
% Zcg located halfway point of fuselage
mZcgi.Fuselage = [Wt.Struc.Fuselage, Dmax/2];

% Engines
D_eng = 54/12; % engine diameter (ft) P&W TF33-P-7
% Side engines placed slightly above halfway point of fuselage
mZcgi.Engine2 = [Wt.Pwr.Engine*(2/3), (Dmax/2) + (0.25 * D_eng)]; 
% Aft engine placed at back of fuselage
mZcgi.Engine1 = [Wt.Pwr.Engine*(1/3), (Dmax/2) + (0.5 * D_eng)]; 

% Nacelle
% http://adg.stanford.edu/aa241/AircraftDesign.html Sect 9.2.2)
nacLength = (2.4077*(constraints.req_Thr^0.3876))/12; % ft
nacDMax = 1.0827*(constraints.req_Thr^0.4134)/12;
mZcgi.Nacelle2 = [Wt.Struc.Nacelle*(2/3), (Dmax/2) + (0.25 * D_eng)];
mZcgi.Nacelle1 = [Wt.Struc.Nacelle*(1/3), (Dmax/2) + (0.5 * D_eng)];

% Nose and Main Gear
% Need estimated strut length (Zcg is 0.50 strut length)
% Reference from bottom of fuselage 
strutLength = 5.84; % Estimated strut length (UPDATE LATER)
mZcgi.NoseGear = [Wt.Struc.NoseGear, -0.6 * strutLength];
mZcgi.MainGear = [Wt.Struc.MainGear, -0.6 * strutLength];

% Fuel System
% Fuel tank should have minimal effect on CG loaded and empty
LengthFuelTank = 12; % ft
% Assume Fuel System at 0.66 Zcg 
mZcgi.FuelSystem = [Wt.Pwr.FuelSystem, Dmax*0.5]; 

% Fixed Equipment Weight
% Flight Control System
z_cockpit = Dmax*0.65; % Based on CAD model; should be roughly 6 ft 
mZcgi.FCsys = [Wt.Feq.FCsysGD, z_cockpit];

% Hydraulic System for Landing Gear
Nose_rat = Wt.Struc.NoseGear/(Wt.Struc.NoseGear+Wt.Struc.MainGear); % weight ratio of nose landing gear
Main_rat = Wt.Struc.MainGear/(Wt.Struc.NoseGear+Wt.Struc.MainGear); % weight ratio of main landing gear
% Slightly above the bottom of fuselage
mZcgi.HydraulicNose = [Wt.Feq.Hydraulic*Nose_rat,0.75]; 
mZcgi.HydraulicMain = [Wt.Feq.Hydraulic*Main_rat,0.75];

% Instrumentation, Avionics
mZcgi.Instrumentation = [Wt.Feq.Iae,Dmax*0.75]; % In cockpit (Roughly 7 ft)

% Electrical System 
mZcgi.ElecSys = [Wt.Feq.ElecSys,3.5]; % Under cockpit

% Air Conditioning, Pressurization, and de-icing systems
mZcgi.Api = [Wt.Feq.ApiGD,4.5]; % located at engine

% Oxygen System
mZcgi.Oxygen = [Wt.Feq.OxygenGD,7]; % CG placed in center of habitable cabin 

% Auxiliary Power Unit
% Placed at Engines
% mXcgi.Apu2 = [Wt.Feq.Apu*(2/3),-XLE_w + cglocAC + x_engine + length_engine*.4];
mZcgi.Apu1 = [Wt.Feq.Apu,(Dmax/2) + (0.25 * D_eng)];

% Aircraft Furnishings
mZcgi.Furn =[Wt.Feq.Furn,(Dmax/2)]; % CG placed in center of habitable cabin

% Operational Weight( Food, Drinks, ETC)
mZcgi.Oper = [Wt.Feq.Oper,(Dmax/2)]; % Center of cabin

% Paint Weight
mZcgi.Paint = [Wt.Feq.Paint,(Dmax/2)]; % Fuselage CG paint

% Empty Weight CG
CG_types = fieldnames(mZcgi);
moment = 0;
wt = 0;
 for i = 1:numel(CG_types)
   moment = moment - mZcgi.(CG_types{i})(1) .* mZcgi.(CG_types{i})(2); %Convert sign NEGATIVE = FLIP OVER
   wt = wt + mZcgi.(CG_types{i})(1);
 end
wtEmpty = wt; 
ZcgEmpty = moment/wt;

% With reference to bottom of fuselage (ft
fprintf('The empty weight Zcg is %0.2f ft from the bottom of the fuselage. \n', abs(ZcgEmpty))

