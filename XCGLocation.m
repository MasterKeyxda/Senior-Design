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

% Hydraulic System for Landing Gear
Nose_rat = Wt.Struc.NoseGear/(Wt.Struc.NoseGear+Wt.Struc.MainGear); % weight ratio of nose landing gear
Main_rat = Wt.Struc.MainGear/(Wt.Struc.NoseGear+Wt.Struc.MainGear); % weight ratio of main landing gear
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

%% CG Excursion Diagram

filename = 'CG_Calc_Copy.xlsx';

% Read in data from 'CG_Calc_Copy.xlsx' spreadsheet
Wt.Excel = xlsread(filename,'G34:G39'); 
cgExcel = xlsread(filename,'I34:I39'); 
close all;
figure()
title('CG Excursion Diagram')
xlabel('CG Location (F.S.), ft')
ylabel('Weight (lbs)')
xMin = min(cgExcel-1);
xMax = max(cgExcel+1);
xlim([xMin, xMax])
hold on;

% Print range of most fwd and aft cg locations
cgRange = max(cgExcel) - min(cgExcel);
fprintf('The cg range is %0.2f ft \n', cgRange)

% Plot CG Diagram
plot(cgExcel, Wt.Excel,'b')
plot(cgExcel, Wt.Excel,'s', 'MarkerSize',5,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])

yMin = 35000; 
yMax = 100000;
ylim([yMin,yMax]) % set y-axis limits
 
% Annotate and plot horizontal line from xMin to cg at MTOW
plot([xMin, cgExcel(4)], [Wt.Excel(4), Wt.Excel(4)],'k--') 
str = 'MTOW';
text(cgExcel(4)-1.5,Wt.Excel(4)+2000,str)

% Annotate and plot horizontal line from xMin to cg at Empty Weight
plot([xMin, cgExcel(1)], [Wt.Excel(1), Wt.Excel(1)],'k--')
str = 'W_E';
text(cgExcel(4)-1.5,Wt.Excel(1)-2000,str)
