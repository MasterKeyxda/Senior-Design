
if ~exist('ctrl', 'var')
    clc;
    clear;
    load('aircraft_vars.mat')
    close all;
end
%% Structural Weight

% Structural Weight includes the weight of the wing , empennage, fuselage,
% nacelles, and landing gear

% Wing
% Torenbeek Method, Eqn 5.7, p.69
WF = Wt.fuel.w_tot; % weight of fuel (lbs)
Wmzf = Wt.WTO - WF; % max zero fuel weight (lbs)
nUlt = 3.75; % ultimate load factor from V-n diagram % GET FROM ANOTHER SCRIPT
percentSub = 0.20; % percentage of subsonic portion of wing
percentSuper = 0.80; % percentage of supersonic portion of wing
SwSub = S_w * percentSub; % subsonic wing area (ft^2) 
SwSuper = S_w * percentSuper; % supersonic wing area (ft^2)

% Subsonic portion of wing (dim from Solidworks model Wing.v2)
cRootSub = 23.475; % root chord length (ft); from Solidworks model
thickToChordSub = 0.10; % thickness to chord ratio subsonic wing
trSub = cRootSub*thickToChordSub; % max thickness of wing root chord; subsonic wing (ft)
sweepSemiChordSub = 60; % wing semi-chord sweep angle (degrees)
bSub = percentSub * WING.geom.span; % subsonic wing span (ft)
% Subsonic wing weight (lbs)
Wt.Struc.WingSub = 0.0017*Wmzf*((bSub/cosd(sweepSemiChordSub))^0.75)*((1 + sqrt((6.3*cosd(sweepSemiChordSub))/bSub)))*(nUlt^0.55)*(((bSub*SwSub) / (trSub*Wmzf*cosd(sweepSemiChordSub)))^0.30);

% Supersonic portion of wing (dim from Solidworks model Wing.v2)
cRootSuper = 19.364; % root chord supersonic wing; 
thickToChordSuper = 0.05; % thickness to chord ration supersonic wing
trSuper = cRootSuper * thickToChordSuper; % max thickness of wing root chord; supersonic wing (ft)
sweepSemiChordSuper = 10; % wing semi-chord sweep angle (degrees)
bSuper = percentSuper * WING.geom.span; % supersonic wing span (ft)
% Supersonic wing weight (lbs)
Wt.Struc.WingSuper = 0.0017*Wmzf*((bSuper/cosd(sweepSemiChordSuper))^0.75)*((1 + sqrt((6.3*cosd(sweepSemiChordSuper)/bSuper))))*(nUlt^0.55)*(((bSuper*SwSuper) / (trSuper*Wmzf*cosd(sweepSemiChordSuper)))^0.30);

% Add Subsonic + Supersonic Wings together
Wt.Struc.Wing = Wt.Struc.WingSub + Wt.Struc.WingSuper; 

% Adjust wing weight for addition of flaps
Kflaps = 1.03; % assume 3% increase in weight
Wt.Struc.Wing = Wt.Struc.Wing * Kflaps; 

% Horizontal Tail
% Torenbeek Method, Eqn 5.19, p.74
semiChordSweepHT = 10; % semi chord sweep angle (degrees); per Keyur
Kh = 1.1; % constant for variable incidence stabilizers (HT)
Wt.Struc.HT = Kh*TAIL.Sh*(3.81*(((TAIL.Sh^0.2)*VD)/(1000*sqrt(cosd(semiChordSweepHT))))-0.287); 

% Vertical Tail
% Torenbeek Method, Eqn(s) 5.20 & 5.21, p.74
semiChordSweepVT = 10; % VT semi chord sweep angle (degrees); per Keyur
bv = TAIL.bv;%9.3846; % tail span; % WON'T UPDATE FROM TAIL.m
zh = bv; % distance from VT root to HT mounting location (ft);
% assume tail mounted at top of VT
Kv = 1 + (0.15*((TAIL.Sh*zh) / (TAIL.Sv*bv))); 
% Wt.Struc.VT = Kv*TAIL.Sv*(3.81*(((TAIL.Sv^0.2)*VD)/(1000*sqrt(cosd(semiChordSweepVT))))-0.287); 
Wt.Struc.VT = 900; % TEMPORARY

% Fuselage
% Sadraey, Eqn 10.7 p. 562
kpf = 0.0025; % Table 10.11 fuselage density factor, transport
densFuse = 2711 * 0.062428; % Table 10.6 density of aerospace aluminum; (lb/ft^3)
Kinlet = 1.25; % constant for inlets on fuselage; 1 if elsewhere
Wt.Struc.Fuselage = L_F * (Dmax^2) * densFuse * kpf * (nUlt^0.25) * Kinlet; 

% Nacelles
% Torenbeek Method, Eqn 5.36, p.80
% For turbojet engine; based on req. take-off thrust
Wt.Struc.Nacelle = 0.055*constraints.req_Thr; % accounts for weight of all nacelles for 3 engines

% Torenbeek Method, Eqn 5.42, p.82
% Applies to business jets with main gear mounted to wing and nose gear on
% the fuselage. For retractable LG
Kgr = 1; % constant for low wing airplanes

% Nose Gear 
AgNose = 12; % constant (Table 5.1)
BgNose = 0.06; % constant (Table 5.1)
Wt.Struc.NoseGear = Kgr*(AgNose + (BgNose*(Wt.WTO^0.75)));

% Main Gear
% Constants (Table 5.1)
AgMain = 33;
BgMain = 0.04;
CgMain = 0.021; 
Wt.Struc.MainGear = Kgr*(AgMain + (BgMain*(Wt.WTO^0.75)) + CgMain*Wt.WTO);

% Total Weight
Wt.Struc.Total = Wt.Struc.Wing + Wt.Struc.HT + Wt.Struc.VT + ...
    Wt.Struc.Fuselage + Wt.Struc.Nacelle + Wt.Struc.NoseGear + ... 
    Wt.Struc.MainGear;
%% Powerplant Weight

% Powerplane weight includes the weight of the engines, fuel system,
% propulsion system, air induction system, starting and ignition system

% Engine
Ne = 3; % number of engines
if strcmp(Wt.enginetype.name, 'Pegasus')
    Wt.Pwr.Engine = Ne*3960; % Pegasus 11-61 dry engine weight (lbs)
elseif strcmp(Wt.enginetype.name, 'JT8D-219')
    Wt.Pwr.Engine = 4741 * Ne;
% JT8D-219 dry engine weight = 4741 lbs
% Air Induction System
elseif strcmp(Wt.enginetype.name, 'PW-TF33-3-7')
    Wt.Pwr.Engine = 4650 * Ne;
end


% Fuel System
Kfsp = 6.71; % Jet-A, specific weight of fuel, lbs/gal (Table 7.9 Sadraey)
% GD Method, Eqn 6.20 and 6.22, p.92
% For wing with self-sealing bladder cells
Wt.Pwr.FuelSystem = (41.6*(((WF/Kfsp)/100)^0.818)) + 7.91*(((WF/Kfsp)/100)^0.854);

% For wing with integral fuel tanks and wet wing
% Torenbeek Method, Eqn 6.24, p.92
% Nt = 3; % number of fuel tanks
% Wt.Pwr.FuelSystem = 80 * (Ne + Nt - 1) + ((15 * (Nt^0.5)) * ((WF / Kfsp)^0.333)); 

% Propulsion System 
% Includes weight of engine controls and weight of engine starting system

% Engine Controls
% GD Method, eqn 6.23, p.93
% Fuselage/wing-root mounted jet engines, 
Kec = 0.686; % non-afterburning engines
Wt.Pwr.EngineControls = Kec * ((L_F * Ne)^0.792); 

% Engine Starting System
% GD Method, Eqn(s) 6.27 & 6.28, p.94
% For airplane with two jet engines with pneumatic starting systems
EngStartSysTwo = 9.33 * ((Wt.Pwr.Engine/1000)^1.078);
% For airplane with four or more jet engines with pneumatic starting systems
EngStartSysFour = 49.19 * ((Wt.Pwr.Engine/1000)^0.541);
% Average results of both equations to apply for 3 engines
Wt.Pwr.EngStartSys = (EngStartSysTwo + EngStartSysFour)/2;

Wt.Pwr.Propulsion = Wt.Pwr.EngineControls + Wt.Pwr.EngStartSys;
% Propulsion weight is part of engine for electric power
% Total Powerplant Weight 
Wt.Pwr.Total = Wt.Pwr.Engine + Wt.Pwr.FuelSystem + Wt.Pwr.Propulsion;
%% Fixed Equipment Weight

% Fixed Equipment weight includes the weight of the flight control system
% hydraulic and pneumatics, electrical system, instrumentation, avionics, 
% and electronics, air conditioning, pressurization, and de-icing, 
% auxilliary power unit (APU), oxygen, furnishings, baggage and cargo
% handling, operational items, and paint

% Flight Control System
% GD Method, Eqn 7.5, p.99
VD = (VD / sqrt(sigma)) * 1.68781; % convert KEAS to KTAS to ft/s 
qD = 0.5*rho_cr*(VD^2); % design dynamic pressure (lbs/ft^2)
Wt.Feq.FCsysGD = 56.01 * (((Wt.WTO*qD) / 100000)^0.576); % flight control system

% Torenbeek Method, Eqn 7.6, p.99
Kfc = 0.64; % constant for airplane with powered flight controls
Wt.Feq.FCsysToren = Kfc * (Wt.WTO^(2/3)); % flight control system

% Hydraulic and Pneumatic Systems
% p.101
Wt.Feq.Hydraulic = 0.0130 * Wt.WTO; % ranges from 0.0070 to 0.0150 of WTO for
% business jets

% Instrumentation, Avionics, Electronics
% GD Method, Eqn 7.23, p.103
Npil = 2; % number of pilots
Wt.Feq.Iae = Npil*(15 + 0.032*(Wt.WTO/1000)) + Ne*(5 + 0.006*(Wt.WTO/1000)) + 0.015*(Wt.WTO/1000) + (0.012*Wt.WTO); 

% Torenbeek Method, Eqn 7.25, p.104
% WE = 60983; % empty weight
% Wt.Feq.IaeToren = 0.575 * (WE^0.556)*(req.range^0.25);

% Electical System
% GD Method, Eqn 7.15, p.102
Wt.Feq.ElecSys = 1163 * (((Wt.Pwr.FuelSystem + Wt.Feq.Iae) / 1000)^0.506);

% Air-conditioning, pressurization, de-icing systems
% GD Method, Eqn 7.29, p.105
Vpax = 609; % volume of passenger cabin (ft^3)
Wt.Feq.ApiGD = 469*((Vpax*(Wt.pld.n_pass + Wt.oew.n_crew)/10000)^0.419); 

% Torenbeek Method, Eqn 7.30 p.105
Wt.Feq.Api = 6.75*(L_C^1.28); %L_C = cabin length

% Oxygen System
% GD Method, Eqn 7.35, p.106
Wt.Feq.OxygenGD = 7 * ((Wt.pld.n_pass + Wt.oew.n_crew)^0.702);

% Auxiliary Power Unit (APU)
% Not sure if needed
% Eqn 7.40
Wt.Feq.Apu = 0.005*Wt.WTO; % can be estimated as 0.004 to 0.013WTO

% Airplane Furnishings 
% Furnishing includes:
% 1. Seats, insulations, trim panels, sound proofing, instrument panels,
% control stands, lighting and wiring
% 2. Pantry structure and provisions
% 3. Lavatory and associated systems
% 4. Overhead luggage contrainers, wardrobes
% 5. Escape provisions, fire fighting equipment
% Wt.Feq.FurnStand = (43.7 - 0.037*Wt.pld.n_pass)*Wt.pld.n_pass + 46*Wt.pld.n_pass; 
% Appendix A - Roskam Part 5 Data
% Learjet 25 --> 8 passenger, Wt.Furnishings = 720 lbs
% Learjet 28 --> 8 passenger, Wt.Furnishings = 768 lbs
% Lockheed Jetstar --> 8-10 passengers, Wt.Furnishings = 1521 lbs
% Cessna Citation --> 8 passengers, Wt.Furnishings = 800 lbs
Wt.Feq.Furn = 756; % chosen based on similar business jets

% Baggage and Cargo Handling Equipment
% GD Method, Eqn 7.48, p.110
Kbc = 0.316; % constant with preload provisions
Wt.Feq.CargoEquip = Kbc * (Wt.pld.n_pass^1.456);

% Operational Items
% Includes food, water, drinks, lavatory supplies, p.108
% Standford http://adg.stanford.edu/aa241/AircraftDesign.html
Wt.Feq.Oper = 28 * Wt.pld.n_pass; 

% Paint Weight
Wt.Feq.Paint = 0.005*Wt.WTO; % ranges from 0.003 to 0.006WTO

% Total Fixed Equipment Weight
Wt.Feq.Total = Wt.Feq.FCsysGD + Wt.Feq.Hydraulic + Wt.Feq.Iae + ...
    Wt.Feq.ElecSys + Wt.Feq.ApiGD + Wt.Feq.OxygenGD + Wt.Feq.Apu + ... 
    Wt.Feq.Furn + Wt.Feq.CargoEquip + Wt.Feq.Oper + Wt.Feq.Paint; 


%% Export to Excel
StrucLabels = {'Wing';'Horizontal Tail';'Vertical Tail';'Fuselage';...
    'Nacelles';'Nose Gear';'Main Gear';'Landing Gear';'Total Structural Weight'};
PwrLabels = {'Engines'; 'Fuel System'; 'Propulsion System'; 'Total Powerplant Weight'};
FeqLabels = {'Flight Control System'; 'Hydraulic Systems'; ...
    'Instrumentation, Avionics, Electronics';'Electical System'; ...
    'Air-conditioning, pressurization, de-icing systems';'Oxygen System'; ...
    'Auxiliary Power Unit (APU)'; 'Furnishings'; ... 
    'Baggage and Cargo Handling Equipment'; 'Operational Items'; 'Paint'; ...
    'Total Fixed Equipment Weight'}; 

% Structural Weight Array
StrucWT = [Wt.Struc.Wing; Wt.Struc.HT; Wt.Struc.VT; ...
    Wt.Struc.Fuselage; Wt.Struc.Nacelle; Wt.Struc.NoseGear; ...
    Wt.Struc.MainGear; Wt.Struc.Total];

% Powerplant Weight Array
PwrWT = [Wt.Pwr.Engine; Wt.Pwr.FuelSystem; Wt.Pwr.Propulsion; Wt.Pwr.Total];

% Fixed Equipment Weight Array
FeqWT = [Wt.Feq.FCsysToren; Wt.Feq.Hydraulic; Wt.Feq.Iae; Wt.Feq.ElecSys; ...
    Wt.Feq.Api; Wt.Feq.OxygenGD; Wt.Feq.Apu; Wt.Feq.Furn; Wt.Feq.CargoEquip; ... 
    Wt.Feq.Oper; Wt.Feq.Paint];

% Total Weight
%WETotal = Wt.Struc.Total + Wt.Pwr.Total + Wt.Feq.Total;
% Concatenate arrays
%WTArray = [StrucWT;PwrWT;FeqWT;WETotal];
%Labels = {StrucLabels(:); PwrLabels; FeqLabels}
%xlswrite('Aircraft_Weight.xlsx',WETotal, 'B3:B27')

% Class 2 Weight Summary
WEnew = Wt.Struc.Total + Wt.Pwr.Total + Wt.Feq.Total;
% WTOnew = WEnew + WF + Wt.pld.w_tot + Wt.oew.crew; 
WTOnew = (WEnew + Wt.pld.w_tot + Wt.oew.crew)/(1 - Wt.fuel.Wf_Wto);

% Percent difference between new and prelminary WTO
WeightDiff = (abs((WTOnew - Wt.WTO)) / Wt.WTO) * 100; 
if exist('ctrl', 'var')
    fprintf('\n\t WEIGHT CORRECTIONS \n');
    fprintf('The discrepancy between the preliminary WTO and the new WTO is: %0.2f percent \n', WeightDiff)
    if WeightDiff > 5
        fprintf('Iteration required. \n')
        Wt.WTO = WTOnew; % overwrite preliminary sizing MTOW (lb)
        Wt.fuel.w_tot = Wt.fuel.Wf_Wto * Wt.WTO;
        fprintf('New weight set to %6.2f lb\n', WTOnew);
    else
        fprintf('No iteration required. \n')
        Wt.WTO = WTOnew; % overwrite preliminary sizing MTOW (lb)
        Wt.fuel.w_tot = Wt.fuel.Wf_Wto * Wt.WTO;
        fprintf('Final weight set to %6.2f lb\n', WTOnew);
        if Wt.fuel.w_tot > Wt.fuel.w_max
            fprintf('Fuel mass exceeds max requirements by %0.5f percent\n', (Wt.fuel.w_tot - Wt.fuel.w_max)/Wt.fuel.w_tot * 100);
        else
            fprintf('Fuel mass meets max requirements by %0.5f percent\n', abs(Wt.fuel.w_tot - Wt.fuel.w_max)/Wt.fuel.w_tot * 100);
        end
    end
end
