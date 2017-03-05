
clc
clear all;
load('aircraft_vars.mat')
%% Structural Weight

% Structural Weight includes the weight of the wing , empennage, fuselage,
% nacelles, and landing gear

%% Powerplant Weight

% Powerplane weight includes the weight of the engines, fuel system,
% propulsion system, air induction system, starting and ignition system

% Engine
Ne = 3; % number of engines
Wt.Pwr.Engine = Ne*3960; % Pegasus 11-61 dry engine weight (lbs)

% Fuel System
% For wing with integral fuel tanks and wet wing
% Torenbeek Method, Eqn 6.24, p.92
Kfsp = 6.71; % Jet-A, specific weight of fuel, lbs/gal
WF = Wt.fuel.w_tot; % weight of fuel (lbs)
Nt = 3; % number of fuel tanks
Wt.Pwr.FuelSystem = 80 * (Ne + Nt - 1) + ((15 * (Nt^0.5)) * ((WF / Kfsp)^0.333)); 

% Propulsion System 
% Includes weight of engine controls and weight of engine starting system

% Engine Controls
% GD Method, eqn 6.23, p.93
% Fuselage/wing-root mounted jet engines, 
Kec = 0.686; % non-afterburning engines
lf = 140.8; % length of fuselage ft
Wt.Pwr.EngineControls = Kec * ((lf * Ne)^0.792); 

% Engine Starting System
% GD Method, Eqn(s) 6.27 & 6.28, p.94
% For airplane with two jet engines with pneumatic starting systems
EngStartSysTwo = 9.33 * ((Wt.Pwr.Engine/1000)^1.078);
% For airplane with four or more jet engines with pneumatic starting systems
EngStartSysFour = 49.19 * ((Wt.Pwr.Engine/1000)^0.541);
% Average results of both equations to apply for 3 engines
Wt.Pwr.EngStartSys = (EngStartSysTwo + EngStartSysFour)/2;

Wt.Pwr.Propulsion = Wt.Pwr.EngineControls + Wt.Pwr.EngStartSys;

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
atm.rho_SL = 0.002378; % SL density slug/ft^3
VD = 543.6; % design dive speed (KEAS), from V-n diagram
VD = VD * sqrt(atm.rho_SL / rho_cr); % convert to KTAS; 
VD = VD * 1.68781; % convert KTAS to ft/s
qD = 0.5*rho_cr*(VD^2); % dynamic pressure (lbs/ft^2)
Wt.Feq.FCsysGD = 56.01 * (((WTO*qD) / 100000)^0.576); % flight control system

% Torenbeek Method, Eqn 7.6, p.99
Kfc = 0.64; % constant for airplane with powered flight controls
Wt.Feq.FCsysToren = Kfc * (WTO^(2/3)); % flight control system

% Hydraulic and Pneumatic Systems
% p.101
Wt.Feq.Hydraulic = 0.0070 * WTO; % ranges from 0.0070 to 0.0150 of WTO for
% business jets

% Instrumentation, Avionics, Electronics
% GD Method, Eqn 7.23, p.103
Npil = 2; % number of pilots
Wt.Feq.Iae = Npil*(15 + 0.032*(WTO/1000)) + Ne*(5 + 0.006*(WTO/1000)) + 0.015*(WTO/1000) + (0.012*WTO); 

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
Wt.Feq.Apu = 0.005*WTO; % can be estimated as 0.004 to 0.013WTO

% Airplane Furnishings 
% Furnishing includes:
% 1. Seats, insulations, trim panels, sound proofing, instrument panels,
% control stands, lighting and wiring
% 2. Pantry structure and provisions
% 3. Lavatory and associated systems
% 4. Overhead luggage contrainers, wardrobes
% 5. Escape provisions, fire fighting equipment
% Appendix A - Roskam Part 5 Data
% Learjet 25 --> 8 passenger, Wt.Furnishings = 720 lbs
% Learjet 28 --> 8 passenger, Wt.Furnishings = 768 lbs
% Lockheed Jetstar --> 8-10 passengers, Wt.Furnishings = 1521 lbs
% Cessna Citation --> 8 passengers, Wt.Furnishings = 800 lbs
% Possibly reference Sadraey
Wt.Feq.Furn = 756;

% Baggage and Cargo Handling Equipment
Kbc = 0.0646; % constant without preload provisions
Wt.Feq.CargoEquip = Kbc * (Wt.pld.n_pass^1.456);

% Operational Items
% Includes food, water, drinks, lavatory supplies, p.108
% GD Method, lav, water, food prov portions of Eqn 7.44
% Klav = 3.90; % constant for business jets
% Kbuf = 1.02; 
% Wt.Feq.Oper = Klav*(Wt.pld.n_pass^1.33) + Kbuf*(Wt.pld.n_pass^1.12);
Wt.Feq.Oper = 17 * Wt.pld.n_pass; % AircraftStructures.pdf 
% (see Weight and Balance folder - OneDrive)

% Paint Weight
Wt.Feq.Paint = 0.003*WTO; % ranges from 0.003 to 0.006WTO

% Total Fixed Equipment Weight
Wt.Feq.Total = Wt.Feq.FCsysToren + Wt.Feq.Hydraulic + Wt.Feq.Iae + ...
    Wt.Feq.ElecSys + Wt.Feq.Api + Wt.Feq.OxygenGD + Wt.Feq.Apu + ... 
    Wt.Feq.Furn + Wt.Feq.CargoEquip + Wt.Feq.Oper + Wt.Feq.Paint; 