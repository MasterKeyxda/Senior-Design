%% Constraint Plots: Roskam

% This script determines an appropriate wing loading (W/S) and
% thrust-to-weight ratio (T/W) using the procedure outlined in Roskam
% Airplane Design Part 1

constraints.wingLoading = 40:0.5:200; % wing loading parameter (lb/ft^2)
% CL for takeoff and landing Roskam Part I, p.91; 
constraints.CLMax = 1.6; 
constraints.CLMax_L = 2; % landing
constraints.CLMax_TO = 2; % takeoff

%% Stall Speed (Section 3.1)

% Flaps UP
constraints.Vstall_Up = 125 * 1.688; % stall speed, knots -> ft/s

% Generate stall speed curve (Flaps Up)
constraints.Vstall_Crv_Up = 0.5*atm.rho_sl*(constraints.Vstall_Up^2)*constraints.CLMax;

% Flaps Down (Landing)
constraints.Vstall_Down = 115 * 1.688; % stall speed, knots -> ft/s


% Generate stall speed curve (Flaps Up)
constraints.Vstall_Crv_Down = 0.5*atm.rho_sl*(constraints.Vstall_Down^2)*constraints.CLMax_L;

% Choose min W/S from flaps up and down scenarios
constraints.Vstall_Crv = min(constraints.Vstall_Crv_Up, constraints.Vstall_Crv_Down);

%% Take-off Distance Sizing (Section 3.2.3 - 3.2.4)

% FAR 25 takeoff parameter (lbs/ft^2)
TOP_25 = req.takeoffRun / 37.5; % Eqn 3.8 Roskam Pt 1

atm.alt_TO = 8; % alt for takeoff sizing (kft)
[~,~,atm.sig_rhoTO,~] = AltTable(atm.alt_TO,'h'); % sigma at takeoff

% Takeoff field length (Eqn 3.7)
constraints.takeoff.FL = TOP_25 * atm.sig_rhoTO; 

% Generate takeoff distance curve
constraints.TO_Curve = constraints.wingLoading / (constraints.CLMax_TO * constraints.takeoff.FL);

%% Landing Distance Sizing (Section 3.3.3 - 3.3.4)

% Landing Field Length (Assume to be the same as takeoff)
req.landFL = 5000; % ft

constraints.Vaprch = sqrt(req.landFL/0.3); % approach speed (knots)
constraints.Vstall_L = (constraints.Vaprch / 1.3) * 1.688; % stall speed landing (ft/s)

% Wing Loading during Landing Eqn 3.1
constraints.wingLoading_Land = (atm.rho_sl*(constraints.Vstall_L^2)*constraints.CLMax_L) / 2;

% ERROR 

