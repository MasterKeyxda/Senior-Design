%% Constraint Plots: Roskam

% This script determines an appropriate wing loading (W/S) and
% thrust-to-weight ratio (T/W) using the procedure outlined in Roskam
% Airplane Design Part 1

constraints.wingLoading = 40:0.5:200; % wing loading parameter (lb/ft^2)
% CL for takeoff and landing Roskam Part I, p.91; 
constraints.CLMax = 1.6; 
constraints.CLMax_L = 2; % landing
constraints.CLMax_TO = 1.8; % takeoff

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

atm.alt_TO = 0; % alt for takeoff sizing (kft)
[~,~,atm.sig_rhoTO,~] = AltTable(atm.alt_TO,'h'); % sigma at takeoff

% Takeoff field length (Eqn 3.7)
constraints.takeoff.FL = TOP_25 * atm.sig_rhoTO; 

% Generate takeoff distance curve
constraints.TO_Curve = constraints.wingLoading / (constraints.CLMax_TO * constraints.takeoff.FL);

%% Landing Distance Sizing (Section 3.3.3 - 3.3.4)

% Landing Field Length (Assume to be the same as takeoff)
req.landFL = 7000; % ft

constraints.Vaprch = sqrt(req.landFL/0.3); % approach speed (knots)
constraints.Vstall_L = (constraints.Vaprch / 1.3) * 1.688; % stall speed landing (ft/s)

% Wing Loading during Landing Eqn 3.1
constraints.wingLoading_Land = (atm.rho_sl*(constraints.Vstall_L^2)*constraints.CLMax_L) / 2;

% Set Landing Weight Ratio based on Table 3.3 in Roskam
% Avg values for business jet and supersonic transports
Wt.LandRatio = 0.815; % Ratio of Landing Weight to WTO

% Relate wing loading during landing to takeoff
constraints.wingLoading_TO = constraints.wingLoading_Land / Wt.LandRatio;

%% Climb Sizing (Section 3.4)

%% Ceiling Sizing (Section 3.4.10 & p.183)

% Cruise air properties
atm.alt_cr = 46; % cruise alt in kft; adjust to match main script
[~,~,atm.sig_rho_cr,atm.asnd_cr] = AltTable(atm.alt_cr,'h'); 
atm.rho_cr = atm.sig_rho_cr * atm.rho_sl;

% Parameters at service ceiling
constraints.Vcruise = atm.asnd_cr * req.cr_M0(1) * 1116; % ft/s
constraints.ROC = 500/60; % Rate of Climb at service ceiling
constraints.L_D = 7; % assumed L/D max

% Required thrust at service ceiling
constraints.thrstService = (constraints.ROC/constraints.Vcruise) + (1/constraints.L_D);

% Convert thrust at service ceiling to takeoff thrust
constraints.thrstRatio = 0.23; % Estimate; need ratio of SL thrust to cruise thrust
constraints.ceilingCurve = constraints.thrstService/constraints.thrstRatio;

%% Cruise Speed Sizing (Section 3.6.4 - 3.6.5)

% Zero Lift Drag Coefficient
CD_0 = 0.0200; % Estimate

AR = 3; % aspect ratio 
e = 0.8; % Oswald efficiency factor

% Dynamic Pressure at specified cruise alt (psf)
atm.dyn_press = 0.5*(atm.rho_cr)*(constraints.Vcruise^2); 

% NOTE: Not required takeoff thrust (Eqn 3.60)
constraints.thrstCruise = ((CD_0*atm.dyn_press)./constraints.wingLoading) + (constraints.wingLoading ./ (atm.dyn_press * pi * AR * e)); 

% Takeoff thrust for Vmax Curve (p.182)
constraints.VmaxCurve = constraints.thrstCruise / constraints.thrstRatio;

%% Generate Constraint Plot

figure();
hold on; 
title('Wing Area (S) and Engine Thrust (T) Constraint Plot')
xlabel('Wing Loading, W/S (lbf / ft^2)')
ylabel('Thrust-to-Weight Ratio, T/W (lbf/lbf)')

% Landing Distance / Stall Curve
plot([constraints.wingLoading_TO,constraints.wingLoading_TO],[0,1.5])

% Takeoff Distance Curve 
plot(constraints.wingLoading, constraints.TO_Curve)

% Climb Curve


% Service Ceiling Curve
plot([constraints.wingLoading(1), constraints.wingLoading(end)],[constraints.ceilingCurve, constraints.ceilingCurve])

% Vmax Curve
plot(constraints.wingLoading, constraints.VmaxCurve)

legend('Landing/Stall', 'Takeoff', 'Service Ceiling','Cruise')
