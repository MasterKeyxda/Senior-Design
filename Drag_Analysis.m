%% Drag Analysis

% This script provides an estimate for CD0
% Drag breakdown procedure is based on the Aircraft Performance Analysis
% by M. Sadraey Ch 3
% Date Created: 3/17/17

% Run iter_weights.m before running this script
%% Common Equation Parameters
Sref = WING.geom.S_area; % Wing planform area chosen as reference area

% Fuselage
% LF_DF refers to fuselage fineness ratio (from Preliminary_Fuselage_Design.m)
fLD = 1 + (60 / (LF_DF^3)) + 0.0025*(LF_DF);
Swet.fuse = 2773; % from Solidworks (ft^2)

% Wing
Swet.wing = WING.geom.S_area - 2*0.5*D_C*(WING.geom.sub.Cr +...
    (WING.geom.sub.Cr+(0.5*D_C)*((WING.geom.sub.Ct-WING.geom.sub.Cr)/WING.geom.sub.span)))/2; % wing wetted area correction for removed area inside fuselage
Swet.wing = 2*Swet.wing; % wing wetted area
wing_TCRatio = 0.0525; % wing max thickness to chord ratio
ftc_w = 1 + (2.7*wing_TCRatio) + 100*((wing_TCRatio)^4); % ftc_w is a function of thickness ratio
Cdmin.wing = min(WING.biconvex(2).Cd);

% Horizontal Tail
Swet.ht = 2*TAIL.Sh; % wetted area of horizontal tail
HT_TCRatio = 0.0525; % horizontal tail max thickness to chord ratio
ftc_ht = 1 + (2.7*HT_TCRatio) + 100*((HT_TCRatio)^4);
Cdmin.ht = min(WING.biconvex(2).Cd);

% Vertical Tail
Swet.vt = 2*TAIL.Sv; % wetted area of vertical tail
VT_TCRatio = 0.0525; % vertical tail max thickness to chord ratio
ftc_vt = 1 + (2.7*VT_TCRatio) + 100*((VT_TCRatio)^4);
Cdmin.vt = min(WING.biconvex(2).Cd); 

nacL_D = 2; % nacelle fineness ratio
fLD_nac = 1 + (60 / (nacL_D^3)) + 0.0025*(nacL_D);
nNacelle = 2;
nacLength = 16.79; % nacelle length ft (Solidworks/keyur)
Swet.nacelle = 47.09; % from Solidworks (ft^2)

%% Supersonic Cruise

h = 42; % kft
M = 1.6; % Mach Number
visc = 0.297E-6; % viscosity (slugs/ft-s)

% Call Function to calculates dynamic pressure, reynolds number, and skin
% friction coefficients
[V, rho, q ,Cf_Turb,Cf_lam,Re_Ratio] = dragCalc(h, M, visc, L_F, WING.geom.MAC, TAIL.ch, TAIL.cv, nacLength);

%-----Fuselage-----%
% Eqn 3.14
fM = 1 - (0.08*M^1.45);
Cf_fuse = Cf_Turb(1) + Re_Ratio(1)*Cf_lam - Re_Ratio(1)*Cf_Turb(1);
CD0.fuse = Cf_fuse*fLD*fM*(Swet.fuse/Sref);

%-----Wing-----%
% Eqn 3.19
Cf_wing = Cf_Turb(2) + Re_Ratio(2)*Cf_lam - Re_Ratio(2)*Cf_Turb(2); % Cf laminar/turb mixed flow
CD0.wing = Cf_wing * ftc_w * fM * (Swet.wing/Sref)*((Cdmin.wing / 0.004)^0.4);

%-----Horizontal Tail-----%
% Eqn 3.20
Cf_ht = Cf_Turb(3) + Re_Ratio(3)*Cf_lam - Re_Ratio(3)*Cf_Turb(3); % Cf laminar/turb mixed flow
CD0.ht = Cf_ht * ftc_ht * fM * (Swet.ht/Sref)*((Cdmin.ht/ 0.004)^0.4);

%-----Vertical Tail-----%
% Eqn 3.21
Cf_vt = Cf_Turb(4) + Re_Ratio(4)*Cf_lam - Re_Ratio(4)*Cf_Turb(4); % Cf laminar/turb mixed flow
CD0.vt = Cf_vt * ftc_vt * fM * (Swet.vt/Sref)*((Cdmin.vt/ 0.004)^0.4);

%-----Nacelles-----%
Cf_nacelle = Cf_Turb(5) + Re_Ratio(5)*Cf_lam - Re_Ratio(5)*Cf_Turb(5); % Cf laminar/turb mixed flow
CD0.nacelle = nNacelle * Cf_nacelle*fLD_nac*fM*(Swet.nacelle/Sref);

%-----Protruding / Miscellaneous Components-----%
% Accounted for with correction factor below.

%-----Trim Drag Contribution-----%
WING.incidence = 3; % wing incidence angle (deg)
Cl.wing = WING.incidence * (pi/180) * (4/sqrt((M^2)-1)); 
Kwing = 0.85; % correction factor for tail
CL.wing = Cl.wing * Kwing; 

% 3D Lift Coefficient entire aircraft
CL.overall = 2*Wt.WTO / (rho * (V^2) * Sref);

% 3D Tail lift coefficient
CL.tail = (CL.overall - CL.wing) * (Sref / TAIL.Sh);

% Trim Drag CD0 --> eqn 3.37
ARTail = (TAIL.bh^2) / TAIL.Sh; % ht aspect ratio
eTail = 0.8; % tail oswald efficiency factor
CD0.trim = (1 / (pi * eTail * ARTail)) * (TAIL.Sh/Sref) * (CL.tail^2);

%-----TOTAL CD0-----%
CD0.correction = 1.1; % correction factor for misc items (Table 3.5, jet transport)
% CD0 for supersonic cruise
CD0.total.super = CD0.correction*(CD0.wing + CD0.ht + CD0.vt + CD0.fuse + CD0.nacelle + CD0.trim);

%------Wave Drag------%
% CD_Wave (3.5.2 Aircraft Wave Drag)

% Lift Dependent Wave drag
CDw.K_wl = 2 * (WING.geom.S_area / (WING.geom.span * L_F));
CDw.lift = CDw.K_wl * WING.geom.S_area * (CL.overall^2) * (M^2 - 1)/(2 * pi * L_F^2);

% Volume Dependent Wave Drag
CDw.beta = sqrt(M^2 - 1);
CDw.K_wv = 1.17 * ((1 + 0.75*CDw.beta * WING.geom.span / L_F)/(1+2*CDw.beta * WING.geom.span / L_F));
CDw.vol_aircraft = 4956.42; %ft^3 -> obtained from solidworks
CDw.volume = 128 * CDw.K_wv * CDw.vol_aircraft^2 / (pi * WING.geom.S_area * L_F^4);

CDw.total = CDw.lift + CDw.volume;

% Wave Drag
CD.wave = CDw.volume + CDw.K_wl .* WING.geom.S_area .* (CL.overall.^2) .* (M^2 - 1)/(2 .* pi .* L_F^2);

%% Subsonic Cruise

h = 37; % kft
M = 0.8; % Mach Number
visc = 0.297E-6; % viscosity (slugs/ft-s)

% Call Function to calculates dynamic pressure, reynolds number, and skin
% friction coefficients
[V, rho, q ,Cf_Turb,Cf_lam,Re_Ratio] = dragCalc(h, M, visc, L_F, WING.geom.MAC, TAIL.ch, TAIL.cv, nacLength);

%-----Fuselage-----%
% Eqn 3.14
fM = 1 - (0.08*M^1.45);
Cf_fuse = Cf_Turb(1) + Re_Ratio(1)*Cf_lam - Re_Ratio(1)*Cf_Turb(1);
CD0.fuse = Cf_fuse*fLD*fM*(Swet.fuse/Sref);

%-----Wing-----%
% Transonic drag predictions for wing are based on Ch 5 of Roskam Airplane
% Aero and Perf.
RNw = 2.9368e7; % RE of wing
Rwf =  1.02; % wing-fuselage interference factor (Figure 5.11 Roskam)
RLS = 1.25; % lifting surface correction factor (Figure 5.12 Roskam) 
Cfw = 0.455 / (((log10(RNw))^2.58)*((1 + 0.144*(M^2))^0.58)); % turbulent flat plate friction coefficient of the wing
Lprime = 1.2; % thickness location parameter; 1.2 for at (t/c)max >= 0.30c
tcWing = 0.0525; % thickness to camber ratio wing (NACA 0012 Airfoil)
CD0.wing = Rwf * RLS * Cfw * (1 + (Lprime * tcWing) + 100* ((tcWing)^4))*(Swet.wing/Sref);  
% Eqn 3.19
% Cf_wing = Cf_Turb(2) + Re_Ratio(2)*Cf_lam - Re_Ratio(2)*Cf_Turb(2); % Cf laminar/turb mixed flow
% CD0.wing = Cf_wing * ftc_w * fM * (Swet.wing/Sref)*((Cdmin.wing / 0.004)^0.4);

%-----Horizontal Tail-----%
% Eqn 3.20
Cf_ht = Cf_Turb(3) + Re_Ratio(3)*Cf_lam - Re_Ratio(3)*Cf_Turb(3); % Cf laminar/turb mixed flow
CD0.ht = Cf_ht * ftc_ht * fM * (Swet.ht/Sref)*((Cdmin.ht/ 0.004)^0.4);

%-----Vertical Tail-----%
% Eqn 3.21
Cf_vt = Cf_Turb(4) + Re_Ratio(4)*Cf_lam - Re_Ratio(4)*Cf_Turb(4); % Cf laminar/turb mixed flow
CD0.vt = Cf_vt * ftc_vt * fM * (Swet.vt/Sref)*((Cdmin.vt/ 0.004)^0.4);

%-----Nacelles-----%
Cf_nacelle = Cf_Turb(5) + Re_Ratio(5)*Cf_lam - Re_Ratio(5)*Cf_Turb(5); % Cf laminar/turb mixed flow
CD0.nacelle = nNacelle * Cf_nacelle*fLD_nac*fM*(Swet.nacelle/Sref);

%-----Protruding / Miscellaneous Components-----%
% Accounted for with correction factor below.

%-----Trim Drag Contribution-----%
%WING.incidence = 3; % wing incidence angle (deg)
%Cl.wing = WING.incidence * (pi/180) * (4/sqrt((M^2)-1)); 
%Kwing = 0.85; % correction factor for tail
%CL.wing = Cl.wing * Kwing; 

% 3D Lift Coefficient entire aircraft
%CL.overall = 2*Wt.WTO / (rho * (V^2) * Sref);

% 3D Tail lift coefficient
%CL.tail = (CL.overall - CL.wing) * (Sref / TAIL.Sh);

% Trim Drag CD0 --> eqn 3.37
%ARTail = (TAIL.bh^2) / TAIL.Sh; % ht aspect ratio
%eTail = 0.8; % tail oswald efficiency factor
%CD0.trim = (1 / (pi * eTail * ARTail)) * (TAIL.Sh/Sref) * (CL.tail^2);

%-----TOTAL CD0-----%
CD0.correction = 1.1; % correction factor for misc items (Table 3.5, jet transport)
% CD0 for subsonic cruise
CD0.total.sub = CD0.correction*(CD0.wing + CD0.ht + CD0.vt + CD0.fuse + CD0.nacelle);

%% Takeoff 
h = 0; % sea level
M = constraints.Vstall * 1.2/1116.5; 
visc = 0.374E-6; % viscosity (slugs/ft-s)
[V, rho, q ,Cf_Turb,Cf_lam,Re_Ratio] = dragCalc(h, M, visc, L_F, WING.geom.MAC, TAIL.ch, TAIL.cv, nacLength);

%-----Fuselage-----%
% Eqn 3.14
fM = 1 - (0.08*M^1.45);
Cf_fuse = Cf_Turb(1) + Re_Ratio(1)*Cf_lam - Re_Ratio(1)*Cf_Turb(1);
CD0.fuse = Cf_fuse*fLD*fM*(Swet.fuse/Sref);

%-----Wing-----%
% Eqn 3.19
Cf_wing = Cf_Turb(2) + Re_Ratio(2)*Cf_lam - Re_Ratio(2)*Cf_Turb(2); % Cf laminar/turb mixed flow
CD0.wing = Cf_wing * ftc_w * fM * (Swet.wing/Sref)*((Cdmin.wing / 0.004)^0.4);

%-----Horizontal Tail-----%
% Eqn 3.20
Cf_ht = Cf_Turb(3) + Re_Ratio(3)*Cf_lam - Re_Ratio(3)*Cf_Turb(3); % Cf laminar/turb mixed flow
CD0.ht = Cf_ht * ftc_ht * fM * (Swet.ht/Sref)*((Cdmin.ht/ 0.004)^0.4);

%-----Vertical Tail-----%
% Eqn 3.21
Cf_vt = Cf_Turb(4) + Re_Ratio(4)*Cf_lam - Re_Ratio(4)*Cf_Turb(4); % Cf laminar/turb mixed flow
CD0.vt = Cf_vt * ftc_vt * fM * (Swet.vt/Sref)*((Cdmin.vt/ 0.004)^0.4);

%-----Nacelles-----%
Cf_nacelle = Cf_Turb(5) + Re_Ratio(5)*Cf_lam - Re_Ratio(5)*Cf_Turb(5); % Cf laminar/turb mixed flow
CD0.nacelle = nNacelle * Cf_nacelle*fLD_nac*fM*(Swet.nacelle/Sref);

%-----Protruding / Miscellaneous Components-----%
% Accounted for with correction factor below.

%-----High Lift Devices-----%
% Trailing Edge HLD (Eqn 3.27)
% Eqn 3.27 is based on a flap with a flap-span-to-wing-span ratio of 70%
% Constants (Table 3.3)
A_TEflap = 0.0011; % double-slotted flap
B_TEflap = 1; 
deltaFlap = 15; % flap deflection angle (deg)
Cf_C_Ratio = 0.3; % ratio between avg flap chord to avg wing chord
CD0.TEflap = Cf_C_Ratio*A_TEflap*(deltaFlap^B_TEflap);

% Leading Edge HLD (Eqn 3.28)
Csl_C_Ratio = 0.15; % ratio between avg extended slat chord and extended wing chord
CD0.LEslat = (Csl_C_Ratio)*CD0.wing;

%-----Landing Gear-----%
nWhlMain = 4; % number of main gear assembly wheels
nWhlNose = 2; % number of nose gear assembly wheels
CDlg = 0.30; % LG with no fairing

% Nose Gear
D.nose = 1.84; % wheel diameter(ft)
w.nose = 0.57; % wheel width (ft)
Slg.nose = D.nose*w.nose; % frontal area of nose gear wheel
CD0.LGNose = nWhlNose*CDlg*(Slg.nose/Sref); 

% Main Gear
D.main = 2.63; % wheel diameter (ft)
w.main = 0.81; % wheel width (ft)
Slg.main = D.main*w.main; % frontal area (ft)
CD0.LGMain = nWhlMain*CDlg*(Slg.main/Sref);

CD0.LGTotal = CD0.LGNose + CD0.LGMain; 

%-----TOTAL CD0-----%
CD0.correction = 1.1; % correction factor for misc items (Table 3.5, jet transport)
% CD0 for takeoff configuration
CD0.total.TO = CD0.correction*(CD0.wing + CD0.ht + CD0.vt + CD0.fuse + CD0.nacelle + CD0.TEflap + CD0.LEslat + CD0.LGTotal);

%% Landing

h = 0; % sea level
M = constraints.Vstall * 1.15/1116.5;
visc = 0.374E-6; % viscosity (slugs/ft-s)
[V, rho, q ,Cf_Turb,Cf_lam,Re_Ratio] = dragCalc(h, M, visc, L_F, WING.geom.MAC, TAIL.ch, TAIL.cv, nacLength);

%-----Fuselage-----%
% Eqn 3.14
fM = 1 - (0.08*M^1.45);
Cf_fuse = Cf_Turb(1) + Re_Ratio(1)*Cf_lam - Re_Ratio(1)*Cf_Turb(1);
CD0.fuse = Cf_fuse*fLD*fM*(Swet.fuse/Sref);

%-----Wing-----%
% Eqn 3.19
Cf_wing = Cf_Turb(2) + Re_Ratio(2)*Cf_lam - Re_Ratio(2)*Cf_Turb(2); % Cf laminar/turb mixed flow
CD0.wing = Cf_wing * ftc_w * fM * (Swet.wing/Sref)*((Cdmin.wing / 0.004)^0.4);

%-----Horizontal Tail-----%
% Eqn 3.20
Cf_ht = Cf_Turb(3) + Re_Ratio(3)*Cf_lam - Re_Ratio(3)*Cf_Turb(3); % Cf laminar/turb mixed flow
CD0.ht = Cf_ht * ftc_ht * fM * (Swet.ht/Sref)*((Cdmin.ht/ 0.004)^0.4);

%-----Vertical Tail-----%
% Eqn 3.21
Cf_vt = Cf_Turb(4) + Re_Ratio(4)*Cf_lam - Re_Ratio(4)*Cf_Turb(4); % Cf laminar/turb mixed flow
CD0.vt = Cf_vt * ftc_vt * fM * (Swet.vt/Sref)*((Cdmin.vt/ 0.004)^0.4);

%-----Nacelles-----%
Cf_nacelle = Cf_Turb(5) + Re_Ratio(5)*Cf_lam - Re_Ratio(5)*Cf_Turb(5); % Cf laminar/turb mixed flow
CD0.nacelle = nNacelle * Cf_nacelle*fLD_nac*fM*(Swet.nacelle/Sref);

%-----Protruding / Miscellaneous Components-----%
% Accounted for with correction factor below.

%-----High Lift Devices-----%
% Trailing Edge HLD (Eqn 3.27)
% Eqn 3.27 is based on a flap with a flap-span-to-wing-span ratio of 70%
% Constants (Table 3.3)
A_TEflap = 0.0011; % double-slotted flap
B_TEflap = 1; 
deltaFlap = 40; % flap deflection angle (deg)
Cf_C_Ratio = 0.3; % ratio between avg flap chord to avg wing chord
CD0.TEflap = Cf_C_Ratio*A_TEflap*(deltaFlap^B_TEflap);

% Leading Edge HLD (Eqn 3.28)
Csl_C_Ratio = 0.15; % ratio between avg extended slat chord and extended wing chord
CD0.LEslat = (Csl_C_Ratio)*CD0.wing;

%-----Landing Gear-----%
nWhlMain = 4; % number of main gear assembly wheels
nWhlNose = 2; % number of nose gear assembly wheels
CDlg = 0.30; % LG with no fairing

% Nose Gear
D.nose = 1.84; % wheel diameter(ft)
w.nose = 0.57; % wheel width (ft)
Slg.nose = D.nose*w.nose; % frontal area of nose gear wheel
CD0.LGNose = nWhlNose*CDlg*(Slg.nose/Sref); 

% Main Gear
D.main = 2.63; % wheel diameter (ft)
w.main = 0.81; % wheel width (ft)
Slg.main = D.main*w.main; % frontal area (ft)
CD0.LGMain = nWhlMain*CDlg*(Slg.main/Sref);

CD0.LGTotal = CD0.LGNose + CD0.LGMain; 

%-----TOTAL CD0-----%
CD0.correction = 1.1; % correction factor for misc items (Table 3.5, jet transport)
% CD0 for takeoff configuration
CD0.total.Land = CD0.correction*(CD0.wing + CD0.ht + CD0.vt + CD0.fuse + CD0.nacelle + CD0.TEflap + CD0.LEslat + CD0.LGTotal);

%% Plot Drag Polars

% Induced drag
e_Oswald = 0.8; % efficiency factor
CL.array = linspace(0, 2, 20);
CD.induced = (CL.array.^2) ./ (pi * e_Oswald * WING.geom.AR);

% Supersonic Cruise Configuration
CD.total.super = CD0.total.super + CD.induced + CD.wave; 

% Subsonic Cruise Configuration
CD.total.sub = CD0.total.sub + CD.induced;

% Takeoff Configuration 
CD.total.TO = CD0.total.TO + CD.induced;

% Landing Configuration
CD.total.Land = CD0.total.Land + CD.induced;

figure()
plot(CL.array, CD.total.super,'b--') % Supersonic Cruise Configuration
hold on;
plot(CL.array, CD.total.sub,'g') % Subsonic Cruise Configuration
plot(CL.array, CD.total.TO) % Takeoff Configuration 
plot(CL.array, CD.total.Land) % Landing Configuration
xlabel('Lift Coefficient, C_L')
ylabel('Drag Coefficient, C_D')
title('Drag Polars')
legend('Supersonic Cruise (M=1.6)', 'Subsonic Cruise (M=0.8)', ...
    'Takeoff Configuration', 'Landing Configuration', 'Location', 'Northwest')