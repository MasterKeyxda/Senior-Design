%% Drag Analysis

% This script provides an estimate for CD0
% Drag breakdown procedure is based on the Aircraft Performance Analysis
% by M. Sadraey Ch 3
% Date Created: 3/17/17

% Run iter_weights.m before running this script
%% Preliminary Info

%-----Standard Day, Cruise Altitude Conditions-----%

% Cruise Altitude = 42 kft
atm.densCr =  atm.rho_sl.*atm.sig_rho; % density (slugs/ft^3)
atm.viscCr = 0.297E-6; % viscosity at 42 kft (slugs/ft-s)
atm.aSoundCr = Wt.fuel.a_snd * 1116.5; % speed of sound (ft/s)

% Sea Level Properties
atm.viscSL = 3.74*10^-7; % viscosity 

%-----Reynolds Numbers, Dynamic Pressure, Skin Friction-----%

% Cruise Velocity at cruise ceiling
Vcr = req.cr_M0(1) * atm.aSoundCr; % ft/s
% Nacelle Diameters (from Solidworks)
nacLength = 16.79; % nacelle length ft (Solidworks/keyur)
% nacDMax; % nacelle ft

% Reynolds numbers (Cruise)
Re.crit = 1E6; % critical reynolds number
Re.fuse = (atm.densCr * Vcr * L_F) / atm.viscCr;   % based on fuselage length
Re.wing = (atm.densCr * Vcr * WING.geom.MAC) / atm.viscCr;  % based on mean geo chord of wetted part of wing
Re.ht = (atm.densCr * Vcr * TAIL.ch) / atm.viscCr; % based on mean geo chord of wetted part of horiz tail
Re.vt = (atm.densCr * Vcr * TAIL.cv) / atm.viscCr; % based on mean geo chord of wetted part of vert tail 
Re.nacelle = (atm.densCr * Vcr * nacLength) / atm.viscCr; % based on nacelle length

% Reynolds numbers (Landing) --> based on actual approach speed; see
% TO_TD_perf.m
Re.fuseLand = (atm.rho_sl * VTDActual * L_F) / atm.viscSL;
Re.wingLand = (atm.rho_sl * VTDActual * WING.geom.MAC) / atm.viscSL;
Re.htLand = (atm.rho_sl * VTDActual * TAIL.ch) / atm.viscSL;
Re.vtSL = (atm.rho_sl * VTDActual * TAIL.cv) / atm.viscSL;
Re.nacelleSL = (atm.rho_sl * VTDActual * nacLength) / atm.viscSL;

% Reynolds numbers (Takeoff) --> based on lift off speed; see TO_TD_perf.m
Re.fuseTO = (atm.rho_sl * V_LOF * L_F) / atm.viscSL;
Re.wingTO = (atm.rho_sl * V_LOF * WING.geom.MAC) / atm.viscSL;
Re.htTO = (atm.rho_sl * V_LOF * TAIL.ch) / atm.viscSL;
Re.vtTO = (atm.rho_sl * V_LOF * TAIL.cv) / atm.viscSL;
Re.nacelleTO = (atm.rho_sl * V_LOF * nacLength) / atm.viscSL;

% Reynolds Number Ratios
Re.fuseRatio = Re.crit / Re.fuse; 
Re.wingRatio = Re.crit / Re.wing;
Re.htRatio = Re.crit / Re.ht;
Re.vtRatio = Re.crit / Re.vt; 
Re.nacelleRatio = Re.crit / Re.nacelle; 

% Skin Friction Coefficients 
Cf.lam = 1.327/sqrt(1E6); % laminar flow (transition at 1E6)
% Assume fully Turbulent Flow --> eqn 3.15a
Cf.fuseTurb = 0.455 / ((log10(Re.fuse))^2.58); % Fuselage
Cf.wingTurb = 0.455 / ((log10(Re.wing))^2.58); % Wing
Cf.htTurb = 0.455 / ((log10(Re.ht))^2.58); % Horizontal Tail
Cf.vtTurb = 0.455 / ((log10(Re.vt))^2.58); % Vertical Tail
Cf.nacelleTurb = 0.455 / ((log10(Re.nacelle))^2.58);
q = 0.5*atm.densCr*(Vcr^2); % dynamic pressure, standard day, cruise altitude

%% CDO Calculations

% Common Equation Parameters
Sref = WING.geom.S_area; % Wing planform area chosen as reference area
fM = 1 - (0.08*req.cr_M0(1)^1.45);

%-----Wing-----%

% Eqn 3.19
Swet.wing = WING.geom.S_area - 2*0.5*D_C*(WING.geom.sub.Cr +...
    (WING.geom.sub.Cr+(0.5*D_C)*((WING.geom.sub.Ct-WING.geom.sub.Cr)/WING.geom.sub.span)))/2; % wing wetted area correction for removed area inside fuselage
Swet.wing = 2*Swet.wing; % wing wetted area
wing_TCRatio = 0.0525; % wing max thickness to chord ratio
ftc_w = 1 + (2.7*wing_TCRatio) + 100*((wing_TCRatio)^4); % ftc_w is a function of thickness ratio
Cdmin.wing = min(WING.biconvex(2).Cd);
Cf.wing = Cf.wingTurb + Re.wingRatio*Cf.lam - Re.wingRatio*Cf.wingTurb; % Cf laminar/turb mixed flow
CD0.wing = Cf.wing * ftc_w * fM * (Swet.wing/Sref)*((Cdmin.wing / 0.004)^0.4);

%-----Horizontal Tail-----%

% Eqn 3.20
Swet.ht = 2*TAIL.Sh; % wetted area of horizontal tail
HT_TCRatio = 0.0525; % horizontal tail max thickness to chord ratio
ftc_ht = 1 + (2.7*HT_TCRatio) + 100*((HT_TCRatio)^4);
Cdmin.ht = min(WING.biconvex(2).Cd);
Cf.ht = Cf.htTurb + Re.htRatio*Cf.lam - Re.htRatio*Cf.htTurb;
CD0.ht = Cf.ht * ftc_ht * fM * (Swet.ht/Sref)*((Cdmin.ht/ 0.004)^0.4);

%-----Vertical Tail-----%

% Eqn 3.21
Swet.vt = 2*TAIL.Sv; % wetted area of vertical tail
VT_TCRatio = 0.0525; % vertical tail max thickness to chord ratio
ftc_vt = 1 + (2.7*VT_TCRatio) + 100*((VT_TCRatio)^4);
Cdmin.vt = min(WING.biconvex(2).Cd); 
Cf.vt = Cf.vtTurb + Re.vtRatio*Cf.lam - Re.vtRatio*Cf.vtTurb;
CD0.vt = Cf.vt * ftc_vt * fM * (Swet.vt/Sref)*((Cdmin.vt/ 0.004)^0.4);

%-----Fuselage-----%

% Eqn 3.14
% LF_DF refers to fuselage fineness ratio (from Preliminary_Fuselage_Design.m)
fLD = 1 + (60 / (LF_DF^3)) + 0.0025*(LF_DF);
Swet.fuse = 2773; % from Solidworks (ft^2) prior 2696.24
Cf.fuse = Cf.fuseTurb + Re.fuseRatio*Cf.lam - Re.fuseRatio*Cf.fuseTurb;
CD0.fuse = Cf.fuse*fLD*fM*(Swet.fuse/Sref);

%-----Nacelles-----%

% "For the purpose of drag calculation, it can be considered that the nacelle is
% similar to the fuselage, except its length-to-diameter ratio is lower. Thus, the nacelle zerolift
% drag coefficient (CDon) will be determined in the same way as in a fuselage. In the
% case where the nacelle length-to-diameter ratio is below 2, assume 2."
nacL_D = 2; % ASSUMPTION nacelle fineness ratio
% nacRadius = nacDMax / 2; % nacelle radius (ft)
fLD_nac = 1 + (60 / (nacL_D^3)) + 0.0025*(nacL_D);
nNacelle = 2;
Swet.nacelle = 47.09; % from Solidworks (ft^2)
CD0.nacelle = nNacelle * Cf.nacelleTurb*fLD_nac*fM*(Swet.nacelle/Sref);
% CD0.nacelle = 0.0012; % APPROXIMATION; based on nacelle CD0 for Gates Learjet 25 business jet

%-----High Lift Devices-----%

% Trailing Edge HLD (Eqn 3.27)
% Eqn 3.27 is based on a flap with a flap-span-to-wing-span ratio of 70%
% Constants (Table 3.3)
A_TEflap = 0.0011; % double-slotted flap
B_TEflap = 1; 
deltaFlap = 40; % flap deflection angle (deg)
Cf_C_Ratio = 0.3; % ratio between avg flap chord to avg wing chord
CD0.TEflap = Cf_C_Ratio*A_TEflap*(deltaFlap^B_TEflap);
% 
% Leading Edge HLD (Eqn 3.28)
Csl_C_Ratio = 0.3; % ratio between avg extended slat chord and extended wing chord
CD0.LEslat = (Csl_C_Ratio)*CD0.wing;

%-----Landing Gear Retracted-----%

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

%-----Protruding / Miscellaneous Components-----%
% Accounted for with correction factor below.

%-----Trim Drag Contribution-----%

WING.incidence = 3; % wing incidence angle (deg)
Cl.wing = WING.incidence * (pi/180) * (4/sqrt((req.cr_M0(1)^2)-1)); 
Kwing = 0.85; % correction factor for tail
CL.wing = Cl.wing * Kwing; 

% 3D Lift Coefficient entire aircraft
CL.overall = 2*Wt.WTO / (atm.densCr * (Vcr^2) * Sref);

% 3D Tail lift coefficient
CL.tail = (CL.overall - CL.wing) * (Sref / TAIL.Sh);

% Trim Drag CD0 --> eqn 3.37
ARTail = (TAIL.bh^2) / TAIL.Sh; % ht aspect ratio
eTail = 0.8; % tail oswald efficiency factor
CD0.trim = (1 / (pi * eTail * ARTail)) * (TAIL.Sh/Sref) * (CL.tail^2);

%-----TOTAL CD0-----%
CD0.correction = 1.1; % correction factor for misc items (Table 3.5, jet transport)
% CD0 for cruise
CD0.total.cruise = CD0.correction*(CD0.wing + CD0.ht + CD0.vt + CD0.fuse + CD0.nacelle + CD0.trim);

% CD0 for takeoff
% Account for HLDs and LG deployed
CD0.total.takeoff = CD0.correction*(CD0.wing + CD0.ht + CD0.vt + CD0.fuse + CD0.nacelle + CD0.trim) + CD0.LEslat + CD0.TEflap + CD0.LGTotal;

% Separate CD0 for takeoff and landing once HLDs calcs are finished; adjust
% control surface deflection angles for each scenario.
% CD0.total.land = CD0.correction*(CD0.wing + CD0.ht + CD0.vt + CD0.fuse + CD0.nacelle + CD0.trim) + CD0.LEslat + CD0.TEflap + CD0.LGTotal;
%% CD_Wave (3.5.2 Aircraft Wave Drag)

% Lift Dependent Wave drag
CDw.K_wl = 2 * (WING.geom.S_area / (WING.geom.span * L_F));
CDw.lift = CDw.K_wl * WING.geom.S_area * (CL.overall^2) * (req.cr_M0(1)^2 - 1)/(2 * pi * L_F^2);

% Volume Dependent Wave Drag
CDw.beta = sqrt(req.cr_M0(1)^2 - 1);
CDw.K_wv = 1.17 * ((1 + 0.75*CDw.beta * WING.geom.span / L_F)/(1+2*CDw.beta * WING.geom.span / L_F));
CDw.vol_aircraft = 4956.42; %ft^3 -> obtained from solidworks
CDw.volume = 128 * CDw.K_wv * CDw.vol_aircraft^2 / (pi * WING.geom.S_area * L_F^4);

CDw.total = CDw.lift + CDw.volume;

%% Plot Drag Polars

% Induced drag
e_Oswald = 0.8; % efficiency factor
CL.overall = linspace(0, 2, 20);
CD.induced = (CL.overall.^2) ./ (pi * e_Oswald * WING.geom.AR); 

% Wave Drag
CD.wave = CDw.volume + CDw.K_wl .* WING.geom.S_area .* (CL.overall.^2) .* (req.cr_M0(1)^2 - 1)/(2 .* pi .* L_F^2);

% Supersonic Cruise Drag Polar
CD.total.super = CD0.total.cruise + CD.induced + CD.wave; 
figure()
plot(CL.overall, CD.total.super)
xlabel('Lift Coefficient, C_L')
ylabel('Drag Coefficient, C_D')
title('Supersonic Cruise Drag Polar')

% Subsonic Drag Polar
CD.total.sub = CD0.total.cruise + CD.induced;
figure()
plot(CL.overall, CD.total.sub)
xlabel('Lift Coefficient, C_L')
ylabel('Drag Coefficient, C_D')
title('Subsonic Cruise Drag Polar')

% Takeoff and Landing Drag Polar
CD.total.takeoff = CD0.total.takeoff + CD.induced;
figure()
plot(CL.overall, CD.total.takeoff)
xlabel('Lift Coefficient, C_L')
ylabel('Drag Coefficient, C_D')
title('Takeoff/Landing Cruise Drag Polar')

% Plot all drag polars on same figure
figure()
hold on;
plot(CL.overall, CD.total.super, 'b--') % supersonic
plot(CL.overall, CD.total.sub, 'r:','LineWidth',1.5) % subsonic
plot(CL.overall, CD.total.takeoff, 'c') % takeoff/landing
xlabel('Lift Coefficient, C_L')
ylabel('Drag Coefficient, C_D')
legend('Supersonic Cruise', 'Subsonic Cruise', 'Takeoff/Landing', 'Location','Northwest')
title('Aircraft Drag Polars')