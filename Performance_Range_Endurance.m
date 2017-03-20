%% Performance: Range and Endurance

% Date Created 3/19/17
% Calculates the range and endurance of the aircraft assuming constant Mach
% number (cruise altitude may vary slightly)

% NOTE: iter_weights.m must be run before this script
% This uses weights and fuel fractions calculated from that script. 
% Reference: Airplane Aerodynamics and Performance Roskam and Lam Section
% 11.2

close all; 
Drag_Analysis;
%% Determine Range (Const Mach Number)

% Fuel fractions from iter_weights.m
Wt.start = Wt.WTO; % starting weight (doesn't necessarily have to be WTO)
Wt.climb = Wt.fuel.w2_1 * Wt.start; % weight at start of climb
Wt.cr.begin = Wt.fuel.w3_2 * Wt.climb; % weight at start of cruise phase
Wt.cr.end = Wt.fuel.w4_3 * Wt.cr.begin; % weight at end of cruise phase

a.soundCr = Wt.fuel.a_snd * 1116.5; % speed of sound at cruise alt (ft/s)
a.soundCr = a.soundCr * 0.5924; % convert to knots 

% Dynamic pressure (psf)
atm.crRange = atm.rho_sl * atm.sig_rho; % density (slugs/ft^3)
V_cr = req.cr_M0(1) * Wt.fuel.a_snd * 1116.5; % cruise speed (ft/s)
q_cr = 0.5*atm.crRange*(V_cr^2); % dynamic pressure (psf)

% Since L/D will vary during the mission an average will used for the
% beginning and end of the cruise phase (see p.546 Roskam and Lam)
CL.cr.begin = Wt.cr.begin / (q_cr * WING.geom.S_area);
CL.cr.end = Wt.cr.end / (q_cr * WING.geom.S_area); 

% Interpolate CD values from drag polar
CD.cr.begin = spline(CL.overall, CD.total.super, CL.cr.begin);
CD.cr.end = spline(CL.overall, CD.total.super, CL.cr.end);
L_D.cr.begin = CL.cr.begin / CD.cr.begin;
L_D.cr.end = CL.cr.end / CD.cr.end;

% Average start and end CL/CD values
L_D.cr.average = (L_D.cr.begin + L_D.cr.end) / 2;
%L_D.cr.average = 7;
% Determine range (nm) of aircraft for constant Mach number (eqn 11.62)
cj = Wt.fuel.sfc_cr * 3600; % lbs/h/lbs GET UPDATED VALUE FROM PARA ANALYSIS
Perf.Rng.ConsM = ((a.soundCr * sqrt(atm.theta) * req.cr_M0(1) * L_D.cr.average) / cj) * log(Wt.cr.begin/Wt.cr.end);
fprintf('The maximum range is %0.2f nm. \n', Perf.Rng.ConsM)
%% Determine Endurance (Const Mach Number)

Wt.ltr.start = Wt.cr.end; % start loiter at end of cruise phase
Wt.ltr.end = Wt.fuel.w5_4 * Wt.ltr.start; % end of loiter weight

% Loiter set to 1h?
% different cj
% Determine Endurance (hours) for constant Mach number (eqn 11.63)
Perf.End.ConsM = (1/cj)*(L_D.cr.average)*log(Wt.ltr.start/Wt.ltr.end);
fprintf('The maximum endurance (const Mach number) is %0.2f hours \n.', Perf.End.ConsM)


%% Determine Range (Const Altitude)

CL.alt.max = sqrt((pi*WING.geom.AR * e_Oswald * CD0.total.cruise)/3); % eqn 11.66
CD.alt.max = (4/3)*CD0.total.cruise; % eqn 11.66
sqrtCL_CD = sqrt(CL.alt.max) / CD.alt.max; 

Perf.Rng.ConsAlt = (1.675 / (cj * sqrt(atm.crRange * WING.geom.S_area)))*(sqrtCL_CD) * (sqrt(Wt.cr.begin) - sqrt(Wt.cr.end));
fprintf('The maximum range (const altitude) is %0.2f nm \n.', Perf.Rng.ConsAlt)