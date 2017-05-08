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
CD.cr.begin = spline(CL.array, CD.total.super, CL.cr.begin);
CD.cr.end = spline(CL.array, CD.total.super, CL.cr.end);
L_D.cr.begin = CL.cr.begin / CD.cr.begin;
L_D.cr.end = CL.cr.end / CD.cr.end;

% Determine L/D from WING.CFD arrays
CL.cr.CFD = spline(WING.CFD.SUP.alpha,WING.CFD.SUP.CL, 2.86); 
CD.cr.CFD = spline(WING.CFD.SUP.alpha,WING.CFD.SUP.CD, 2.86);
L_D.cr.CFD = CL.cr.CFD / CD.cr.CFD;

% Average start and end CL/CD values
%L_D.cr.average = (L_D.cr.begin + L_D.cr.end) / 2;
L_D.cr.average = 7.5; 

% Determine range (nm) of aircraft for constant Mach number (eqn 11.62)
cj = Wt.fuel.sfc_cr * 3600; % lbs/h/lbs GET UPDATED VALUE FROM PARA ANALYSIS
Perf.Rng.ConsM = ((a.soundCr * sqrt(atm.theta) * req.cr_M0(1) * L_D.cr.CFD) / cj) * log(Wt.cr.begin/Wt.cr.end);
fprintf('The maximum range (Constant Mach Number) is %0.2f nm. \n', Perf.Rng.ConsM)
%% Determine Endurance

%Wt.ltr.start = Wt.cr.end; % start loiter at end of cruise phase
%Wt.ltr.end = Wt.fuel.w5_4 * Wt.ltr.start; % end of loiter weight
Wt.end = Wt.start - Wt.fuel.w_tot; 

% Transonic Endurance
cj_trans = 0.484/3600; 
Kend = 1 / (pi*0.8*3); 
CL.trans.ed = sqrt(CD0.total.sub / Kend);
CD.trans.ed = 2*CD0.total.sub;
L_D.trans.ed = CL.trans.ed / CD.trans.ed;
Perf.End.ConsAOA.sub = (1 / cj_trans) * (L_D.trans.ed) * (log(Wt.start / Wt.end));
fprintf('The maximum endurance (cruise climb) is %0.2f hours. \n', Perf.End.ConsAOA.sub/3600)

% Supersonic Endurance
CL.super.ed = sqrt(CD0.total.super / Kend);
CD.super.ed = 2*CD0.total.super;
L_D.super.ed = CL.super.ed / CD.super.ed;
Perf.End.ConsAOA.super = (1 / Wt.fuel.sfc_cr) * (L_D.super.ed) * (log(Wt.start / Wt.end));
fprintf('The maximum endurance (cruise climb) is %0.2f hours. \n', Perf.End.ConsAOA.super/3600)

% % Determine Endurance (hours) for constant AOA
% Perf.End.ConsM = (1/cj_trans)*(L_D.cr.average)*log(Wt.ltr.start/Wt.ltr.end);

%% Determine Range (Const Altitude)

CL.alt.max = sqrt((pi*WING.geom.AR * e_Oswald * CD0.total.super)/3); % eqn 11.66
CD.alt.max = (4/3)*CD0.total.super; % eqn 11.66
sqrtCL_CD = sqrt(CL.alt.max) / CD.alt.max; 

Perf.Rng.ConsAlt = (1.675 / (cj * sqrt(atm.crRange * WING.geom.S_area)))*(sqrtCL_CD) * (sqrt(Wt.cr.begin) - sqrt(Wt.cr.end));
fprintf('The maximum range (const altitude) is %0.2f nm. \n', Perf.Rng.ConsAlt)

%% Payload Range Diagram

% Determine max range from const alt and const mach number methods
Perf.Rng.max = max(Perf.Rng.ConsAlt,Perf.Rng.ConsM);
% Useful Weight
Wt.Useful = Wt.WTO - Wt.WOEW;
% Linearly decrease from max payload to 0
Wt.pld.linear = linspace(Wt.pld.w_tot, 0, 10); 
% Linearly increase fuel to max useful weight
Wt.fuel.linear = linspace(Wt.fuel.w_tot, Wt.Useful, 10); 
% See Roskam and Lam p.560




% Pt C --> Rmax,fuel
% Mission fuel fraction for various fuel weights
Wt.mff.array = 1 - (Wt.fuel.linear / Wt.WTO);

% Cruise Fuel Fraction Array
Wt.fuel.cr_array = Wt.mff.array*(1/Wt.fuel.w2_1)*(1/Wt.fuel.w3_2)*(1/Wt.fuel.w5_4)*(1/Wt.fuel.w6_5)*(1/Wt.fuel.w7_6); 

% Weight at end of cruise for varying cruise fuel fraction
Wt.cr.end_array = Wt.fuel.cr_array * Wt.cr.begin; % weight at end of cruise phase

Perf.Rng.ConsM_array = ((a.soundCr .* sqrt(atm.theta) .* req.cr_M0(1) .* L_D.cr.average) ./ cj) .* log(Wt.cr.begin./Wt.cr.end_array);


% Pt B --> R = Rmax, Max payload
Perf.Rng.Array = [0,Perf.Rng.max];
Wt.pld.Array = [Wt.pld.w_tot,Wt.pld.w_tot];
figure()
title('Payload-Range Diagram')
xlabel('Range, nm')
ylabel('Payload Weight (lbs)')
hold on;

% Pt A --> R = 0 nm, Max payload
%plot(Perf.Rng.Array(1), Wt.pld.w_tot, 'o', 'Color', [1 0.5 0], 'MarkerSize', 8)
plot(Perf.Rng.Array(2), Wt.pld.w_tot, 's', 'Color', [1 0.5 0], 'MarkerSize', 8)


% Pt D --> Rmax
plot(Perf.Rng.ConsM_array(end), Wt.pld.linear(end),'x', 'Color', [1 0.5 0], 'MarkerSize', 8)
% Linearly increase fuel, decrease payload linearly
plot(Perf.Rng.ConsM_array, Wt.pld.linear,'b')
legend('R_{design}', 'R_{ferry}','Location', 'Southwest')
plot(Perf.Rng.Array, Wt.pld.Array,'b')
