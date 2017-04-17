%% Payload Range 

% Creates payload range diagram based on the document posted by Professor
% Van Dam

% Run iter_weights.m --> TO_TD_perf.m --> Drag_Analysis.m before this script

%% Range - Constant Altitude
% Range at constant altitude (constant air density, ?, and constant air pressure, p)
% and Mach number, M, is (assuming thrust specific fuel consumption, c, remains constant throughout
% the mission)

% Dynamic pressure (psf)
atm.crRange = atm.rho_sl * atm.sig_rho; % density (slugs/ft^3)
V_cr = req.cr_M0(1) * Wt.fuel.a_snd * 1116.5; % cruise speed (ft/s)
q_cr = 0.5*atm.crRange*(V_cr^2); % dynamic pressure (psf)

% Drag Forces
D_0 = CD0.total.cruise * q_cr * WING.geom.S_area;

% Constants 
K_perf = 1 / (pi * e * WING.geom.AR);
b1 = sqrt(D_0);
b2 = sqrt((K_perf / WING.geom.S_area) * (1/q_cr));

% Pt. B
Wt.fuel.reserve = Wt.fuel.w_max - Wt.fuel.w_tot; % Amount of reserve fuel
Wt.W1.B = Wt.WOEW + Wt.pld.w_tot + Wt.fuel.reserve; % W1 = OEW + Max Pld + Fuel reserve
Range.c_alt.B = (V_cr / Wt.fuel.sfc_cr) * (1 / (b1 * b2)) * (atan((0.97 * Wt.WTO) * (b2/b1)) - atan(Wt.W1.B *(b2/b1)));

% Pt. C
Wt.W1.C = Wt.WTO - Wt.fuel.w_max + Wt.fuel.reserve;
Range.c_alt.C = (V_cr / Wt.fuel.sfc_cr) * (1 / (b1 * b2)) * (atan((0.97 * Wt.WTO) * (b2/b1)) - atan(Wt.W1.C *(b2/b1)));

% Pt. D
Wt.W0.D = Wt.WOEW + Wt.fuel.w_max;
Wt.W1.D = Wt.WOEW;
Range.c_alt.D = (V_cr / Wt.fuel.sfc_cr) * (1 / (b1 * b2)) * (atan((Wt.W0.D) * (b2/b1)) - atan(Wt.W1.D *(b2/b1)));

% Store ranges in array
Range.c_alt.array = [Range.c_alt.B, Range.c_alt.C, Range.c_alt.D];

% Convert ft to nmi
Range.c_alt.array = Range.c_alt.array * 0.000164579; 

% Plot Payload Range Diagram
figure()
plot(0, Wt.pld.w_tot, 'o') % pt. A
hold on;
plot(Range.c_alt.array(1), Wt.pld.w_tot, 'o') % pt. B
plot(Range.c_alt.array(2), Wt.WTO - Wt.fuel.w_max, 'o') % pt. C
plot(Range.c_alt.array(3), 0, 'o') % pt. D
title('Supersonic Cruise Payload Range Diagram')
xlabel('Range (nmi)')
ylabel('Payload (lb)')
%% Range at Cruise Climb Conditions
% Range at cruise climb conditions (constant lift coefficient, CL, and Mach number, M) and
% assuming thrust specific fuel consumption, c, remains constant throughout the mission

L_D.cr.average = 7.5;

% Pt. B
Range.cr_climb.B = ((Wt.fuel.a_snd * 1116.5) / Wt.fuel.sfc_cr) * req.cr_M0(1) * L_D.cr.average * log((0.97 * Wt.WTO) / Wt.W1.B);

% Pt. C
Range.cr_climb.C = ((Wt.fuel.a_snd * 1116.5) / Wt.fuel.sfc_cr) * req.cr_M0(1) * L_D.cr.average * log((0.97 * Wt.WTO) / Wt.W1.C);

% Pt. D
Range.cr_climb.D = ((Wt.fuel.a_snd * 1116.5) / Wt.fuel.sfc_cr) * req.cr_M0(1) * L_D.cr.average * log((Wt.W0.D) / Wt.W1.D);

% Store ranges in array
Range.cr_climb.array = [Range.cr_climb.B, Range.cr_climb.C, Range.cr_climb.D];


