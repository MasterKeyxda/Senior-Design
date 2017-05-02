%% Payload Range 

% Creates payload range diagram based on the document posted by Professor
% Van Dam

% Run iter_weights.m --> TO_TD_perf.m --> Drag_Analysis.m before this script

%% Range - Constant Altitude (Supersonic)
% Range at constant altitude (constant air density, rho, and constant air pressure, p)
% and Mach number, M, is (assuming thrust specific fuel consumption, c, remains constant throughout
% the mission)

% Dynamic pressure (psf)
atm.crRange = atm.rho_sl * atm.sig_rho; % density (slugs/ft^3)
V_cr = req.cr_M0(1) * Wt.fuel.a_snd * 1116.5; % cruise speed (ft/s)
q_cr = 0.5*atm.crRange*(V_cr^2); % dynamic pressure (psf)

% Drag Forces
D_0 = CD0.total.super * q_cr * WING.geom.S_area;

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

%% Range at Constant Altitude (Transonic M = 0.8)
h = 37; % transonic cruise altitude (kft)
[~,~,atm.sig_trans,atm.snd_trans] = AltTable(h,'h');

% Dynamic pressure (psf)
atm.transCr = atm.rho_sl * atm.sig_trans; % density (slugs/ft^3)
V_transCr = 0.8 * atm.snd_trans * 1116.5; % transonic cruise speed (ft/s)
q_trans = 0.5*atm.transCr*(V_transCr^2); % dynamic pressure (psf)

% Drag Forces
D_0 = CD0.total.sub * q_trans * WING.geom.S_area;

% Constants 
K_perf = 1 / (pi * e * WING.geom.AR);
b1 = sqrt(D_0);
b2 = sqrt((K_perf / WING.geom.S_area) * (1/q_trans));

% TSFC
Wt.trans.sfc_cr = 0.484/3600; % From Mattingly (100% Thrust, 37 kft)

% Pt. B 
Range.trans_alt.B = (V_transCr / Wt.trans.sfc_cr) * (1 / (b1 * b2)) * (atan((0.97 * Wt.WTO) * (b2/b1)) - atan(Wt.W1.B *(b2/b1)));

% Pt. C
Range.trans_alt.C = (V_transCr / Wt.trans.sfc_cr) * (1 / (b1 * b2)) * (atan((0.97 * Wt.WTO) * (b2/b1)) - atan(Wt.W1.C *(b2/b1)));

% Pt. D
Range.trans_alt.D = (V_transCr / Wt.trans.sfc_cr) * (1 / (b1 * b2)) * (atan((Wt.W0.D) * (b2/b1)) - atan(Wt.W1.D *(b2/b1)));

% Store ranges in array
Range.trans_alt.array = [Range.trans_alt.B, Range.trans_alt.C, Range.trans_alt.D];

% Convert ft to nmi
Range.trans_alt.array = Range.trans_alt.array * 0.000164579; 

%% Payload Range Diagram (Const Mach Number and Altitude)

% Plot Payload Range Diagram
figure()

%-----Supersonic:Const Mach # and Alt-----%
% pt. A
plot(0, Wt.fuel.w_tot + Wt.pld.w_tot, 's', 'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
hold on;
% pt. B
plot(Range.c_alt.array(1), Wt.fuel.w_tot + Wt.pld.w_tot,'s', 'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
% pt. C
plot(Range.c_alt.array(2), Wt.WTO - Wt.fuel.w_max, 's', 'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
% pt. D
plot(Range.c_alt.array(3), 0, 's', 'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
% Connect points with lines
hLine1 = plot([0,Range.c_alt.array(1)], [Wt.fuel.w_tot + Wt.pld.w_tot,Wt.fuel.w_tot + Wt.pld.w_tot], 'b');
plot([Range.c_alt.array(1),Range.c_alt.array(2)], [Wt.fuel.w_tot + Wt.pld.w_tot, Wt.WTO - Wt.fuel.w_max], 'b')
plot([Range.c_alt.array(2),Range.c_alt.array(3)], [Wt.WTO - Wt.fuel.w_max,0],'b')

%-----Transonic:Const Mach # and Alt-----%
% pt. A
plot(0, Wt.fuel.w_tot + Wt.pld.w_tot, 'd', 'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
hold on;
% pt. B
plot(Range.trans_alt.array(1), Wt.fuel.w_tot + Wt.pld.w_tot,'d', 'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
% pt. C
plot(Range.trans_alt.array(2), Wt.WTO - Wt.fuel.w_max, 'd', 'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
% pt. D
plot(Range.trans_alt.array(3), 0, 'd', 'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5])
% Connect points with lines
hLine2 = plot([0,Range.trans_alt.array(1)], [Wt.fuel.w_tot + Wt.pld.w_tot,Wt.fuel.w_tot + Wt.pld.w_tot], 'r');
plot([Range.trans_alt.array(1),Range.trans_alt.array(2)], [Wt.fuel.w_tot + Wt.pld.w_tot, Wt.WTO - Wt.fuel.w_max], 'r')
plot([Range.trans_alt.array(2),Range.trans_alt.array(3)], [Wt.WTO - Wt.fuel.w_max,0],'r')

title('Payload Range Diagram (Constant M and Alt.)')
xlabel('Range (nmi)')
ylabel('Payload (lb)')
legend([hLine1,hLine2],'Supersonic Cruise (M=1.6)','Transonic Cruise (M=0.8)', 'Location', 'Southwest')

%% Range at Cruise Climb Conditions (Supersonic)
% Range at cruise climb conditions (constant lift coefficient, CL, and Mach number, M) and
% assuming thrust specific fuel consumption, c, remains constant throughout the mission

% See Roskam Airplane Aero and Perf Sect 11.1.2.1 (eqns 11.10 & 11.11)
%L_D.cr_super.avg = 0.5*sqrt((pi*e*WING.geom.AR)/CD0.total.super); 
L_D.cr_super.avg = 7.5;

% Pt. B
Range.cr_climb.B = ((Wt.fuel.a_snd * 1116.5) / Wt.fuel.sfc_cr) * req.cr_M0(1) * L_D.cr_super.avg * log((0.97 * Wt.WTO) / Wt.W1.B);

% Pt. C
Range.cr_climb.C = ((Wt.fuel.a_snd * 1116.5) / Wt.fuel.sfc_cr) * req.cr_M0(1) * L_D.cr_super.avg * log((0.97 * Wt.WTO) / Wt.W1.C);

% Pt. D
Range.cr_climb.D = ((Wt.fuel.a_snd * 1116.5) / Wt.fuel.sfc_cr) * req.cr_M0(1) * L_D.cr_super.avg * log((Wt.W0.D) / Wt.W1.D);

% Store ranges in array
Range.cr_climb.array = [Range.cr_climb.B, Range.cr_climb.C, Range.cr_climb.D];

% Convert ft to nmi
Range.cr_climb.array = Range.cr_climb.array * 0.000164579; 

%% Range at Cruise Climb Conditions (Transonic M = 0.8)

% L/D for transonic cruise climb conditions 
% See Roskam Airplane Aero and Perf Sect 11.1.2.1 (eqns 11.10 & 11.11)
L_D.cr_sub.avg = 0.5*sqrt((pi*e*WING.geom.AR)/CD0.total.sub); % this is max possible, what do I actually set this to?

% Pt. B
Range.climb_trans.B = ((atm.snd_trans * 1116.5) / Wt.trans.sfc_cr) * 0.8 * L_D.cr_sub.avg * log((0.97 * Wt.WTO) / Wt.W1.B);

% Pt. C
Range.climb_trans.C = ((atm.snd_trans * 1116.5) / Wt.trans.sfc_cr) * 0.8 * L_D.cr_sub.avg  * log((0.97 * Wt.WTO) / Wt.W1.C);

% Pt. D
Range.climb_trans.D = ((atm.snd_trans * 1116.5) / Wt.trans.sfc_cr) * 0.8 * L_D.cr_sub.avg * log((Wt.W0.D) / Wt.W1.D);

% Store ranges in array
Range.climb_trans.array = [Range.climb_trans.B, Range.climb_trans.C, Range.climb_trans.D];

% Convert ft to nmi
Range.climb_trans.array = Range.climb_trans.array * 0.000164579; 

%% Payload Range Diagram (Cruise Climb Conditions)
% Cruise climb conditions indicate constant lift coefficient and Mach
% number

%-----Supersonic:Cruise Climb-----%

% Plot Payload Range Diagram
figure()
% pt. A
plot(0, Wt.fuel.w_tot + Wt.pld.w_tot, 's', 'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
hold on;
% pt. B
plot(Range.cr_climb.array(1), Wt.fuel.w_tot + Wt.pld.w_tot,'s', 'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
% pt. C
plot(Range.cr_climb.array(2), Wt.WTO - Wt.fuel.w_max, 's', 'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
% pt. D
plot(Range.cr_climb.array(3), 0, 's', 'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
% Connect points with lines
hline3 = plot([0,Range.cr_climb.array(1)], [Wt.fuel.w_tot + Wt.pld.w_tot,Wt.fuel.w_tot + Wt.pld.w_tot], 'b');
plot([Range.cr_climb.array(1),Range.cr_climb.array(2)], [Wt.fuel.w_tot + Wt.pld.w_tot, Wt.WTO - Wt.fuel.w_max], 'b')
plot([Range.cr_climb.array(2),Range.cr_climb.array(3)], [Wt.WTO - Wt.fuel.w_max,0],'b')

%-----Transonic:Cruise Climb-----%

% Plot Payload Range Diagram
% pt. A
plot(0, Wt.fuel.w_tot + Wt.pld.w_tot, 'd', 'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
hold on;
% pt. B
plot(Range.climb_trans.array(1), Wt.fuel.w_tot + Wt.pld.w_tot,'d', 'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
% pt. C
plot(Range.climb_trans.array(2), Wt.WTO - Wt.fuel.w_max, 'd', 'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
% pt. D
plot(Range.climb_trans.array(3), 0, 'd', 'MarkerSize',6,...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor',[0.5,0.5,0.5]) 
% Connect points with lines
hline4 = plot([0,Range.climb_trans.array(1)], [Wt.fuel.w_tot + Wt.pld.w_tot,Wt.fuel.w_tot + Wt.pld.w_tot],'r');
plot([Range.climb_trans.array(1),Range.climb_trans.array(2)], [Wt.fuel.w_tot + Wt.pld.w_tot, Wt.WTO - Wt.fuel.w_max], 'r')
plot([Range.climb_trans.array(2),Range.climb_trans.array(3)], [Wt.WTO - Wt.fuel.w_max,0],'r')

% Plot vertical lines for pt B and pt C
% plot([Range.climb_trans.array(1),Range.climb_trans.array(1)], [Wt.fuel.w_tot + Wt.pld.w_tot,0],'k:')
% plot([Range.climb_trans.array(2),Range.climb_trans.array(2)], [Wt.WTO - Wt.fuel.w_max,0],'k:')
% plot([Range.cr_climb.array(1),Range.cr_climb.array(1)], [Wt.fuel.w_tot + Wt.pld.w_tot,0],'k:')
% plot([Range.cr_climb.array(2),Range.cr_climb.array(2)], [Wt.WTO - Wt.fuel.w_max,0],'k:')

title('Payload Range Diagram (Cruise-Climb)')
xlabel('Range (nmi)')
ylabel('Payload (lb)')
legend([hline3, hline4],'Supersonic Cruise (M=1.6)','Transonic Cruise (M=0.8)', 'Location', 'Southwest')

