clc;
clear;
close all;

addpath(genpath('power_curves\'));

%% Load Data

temp_sub = importdata('DifferentAltitude-FullThrottle-stdDay-0.8M.xlsx');
temp_sup = importdata('AltitudeData-FullThrottle-stdDay.xlsx');
temp_climb = importdata('SpeedData-FullThrottle-stdDay-42ft.xlsx');

thrust_data = struct('sub', [], 'sup', [], 'climb', []);
thrust_data.sub = struct('alt', temp_sub.data(:,1), 'thr', temp_sub.data(:,2).*3);
thrust_data.sup = struct('alt', temp_sup.data(:,1), 'thr', temp_sup.data(:,2).*3);
thrust_data.climb = struct('mach', temp_climb.data(:,1), 'thr', temp_climb.data(:,2).*3);

clear temp_sub temp_sup
load('aircraft_vars.mat');
Drag_Analysis;
close all;
clearvars -except WING CD0 thrust_data atm Wt

%% Total Drag Calculations

% Need V-stall
M_min = 1.0; % set mach number to transonic speed for supersonic climb
M_num = M_min:0.01:1.6;
V_inf = M_num.*1116.*sqrt(atm.theta);
C_L = 2* Wt.WTO./(atm.rho_sl * atm.sig_rho .*(V_inf).^2 * WING.geom.S_area);

CDi = C_L.^2 ./ (0.7 * WING.geom.AR * pi);
CD_tot = CD0.total.super + CDi;
THR_climb = spline(thrust_data.climb.mach, thrust_data.climb.thr, M_num);

P_avail = THR_climb .* V_inf;
P_req = Wt.WTO .* V_inf .* CD_tot ./ C_L;

RC = ((P_avail - P_req)./Wt.WTO);
V_RCmax = V_inf(RC == max(RC));

figure();plot(V_inf, P_req./(10^6),'LineWidth',2);
hold on;
plot(V_inf, P_avail./(10^6),'LineWidth',2);

set(gca,'FontSize',18);
title('Climb Curve at 42kft (Std. Day)','FontSize', 24);
xlabel('ft/s','FontSize', 24);
ylabel('10^6 ft \cdot lbf / s', 'FontSize', 24);
legend('P_{req}', 'P_{avail}', 'Location', 'Best');
<<<<<<< HEAD
plot([V_inf(1), V_inf(1)], [min(ylim), max(P_avail./(10^6))], 'linestyle', ':', 'color', 'k');
plot([V_RCmax, V_RCmax], [P_req(RC == max(RC))./(10^6), P_avail(RC == max(RC))./(10^6)], 'linestyle', '--', 'color', 'k'); 
=======
title('Climb Curve at 42kft (Std. Day)');
xlabel('ft/s');
ylabel('ft \cdot lbf/s');
plot([V_inf(1), V_inf(1)], [min(ylim), max(P_avail)], 'linestyle', ':', 'color', 'k');
plot([V_RCmax, V_RCmax], [P_req(RC == max(RC)), P_avail(RC == max(RC))], 'linestyle', '--', 'color', 'k'); 
>>>>>>> refs/remotes/origin/master
strmax = sprintf('\tMax RC = %0.2f fpm',(max(RC)*60));
text(V_inf(RC== max(RC)), 0.85*max(ylim), strmax, 'HorizontalAlignment', 'center', 'FontSize', 18);
saveas(gcf, 'climb_curve_42.png');

%% Thrust v. Altitude Plots

figure();
plot(thrust_data.sub.alt./1000, thrust_data.sub.thr);
hold on;
plot(thrust_data.sup.alt./1000, thrust_data.sup.thr);
title('Thrust vs. Altitude');
ylabel('Thrust (lbf)');
xlabel('Alt (kft)');
legend('Subsonic Cruise (Mach 0.8)', 'Supersonic Cruise (Mach 1.6)');
saveas(gcf, 'alt_thr.png');
