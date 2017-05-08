clc;
clear;

load('aircraft_vars');
hld_calcs;
Drag_Analysis;
close all;
addpath('TO_TD_perf\');

%% Notes (from Roskam and Lam)

% takeoff ends when aircraft reaches height of 1500ft above terrain (climb
% to cruise occurs here)

% Important reference parameters
    % V_s = calibrated stall speed or minimum steady flight speed
    % V_MCG = min grond ctrl speed
    % V_EF = engine-failure speed (V_EF > V_MCG)
    % V_MC = min ctrl air speed (FAR 25.149), may not exceed 1.2 V_s
    % V_1 = takeoff decision speed; takeoff must continued even if crit.
        % eng failure occurs, (V_1 > V_EF > V_MCG)
    % V_R = speed at which pilot initiates roatation (V_R >= V_1, V_R >= 1.05V_MC)
    % V_MU = speed to continue with lift off and takeoff (V_MU >= V_2min, V_MU >= V_R)
    % V_LOF = calibrated speed at which airplane becomes airborne (V_LOF > V_R)
    % V_2min = min safety speed (V_2min >= 1.2 V_s [we are twin-engined?])
    % V_2 = provides climb gradient required by FAR 25.121 (V_2 >= V_R + sped increment)

%% Equations of Motion During Take-Off

% hld.gamma_TO = atan(35/100); % FAR 25 requirement for 35 ft obstacle, and with 100 flight path to clear it
% hld.alpha_TO = 10 * pi / 180; % take-off AOA, deg -> radians
% hld.V_TO = constraints.Vstall * 1.2; % fly at 20% safety margin 
% hld.Cl_base = 0.535; % generated w/ OpenVSP
% hld.sf = 1.2; % safety factor to obtain Cl_max
% hld.dCl_req = (hld.sf*Wt.WTO * cos(hld.gamma_TO) - constraints.req_Thr * sin(hld.alpha_TO))/(0.5*atm.rho_sl*(hld.V_TO^2)*WING.geom.S_area) - hld.Cl_base;
% fprintf('Cl deficit: %0.5f\n', hld.dCl_req);

% Takeoff Equation:
% S_TO = S_G + S_A = S_NGR + S_R + S_TR + S_CL

%% Take-off Ground Roll (S_NGR + S_R)

fprintf('\nTakeoff Calculations\n');

% Flight Conditions
P_ref = 2116; % lbf/ft^2
T_ref = 518.69; % Rankine
a0 = 1116; % ft/s
rho = 0.0023769; % slugs/ft^2
gc = 32.174;
mu_g = 0.25; % typical constant for asphalt/tarmac
phi = 0.0; % runway gradient

Thr = constraints.req_Thr; % ASSUMPTION HERE... CONSTANT THRUST! 

% Figure 10.11 on Roskam and Lam, x-axis is 2h/b, y-axis is AR/AReff
figs.fig11 = dlmread('effective_AR_roskam_lam_fig10-11.txt', '%0.5f\t%0.5f\n', 2, 0);
figs.fig11 = figs.fig11(find([diff(figs.fig11(:,1)); 1]),:);

% Aircraft Properties
h_start = 5.54; % ft above ground
AR_eff = @(h) WING.geom.AR / (spline(figs.fig11(:,1)', figs.fig11(:,2)', 2*h/WING.geom.span));
alf0 = 2*pi/180; % 2 radians
CD0_TO = CD0.total.TO;
N_nose = 0.25;
Nose_arm = 50; % ft from CG
% Main_arm = 
L_arm = 0.3*WING.geom.sup.Cr;

% Define some temporary functions
CL_alf = @(AR, Ma) (2.*pi .* AR)./(2 + sqrt(AR^2 .* sqrt(abs(Ma.^2 - 1)) .* (1 + (tand(10)^2)./abs(Ma.^2 - 1) + 4)));
% d_Alf_0g = @(h) (0.5*(WING.biconvex.tc_u + WING.biconvex.tc_l)) * (-0.1177/((h/WING.geom.MAC)^2) + 3.5655/((h / WING.geom.MAC)^2));
d_CDi_g = @(h, Cl, AR) -1.* (1 - 1.32.*(h./WING.geom.span))./(1.05 + 7.4.*(h./WING.geom.span)) .* Cl.^2 ./ (pi .* AR);
L_g = @(Cl, V) 0.5*rho*WING.geom.S_area*V^2 * Cl;
D_g = @(Cd, V) 0.5*rho*WING.geom.S_area*V^2 * Cd;
d_Alf_0g = @(h)(WING.biconvex(1).tc_u + WING.biconvex(1).tc_l) * (-0.1177/((h/WING.geom.MAC)^2) + 3.5655/((h / WING.geom.MAC)^2)) * pi / 180;

% Initialize boundary conditions
% dt = 0.05; % seconds
ii = 1; % iterate with this
dV_ii = 0.1;
V_TO = (1:3).*dV_ii;
Vw = 0 * 1.467; % mph -> ft/s, constant head-wind speed

% Lift  
CL_ag = CL_alf(AR_eff(h_start), V_TO./a0); % Ground effect lift coefficient
dAlf_0g_i = d_Alf_0g(h_start);
CLg = CL_ag*(alf0-dAlf_0g_i) + dCL_tot; 

% Drag 
Cdg = CD0_TO + (CLg.^2) ./ (pi * AR_eff(h_start) * 0.9) + d_CDi_g(h_start, CLg, AR_eff(h_start));

q_bar = (0.5*atm.rho_sl.*(V_TO(ii:(ii+2))+Vw).^2); 
Ag_TO = gc.*((Thr/Wt.WTO - mu_g) - (Cdg - mu_g .* CLg).*q_bar./(Wt.WTO./WING.geom.S_area) - phi);
f_ii = (V_TO)./Ag_TO;
Sg_TO(ii) = (dV_ii./3).*(f_ii(ii) + 4*f_ii(ii+1) + f_ii(ii+2));

while (L_g(CLg(end), V_TO(end)+Vw) < Wt.WTO)% || (Xpos(ii,2) < 1.2*constraints.Vstall)
    % renew consistent variables
    ii = ii + 1;
    V_TO(end+1) = V_TO(end) + dV_ii;
    
    % Lift coefficient model
    CL_ag = CL_alf(AR_eff(h_start), (V_TO(ii:(ii+2))+Vw)./a0); % Ground effect lift coefficient
    dAlf_0g_i = d_Alf_0g(h_start);
    CLg = CL_ag*(alf0-dAlf_0g_i) + dCL_tot; 
    
    % missing dCl for flaps at various AOA
    % Also missing loss of ground effect due to flaps
    
    % Drag Coefficient Model
    Cdg = CD0_TO + (CLg.^2) ./ (pi * AR_eff(h_start) * 0.9) + d_CDi_g(h_start, CLg, AR_eff(h_start));
%     D_i = D_g(Cdi, Xpos(ii-1, 2) + 1.5*Vw);
    
    q_bar = (0.5*atm.rho_sl.*(V_TO(ii:(ii+2))+Vw).^2); % dynamic pressure should be greater if there is headwind
    Ag_TO = gc.*((Thr/Wt.WTO - mu_g) - (Cdg - mu_g .* CLg).*q_bar./(Wt.WTO./WING.geom.S_area) - phi);
    f_ii = (V_TO(ii:(ii+2)))./Ag_TO;
    Sg_TO(ii) = (dV_ii./3).*(f_ii(1) + 4*f_ii(2) + f_ii(3)) + Sg_TO(ii-1);

end

V_LOF = V_TO(end);

fprintf('Take-Off Ground Roll Distance: %0.3f ft\n', Sg_TO(end));
fprintf('V_{LOF}: %0.3f ft/s\n', V_LOF);

%% Takeoff Transition - bring vehicle to the right gamma

% gamma = 0.0;
% Thr = constraints.req_Thr*0.8; % ASSUMPTION HERE... CONSTANT THRUST! 
% 
% % Change alpha scheduling to clear obstacle... 1deg / 1 sec
% alpha_schedule = 0:(dt*pi/180):hld.alpha_TO;
% 
% dS_trans = 0.0; % calculate incremental distance for transition period
% h_trans = h_start;
% 
% iii = ii + 1;
% 
% if (iii - ii) > length(alpha_schedule)
%     alpha_i1 = alpha_schedule(end);
% else
%     alpha_i1 = alpha_schedule(iii-ii);
% % else
% %     alpha_i1 = 0;
% end
% 
% % Lift coefficient model
% CL_ag = CL_alf(AR_eff(h_trans), Xpos(iii-1,2)/a0); % Ground effect lift coefficient
% dAlf_0g_i = d_Alf_0g(h_trans);
% 
% CLg = CL_ag*(alpha_i1 + alf0 -dAlf_0g_i) + dCL; 
% % missing dCl for flaps at various AOA
% % Also missing loss of ground effect due to flaps
% L_i = L_g(CLg, Xpos(iii-1, 2));
% 
% % Drag Coefficient Model
% Cdg = CD0_TO + CLg^2 / (pi * AR_eff(h_trans) * 0.9) + d_CDi_g(h_trans, CLg, AR_eff(h_trans));
% D_i = D_g(Cdg, Xpos(iii-1, 2));
% 
% gamma(2) = gc / (Xpos(iii - 1,2)*Wt.WTO) * (Thr * sin(alpha_i1) + L_i - Wt.WTO*cos(gamma(1)));
% dV = dt*(gc/Wt.WTO)*(Thr*cos(alpha_i1) - D_i - Wt.WTO*sin(gamma(1)));
% V_ii = Xpos(iii-1,2) + dV;
% t_ii = dt + Xpos(iii-1, 3);
% X_ii = Xpos(iii-1, 1) + dt * Xpos(iii-1, 2);
% 
% Xpos(iii, :) = [X_ii, V_ii, t_ii];
% h_trans = h_trans + (dt*V_ii)*sin(gamma(2));
% dS_trans = dS_trans + (dt*V_ii)*cos(gamma(2));
% cleared = 0;
% ii_a = 1;
% 
% while h_trans < 35 
%     iii = iii + 1;
%     
%     if (iii - ii) > length(alpha_schedule) && cleared
%         alpha_i1 = 0.0;
%     elseif (iii - ii) > length(alpha_schedule)
%         alpha_i1 = alpha_schedule(end);
%     else
%         alpha_i1 = alpha_schedule(ii_a);
%         ii_a = ii_a + 1;
%     end
%     
%     % Lift coefficient model
%     CL_ag = CL_alf(AR_eff(h_trans), Xpos(iii-1,2)/a0); % Ground effect lift coefficient
%     dAlf_0g_i = d_Alf_0g(h_trans);
%     
%     CLg = CL_ag*((alpha_i1-gamma(iii-ii-1)) + alf0 -dAlf_0g_i) + dCL; 
%     % missing dCl for flaps at various AOA
%     % Also missing loss of ground effect due to flaps
%     L_i = L_g(CLg, Xpos(iii-1, 2));
%     
%     % Drag Coefficient Model
%     Cdg = CD0_TO + CLg^2 / (pi * AR_eff(h_trans) * 0.9) + d_CDi_g(h_trans, CLg, AR_eff(h_trans));
%     D_i = D_g(Cdg, Xpos(iii-1, 2));
%     
%     gamma(iii - ii) = gc / (Xpos(iii - 1,2)*Wt.WTO) * (Thr * sin((alpha_i1-gamma(iii-ii-1))) + L_i - Wt.WTO*cos(gamma(iii-ii-1)));
%     dV = dt*(gc/Wt.WTO)*(Thr*cos((alpha_i1-gamma(iii-ii-1))) - D_i - Wt.WTO*sin(gamma(iii-ii-1)));
%     V_ii = Xpos(iii-1,2) + dV;
%     t_ii = dt + Xpos(iii-1, 3);
%     X_ii = Xpos(iii-1, 1) + dt * Xpos(iii-1, 2);
%     
%     Xpos(iii, :) = [X_ii, V_ii, t_ii];
%     
%     h_trans = h_trans + (dt*V_ii)*sin(gamma(iii-ii));
%     dS_trans = dS_trans + (dt*V_ii)*cos(gamma(iii-ii));
%     
%     if h_trans > 35 && ((h_trans - (dt*V_ii)*sin(gamma(iii-ii))) < 35)
%        fprintf('Vehicle clears 35ft screen at Time: %0.5f\n', Xpos(iii,3));
%        XTO = Xpos(ii,1) + dS_trans;
%        fprintf('Total Take-off Distance Covered: %0.5f ft\n', Xpos(ii,1) + dS_trans);
%        cleared = 1;
%     end
%     
%     if (iii - ii - 1) == 1
%        fprintf('Lift-off Flight Path Angle: %0.5f deg \n', gamma(iii-ii)*180/pi); 
%        fprintf('Pitch Angle: %0.5f deg\n', (gamma(iii-ii) + alpha_i1)*180/pi);
%        fprintf('Wing Height: %0.3f ft\n', h_trans);
%     end
%     
% end
% 
% fprintf('Time to Climb-Out: %0.5f\n', Xpos(end,3) - Xpos(ii, 3));
% fprintf('V_INCL: %0.5f\n', Xpos(end, 2));
% fprintf('Transition Distance: %0.5f\n', dS_trans);
% fprintf('Final flight path angle: %0.5f deg\n', gamma(end)*180/pi);


%% Landing Calcs -> Approach Air Distance

fprintf('\nLanding Calculations\n');
% Requirement velocities according to FAR 25
V_A = 1.3*constraints.Vstall;
% V_TD = 1.1 * V_A;
V_FL = 0.95*V_A;
gam_A = -3.0*pi/180; % degrees, approach angle
h_screen = 50;
n_FL = 1.04;

% Approach Calculations
WT_land = Wt.WTO * Wt.fuel.w2_1 * Wt.fuel.w3_2 * Wt.fuel.w4_3 * Wt.fuel.w5_4 * Wt.fuel.w6_5;
CL_A = 2* WT_land /(rho * WING.geom.S_area * V_A^2);
CD_A = CD0.total.Land + CL_A^2 / (pi * WING.geom.AR * 0.9);
Thr_land = (CD_A/CL_A - gam_A ) * WT_land;
R_flare = V_FL^2 / (gc * (n_FL - 1));
h_flare = R_flare * (1 - cos(gam_A));

S_LTR = -R_flare * gam_A;
S_LDES = max((h_screen - h_flare)/tan(gam_A), 2*h_flare / gam_A);

V_TD = sqrt(V_FL^2 + ((2*gc)/WT_land)*(WT_land*h_flare - 0.5*WT_land*(gam_A + (CD_A/CL_A))*S_LTR)); % error in Roskam equation

fprintf('Required Landing Thrust: %0.2f lbf\n', Thr_land);
fprintf('Approach Velocity: %0.2f ft/s\n', V_A);
fprintf('Flare Velocity: %0.2f ft/s\n', V_FL); 
fprintf('Touchdown Velocity: %0.2f ft/s\n', V_TD);
fprintf('Approach Angle: %0.2f deg\n', -gam_A*180/pi);
fprintf('Total Air Landing Distance: %0.2f ft\n\n', S_LTR + S_LDES);

%% Landing -> Ground Roll Calculations

X_LR = [0.0, V_TD, 0.0];
dt = 0.1;
ii = 1;
mu_gBrakes = 0.5*(0.65+0.8);
Vw = 0 * 1.467; % mph -> ft/s, constant head-wind speed
N_n = 0.0;

while X_LR(end, 2) > 0
    ii = ii + 1;
    
    % Lift coefficient model
    CL_ag = CL_alf(AR_eff(h_start), (X_LR(ii-1,2)+Vw)/a0); % Ground effect lift coefficient
    dAlf_0g_i = d_Alf_0g(h_start);
    
    CLg = CL_ag*(alf0-dAlf_0g_i) + dCL_tot; 
    L_i = L_g(CLg, (X_LR(ii-1,2)+Vw));
    
    % Drag Coefficient Model
    Cdg = CD0_TO + CLg^2 / (pi * AR_eff(h_start) * 0.9) + d_CDi_g(h_start, CLg, AR_eff(h_start));
    D_i = D_g(Cdg, (X_LR(ii-1,2)+Vw));
    
%     if (L_g(CLg(1), (X_LR(ii-1,2)+Vw)) < Wt.WTO) && (N_n == 0)
%        N_n = 0.18*WT_land; % switch to activate LR sequence
%        V_LR = X_LR(ii-1,2);
%        fprintf('LR Velocity: %0.3f ft/s\n', V_LR);
%        fprintf('LR Distance: %0.3f ft\n', Sg_TD(end));
%        fprintf('LR Decceleration: %0.3f ft/s^2\n', Ag_TD(ii));
%     end
    
    if L_i > WT_land
        GF = 0;
    else
        GF = WT_land - L_i; 
    end
    
    dV_ii = dt*(gc/WT_land)*( - D_i - (mu_gBrakes)*(GF) - WT_land*sin(phi) + N_n*(mu_gBrakes - mu_g));
    V_ii = X_LR(ii-1,2) + dV_ii;
    t_ii = dt + X_LR(ii-1, 3);
    X_ii = X_LR(ii-1, 1) + dt * X_LR(ii-1, 2);
    
    X_LR(ii, :) = [X_ii, V_ii, t_ii];
    
    if (L_i < WT_land) && (N_n == 0)
        fprintf('Main Landed for Touchdown at Time: %0.2fs\n', X_LR(ii, 3));
        fprintf('Main Landing Position from TD Point: %0.2f ft\n', X_LR(ii,1));
        S_LR = X_LR(ii,1);
        fprintf('Actual Touchdown Velocity: %0.2f ft/s\n', X_LR(ii, 2));
        V_LR = X_LR(ii, 2);
        N_n = 0.08 * WT_land;
    end
end

fprintf('Actual Ground Landing Distance: %0.2f ft\n', X_LR(end,1)-S_LR);
S_L_tot = X_LR(end,1) + S_LTR + S_LDES;
fprintf('Total Landing Distance: %0.2f ft\n', S_L_tot);

%% Engine Out Analysis

addpath(genpath('power_curves\'));
Eout_data = importdata('Engine_out_landing.xlsx');
thrust_data.eout = struct('mach', Eout_data.data(:,1), 'thr', Eout_data.data(:,2).*2);

V_inf = (0.01*a0):0.01:(V_LOF*1.5);
THR_climb = spline(thrust_data.eout.mach, thrust_data.eout.thr, V_inf./a0);
dAlf_0g_i = d_Alf_0g(h_start);
CL_eo = CL_alf(AR_eff(h_start), V_inf./a0).*(alf0-dAlf_0g_i) + dCL_tot;

Cd_eo = CD0_TO + (CL_eo.^2) / (pi * AR_eff(h_start) * 0.9) + d_CDi_g(h_start, CL_eo, AR_eff(h_start));

P_avail = THR_climb .* V_inf;
P_req = WT_land .* V_inf .* Cd_eo ./ CL_eo;
RC = ((P_avail - P_req)./Wt.WTO);
V_RCmax = V_inf(RC == max(RC));

figure();plot(V_inf, P_req);
hold on;
plot(V_inf, P_avail);
legend('P_{req}', 'P_{avail}', 'Location', 'Best');
title('Aborted Landing (1 Engine Out) - Sea Level, Std. Day');
xlabel('ft/s');
ylabel('ft \cdot lbf');
plot([V_TD, V_TD], [min(ylim), max(P_avail)], 'linestyle', ':', 'color', 'k');
strmax = sprintf('V_{Touchdown} = %0.2f ft/s\nRC = %0.2f fpm',V_TD, RC(sum(V_inf <= V_TD))*60);
text(V_inf(sum(V_inf<= V_TD)), mean(ylim), strmax, 'HorizontalAlignment', 'left');
plot([V_RCmax, V_RCmax], [P_req(RC == max(RC)), P_avail(RC == max(RC))], 'linestyle', '--', 'color', 'k'); 
strmax = sprintf('\tMax RC = %0.2f fpm',(max(RC)*60));
text(V_inf(RC== max(RC)), 0.5*mean(ylim), strmax, 'HorizontalAlignment', 'right');
saveas(gcf, 'eout_climb_curve.png');