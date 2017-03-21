clc;
clear;
close all;

% load('aircraft_vars');
hld_calcs;
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
h_start = 5; % ft above ground
AR_eff = @(h) WING.geom.AR / (spline(figs.fig11(:,1)', figs.fig11(:,2)', 2*h/WING.geom.span));
alf0 = 2*pi/180; % 2 radians
CD0 = 0.0204;
N_nose = 0.25;
Nose_arm = 50; % ft from CG
% Main_arm = 
L_arm = 0.3*WING.geom.sup.Cr;

% Define some temporary functions
CL_alf = @(AR, Ma) (2*pi * AR)/(2 + sqrt(AR^2 * sqrt(abs(Ma^2 - 1)) * (1 + tand(10)^2/abs(Ma^2 - 1) + 4)));
% d_Alf_0g = @(h) (0.5*(WING.biconvex.tc_u + WING.biconvex.tc_l)) * (-0.1177/((h/WING.geom.MAC)^2) + 3.5655/((h / WING.geom.MAC)^2));
d_CDi_g = @(h, Cl, AR) -1* (1 - 1.32*(h/WING.geom.span))/(1.05 + 7.4*(h/WING.geom.span)) * Cl^2 / (pi * AR);
L_g = @(Cl, V) 0.5*rho*WING.geom.S_area*V^2 * Cl;
D_g = @(Cd, V) 0.5*rho*WING.geom.S_area*V^2 * Cd;
d_Alf_0g = @(h)(WING.biconvex(1).tc_u + WING.biconvex(1).tc_l) * (-0.1177/((h/WING.geom.MAC)^2) + 3.5655/((h / WING.geom.MAC)^2)) * pi / 180;

% Initialize boundary conditions
dt = 0.05; % seconds
Xpos(1,:) = [0.0, 0.0, 0.0]; % [X0, V0, t]
ii = 1; % iterate with this
L_i = 0.0;

while (L_i < Wt.WTO)
    ii = ii + 1;
    
    % Lift coefficient model
    CL_ag = CL_alf(AR_eff(h_start), Xpos(ii-1,2)/a0); % Ground effect lift coefficient
    dAlf_0g_i = d_Alf_0g(h_start);
    
    CLi = CL_ag*(alf0-dAlf_0g_i) + dCL; 
    % missing dCl for flaps at various AOA
    % Also missing loss of ground effect due to flaps
    L_i = L_g(CLi, Xpos(ii-1, 2));
    
    % Drag Coefficient Model
    Cdi = CD0 + CLi^2 / (pi * AR_eff(h_start) * 0.9) + d_CDi_g(h_start, CLi, AR_eff(h_start));
    D_i = D_g(Cdi, Xpos(ii-1, 2));
    
    % Ground friction force
    if Wt.WTO > L_i
        GF = Wt.WTO - L_i;
    else
        GF = 0; % want zero if not touching the ground 
    end
    
    dV = dt*(gc/Wt.WTO)*(Thr - D_i - (mu_g * GF) - Wt.WTO*sin(phi));
    V_ii = Xpos(ii-1,2) + dV;
    t_ii = dt + Xpos(ii-1, 3);
    X_ii = Xpos(ii-1, 1) + dt * Xpos(ii-1, 2);
    
    Xpos(ii, :) = [X_ii, V_ii, t_ii];
    
    if (L_i * L_arm > (N_nose * Wt.WTO * Nose_arm)) || (L_i > Wt.WTO)
        fprintf('Nose Lifted for Take-Off at Time %0.5f\n', Xpos(ii, 3));
    end
end

fprintf('Take-Off Ground Roll Distance: %0.3f ft\n', Xpos(end, 1));
fprintf('V_{LOF}: %0.3f ft/s\n', Xpos(end, 2));

%% Takeoff Transition - bring vehicle to the right gamma

gamma = 0.0;
Thr = constraints.req_Thr*0.8; % ASSUMPTION HERE... CONSTANT THRUST! 

% Change alpha scheduling to clear obstacle... 1deg / 1 sec
alpha_schedule = 0:(dt*pi/180):hld.alpha_TO;

dS_trans = 0.0; % calculate incremental distance for transition period
h_trans = h_start;

iii = ii + 1;

if (iii - ii) > length(alpha_schedule)
    alpha_i1 = alpha_schedule(end);
else%if %h_trans >= 7
    alpha_i1 = alpha_schedule(iii-ii);
% else
%     alpha_i1 = 0;
end

% Lift coefficient model
CL_ag = CL_alf(AR_eff(h_trans), Xpos(iii-1,2)/a0); % Ground effect lift coefficient
dAlf_0g_i = d_Alf_0g(h_trans);

CLi = CL_ag*(alpha_i1 + alf0 -dAlf_0g_i) + dCL; 
% missing dCl for flaps at various AOA
% Also missing loss of ground effect due to flaps
L_i = L_g(CLi, Xpos(iii-1, 2));

% Drag Coefficient Model
Cdi = CD0 + CLi^2 / (pi * AR_eff(h_trans) * 0.9) + d_CDi_g(h_trans, CLi, AR_eff(h_trans));
D_i = D_g(Cdi, Xpos(iii-1, 2));

gamma(2) = gc / (Xpos(iii - 1,2)*Wt.WTO) * (Thr * sin(alpha_i1) + L_i - Wt.WTO*cos(gamma(1)));
dV = dt*(gc/Wt.WTO)*(Thr*cos(alpha_i1) - D_i - Wt.WTO*sin(gamma(1)));
V_ii = Xpos(iii-1,2) + dV;
t_ii = dt + Xpos(iii-1, 3);
X_ii = Xpos(iii-1, 1) + dt * Xpos(iii-1, 2);

Xpos(iii, :) = [X_ii, V_ii, t_ii];
h_trans = h_trans + (dt*V_ii)*sin(gamma(2));
dS_trans = dS_trans + (dt*V_ii)*cos(gamma(2));
cleared = 0;

while gamma(iii - ii) < hld.gamma_TO
    iii = iii + 1;
    
    if (iii - ii) > length(alpha_schedule) && cleared
        alpha_i1 = 0.0;
    elseif (iii - ii) > length(alpha_schedule)
        alpha_i1 = alpha_schedule(end);
    else
        alpha_i1 = alpha_schedule(iii-ii);
    end
    
    % Lift coefficient model
    CL_ag = CL_alf(AR_eff(h_trans), Xpos(iii-1,2)/a0); % Ground effect lift coefficient
    dAlf_0g_i = d_Alf_0g(h_trans);
    
    CLi = CL_ag*(alpha_i1 + alf0 -dAlf_0g_i) + dCL; 
    % missing dCl for flaps at various AOA
    % Also missing loss of ground effect due to flaps
    L_i = L_g(CLi, Xpos(iii-1, 2));
    
    % Drag Coefficient Model
    Cdi = CD0 + CLi^2 / (pi * AR_eff(h_trans) * 0.9) + d_CDi_g(h_trans, CLi, AR_eff(h_trans));
    D_i = D_g(Cdi, Xpos(iii-1, 2));
    
    gamma(iii - ii) = gc / (Xpos(iii - 1,2)*Wt.WTO) * (Thr * sin(alpha_i1) + L_i - Wt.WTO*cos(gamma(iii-ii-1)));
    dV = dt*(gc/Wt.WTO)*(Thr*cos(alpha_i1) - D_i - Wt.WTO*sin(gamma(iii-ii-1)));
    V_ii = Xpos(iii-1,2) + dV;
    t_ii = dt + Xpos(iii-1, 3);
    X_ii = Xpos(iii-1, 1) + dt * Xpos(iii-1, 2);
    
    Xpos(iii, :) = [X_ii, V_ii, t_ii];
    
    h_trans = h_trans + (dt*V_ii)*sin(gamma(iii-ii));
    dS_trans = dS_trans + (dt*V_ii)*cos(gamma(iii-ii));
    
    if h_trans > 35 && ((h_trans - (dt*V_ii)*sin(gamma(iii-ii))) < 35)
       fprintf('Vehicle clears 35ft screen at Time: %0.5f\n', Xpos(iii,3));
       XTO = Xpos(ii,1) + dS_trans;
       fprintf('Total Take-off Distance Covered: %0.5f ft\n', Xpos(ii,1) + dS_trans);
       cleared = 1;
    end
    
end

fprintf('Time to Climb-Out: %0.5f\n', Xpos(end,3) - Xpos(ii, 3));
fprintf('V_INCL: %0.5f\n', Xpos(end, 2));
fprintf('Transition Distance: %0.5f\n', dS_trans);

