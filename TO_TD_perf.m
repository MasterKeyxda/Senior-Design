clc;
clear;
close all;

load('aircraft_vars');

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

    