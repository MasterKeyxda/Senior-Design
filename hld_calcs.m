clear;
clc;
close all;

load('aircraft_vars.mat');

%% FLAP/SLAT CALCULATIONS FOR TAKEOFF

% Will be using fowler flap and slats to get a wing that will allow us to
% takeoff

%% LIFT REQUIREMENTS & DEFICIT
% Determine the lift coefficient needed

% WEIGHT: Wt.WTO
% REQUIRED THRUST: constraints.req_Thr

hld.gamma_TO = atan(35/100); % FAR 25 requirement for 35 ft obstacle, and with 100 flight path to clear it
hld.alpha_TO = 10 * pi / 180; % take-off AOA, deg -> radians
hld.V_TO = constraints.Vstall * 1.2; % fly at 20% safety margin 
hld.Cl_base = 0.535; % generated w/ OpenVSP
hld.sf = 1.2; % safety factor to obtain Cl_max
hld.dCl_req = (hld.sf*Wt.WTO * cos(hld.gamma_TO) - constraints.req_Thr * sin(hld.alpha_TO))/(0.5*atm.rho_sl*(hld.V_TO^2)*WING.geom.S_area) - hld.Cl_base;
fprintf('Cl deficit: %0.5f\n', hld.dCl_req);

%% DETERMINE THE HLD CONFIGURATION 