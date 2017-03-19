clc;
clear;
close all;

load('aircraft_vars.mat');

%% Reverse calculate the pi_c and pi_t

% GIVENS
mdot = 498; % lb/s
SFC = 0.56;
Tt4 = 2209.67; % R
Tt0 = 518.69; % R
Force = 21e3; %lbf
BPR = 1.21; % bypass ratio
OPR = 16;

% Air Properties
cp_c = 0.24; %BTU/lbm R
cp_tAIR = 0.295;
gamma_c = 1.4;

% Kerosene Properties
cp_kerosene = 0.48;
hPR_kerosene = 18300; % BTU/lbm

% Gasoline Properties
cp_gas = 0.53;
hPR_gas = 18443.6801; % BTU/lbm

% TL level 4
eta_b = 0.999;
e_c = 0.9;
e_f = 0.89;
e_t = 0.9; % uncooled
pi_b = 0.95;
pi_n = 0.98;
pi_d = 0.96; % type C
eta_m = 0.995;  % shaft only

f_ratio = SFC/3600 * (1+BPR) * Force / mdot;

% Get cp_turbine
k_ratio = 0.3;
g_ratio = 1 - k_ratio;
cp_t = cp_tAIR;%(cp_tAIR + (g_ratio * cp_gas + k_ratio * cp_kerosene)*f_ratio)/(1 + f_ratio);
R = 54.09/ 778.16; % BTU/lbm
gamma_t = 1 / (1 - R/cp_t);
fprintf('Turbine cp: %0.5f\n', cp_t);

% Get Tt3
hPR_fuel = (k_ratio * hPR_kerosene + g_ratio*hPR_gas);
Tt3 = (cp_t / cp_c) * Tt4 - (f_ratio/cp_c)*(eta_b*hPR_fuel - cp_t * Tt4);
fprintf('Tt3 = %4.5f\n', Tt3);
tau_c = Tt3 / Tt0;
pi_c = tau_c ^ (gamma_c * e_c / (gamma_c - 1));
fprintf('pi_c: %0.5f\n', pi_c); 

% pi_t and tau_t
pi_t = OPR/ ( pi_c * pi_b * pi_n * pi_d);
tau_t = pi_t ^ ((gamma_t - 1) * e_t/gamma_t);
fprintf('pi_t: %0.5f\n', pi_t);

% Fan performance (Tau_f)
tau_f = 1 + (1/BPR) * ((tau_t - 1)* (-(cp_t*Tt4/(cp_c * Tt0))*eta_m * (1 + f_ratio)) - (tau_c - 1));
pi_f = (gamma_c * e_f / (gamma_c - 1));
fprintf('pi_f: %0.5f\n', pi_f);

% area = pi * 0.25 * (54/12)^2;
% rho = 0.000531556; % slugs/ft^3
% press = 355.787; % lbf/ft^2
% M0 = 1.6;
% mdot_alt = rho * area * (1.6 * 968.076) * 32.174;
% mdot_alt = 

%% Get New Engine Power Curves?

% FC: flight conditions (P0, T0, M0, gamma, cp)
% ENG: engine parameters (pi_c, pi_f, Tt4_max, hPR, ENG.p0_9, ENG.p0_19, d (diameter, ft), BPR)
% THR: required thrust
% TL: tech level

FC.P0 = 2116.23*atm.delta; % lbf/ft^2
FC.T0 = 518.67 * atm.theta; % Rankine
FC.M0 = req.cr_M0(1);
FC.gamma = 1.4;  % air
FC.cp = 0.24; % BTU/lbm R
FC.rho = 0.000531556; % slugs/ft^3

ENG.pi_c = pi_c;
ENG.pi_f = pi_f;
ENG.Tt4_sl = Tt4;
ENG.T0 = Tt0;
ENG.Tt4_max = 4000; % TECH LEVEL 5
ENG.Tt4_min = Tt4;
ENG.hPR = hPR_kerosene; % replace Jet Fuel B with Jet Fuel A
ENG.p0_9 = 1.0;
ENG.p0_19 = 1.0;
ENG.d = 54/12; % in -> ft
ENG.BPR = BPR;
ENG.mdot_sl = mdot;


n_eng = 1; % three engines
thrust = (Wt.WTO * Wt.fuel.w2_1 * Wt.fuel.w3_2) / (LD_Max * n_eng);

% TL level 5
TL.d_max = 0.97; % supersonic aircraft
TL.e_c = 0.91;
TL.e_f = 0.92;
TL.pi_b = 0.96;
TL.eta_b = 0.999;
TL.e_t = 0.91; % uncooled
TL.pi_n = 0.99; % variable area convergent/divergent nozzle
TL.eta_m = 0.996; % shaft only
% TL

[Tts4, f_ratio_T, SFC_T] = throttle_calcs(FC, ENG, TL, thrust);

