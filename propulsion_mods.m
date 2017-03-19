clc;
clear;
close all;

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
eta_b = 0.99;
e_c = 0.9;
e_t = 0.9;
pi_b = 0.95;
pi_n = 0.98;
pi_d = 0.96;

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

area = pi * 0.25 * (54/12)^2;
rho = 0.000531556; % slugs/ft^3
M0 = 1.6;
mdot_alt = 355.787 * 1.6^2 * area * gamma_c / (1.6 * 968.076);