clc;
clear;
close all;

load('aircraft_vars.mat');

%% Reverse calculate the pi_c and pi_t

% GIVENS - JT8D-7B
mdot = 248; % lb/s
SFC = 0.74;
Tt4 = 2700; % R
Tt0 = 518.69; % R
Force = 17.8e3; %lbf
BPR = 0.36; % bypass ratio
FPR = 3.8; % Fan Pressure Ratio
OPR = 16.0;

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

% TL level 2
eta_b = 0.94;
e_c = 0.84;
e_f = 0.82;
e_t = 0.89; % uncooled
pi_b = 0.92;
pi_n = 0.97;
pi_d = 0.95; % type C
eta_m = 0.97;  % shaft only

f_ratio = SFC/3600 * (1+BPR) * Force / mdot;

% Get cp_turbine
k_ratio = 0.3;
g_ratio = 1 - k_ratio;
cp_t = cp_tAIR;%(cp_tAIR + (g_ratio * cp_gas + k_ratio * cp_kerosene)*f_ratio)/(1 + f_ratio);
Rt = 54.09/ 778.16; % BTU/lbm
gamma_t = 1 / (1 - Rt/cp_t);
fprintf('Turbine cp: %0.5f\n', cp_t);

% Get Tt3
hPR_fuel = (k_ratio * hPR_kerosene + g_ratio*hPR_gas);
Tt3 = (cp_t / cp_c) * Tt4 - (f_ratio/cp_c)*(eta_b*hPR_fuel - cp_t * Tt4);
fprintf('Tt3 = %4.5f\n', Tt3);
tau_c = Tt3 / Tt0;
pi_c = tau_c ^ (gamma_c * e_c / (gamma_c - 1));
fprintf('pi_c: %0.5f\n', pi_c/FPR); 

% Get Tau_T from pi_c and the pi_f
tau_f = FPR ^ ((gamma_c - 1)/(gamma_c * e_c));
tau_t = 1 - (cp_c * Tt0)/(cp_t * Tt4) * 1./(eta_m * (1 + f_ratio))*(tau_c -1 + BPR * (tau_f - 1));
fprintf('tau_t: %0.5f\n', tau_t);
pi_t = tau_t ^ (gamma_t/((gamma_t - 1) * e_t));
fprintf('pi_t: %0.5f\n', pi_t);

% Get exit conditions
pi_t19_s = FPR * pi_n;
M19 = sqrt((2/(gamma_c - 1)) * (pi_t19_s ^ ((gamma_c - 1)/gamma_c) - 1));
tau_19_0 = tau_f / (pi_t19_s ^ ((gamma_c - 1)/gamma_c));

pi_t9_s = pi_c * pi_b * pi_t * pi_n;
M9 = sqrt((2/(gamma_t - 1)) * (pi_t9_s ^ ((gamma_t - 1)/gamma_t) - 1));

% To use with EngineSimU 1.8a


%% Get New Engine Power Curves?
fprintf('Thrust Specific Fuel Consumtion = %f \n',SFC_T*3600)

% FC.P0 = 2116.23 * atm.delta; % lbf/ft^2
% FC.T0 = 518.67 * atm.theta; % Rankine
% FC.M0 = 1.6;%req.cr_M0(1);
% FC.gamma = 1.4;  % air
% FC.cp = 0.24; % BTU/lbm R
% % FC.rho = 0.000531556; % slugs/ft^3
% 
% ENG.pi_c = pi_c;
% ENG.pi_f = FPR;
% ENG.tau_t = tau_t;
% ENG.tau_f = tau_f;
% ENG.Tt4_sl = Tt4;
% ENG.Tt0 = Tt0;
% ENG.tau_lam = cp_t * Tt4 / (cp_c * Tt0);
% ENG.tau_r = 1;
% ENG.Tt4_max = 4000; % TECH LEVEL 5
% ENG.Tt4_min = Tt3;
% ENG.hPR = hPR_kerosene; % replace Jet Fuel B with Jet Fuel A
% ENG.p0_9 = 1.0;
% ENG.p0_19 = 1.0;
% ENG.d = 54/12; % in -> ft
% ENG.BPR = BPR;
% ENG.mdot_sl = mdot;
% ENG.gam_t = gamma_t;
% ENG.M19 = M19;
% ENG.M9 = M9;
% ENG.Rt = Rt * 778.16;
% ENG.P0 = 2116.23;
% 
% n_eng = 3; % three engines
% thrust = (Wt.WTO * Wt.fuel.w2_1 * Wt.fuel.w3_2) / (LD_Max * n_eng);
% 
% % TL level 5
% TL.d_max = 0.97; % supersonic aircraft
% TL.e_c = 0.91;
% TL.e_f = 0.92;
% TL.pi_b = 0.96;
% TL.eta_b = 0.999;
% TL.e_t = 0.91; % uncooled
% TL.pi_n = 0.99; % variable area convergent/divergent nozzle
% TL.eta_m = 0.996; % shaft only
% 
% % % % TL Level 2
% % TL.eta_b = 0.94;
% % TL.e_c = 0.84;
% % TL.e_f = 0.82;
% % TL.e_t = 0.89; % uncooled
% % TL.pi_b = 0.92;
% % TL.pi_n = 0.97;
% % TL.d_max = 0.95; % type C
% % TL.eta_m = 0.97;  % shaft only
% 
% % Ma = linspace(0, 1.6, 100);
% % Tts4 = zeros(size(Ma));
% % f_ratio_T = zeros(size(Ma));
% % SFC_T = f_ratio_T;
% % for i = 1:length(Ma)
% %     FC.M0 = Ma(i);
% [Tts4, f_ratio_T, SFC_T] = throttle_calcs(FC, ENG, TL, thrust);
% % end
% 
% figure();
% plot(Ma, SFC_T);
% figure();plot(Ma, Tts4);