clc;
clear;
close all;

% load('aircraft_vars.mat');

%% Reverse calculate the pi_c and pi_t

engine_list;

% Air Properties
air.cp_c = 0.24; %BTU/lbm R
air.gamma_c = 1.4;

% Kerosene Properties
k_ratio = 0;
cp_kerosene = 0.48;
hPR_kerosene = 18300; % BTU/lbm

% Gasoline Properties
g_ratio = 1;
cp_gas = 0.53;
hPR_gas = 18443.6801; % BTU/lbm

% Fuel Properties
fuel.cp_tAIR = 0.295;
fuel.hPR = (k_ratio * hPR_kerosene + g_ratio*hPR_gas);

% TL level 2
TL2.eta_b = 0.94;
TL2.e_c = 0.84;
TL2.e_f = 0.82;
TL2.e_t = 0.89; % uncooled
TL2.pi_b = 0.92;
TL2.pi_n = 0.97;
TL2.pi_d = 0.95; % type C
TL2.eta_m = 0.97;  % shaft only

% TL level 5
TL.d_max = 0.97; % supersonic aircraft
TL.e_c = 0.91;
TL.e_f = 0.92;
TL.pi_b = 0.96;
TL.eta_b = 0.999;
TL.e_t = 0.91; % uncooled
TL.pi_n = 0.99; % variable area convergent/divergent nozzle
TL.eta_m = 0.996; % shaft only

TF33_parts = reverse_engineer_engine(TF33, fuel, TL2, air);
F100_parts = reverse_engineer_engine(F100_PW_229, fuel, TL2, air);


%% Engine Throttle Params & Get New Engine Power Curves?
% fprintf('Thrust Specific Fuel Consumtion = %f \n',SFC_T*3600)

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