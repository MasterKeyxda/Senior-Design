function [OUT] = reverse_engineer_engine(specs, fuel, TL, air)


% GIVENS - JT8D-7B
% specs.mdot = 248; % lb/s
% specs.SFC = 0.74;
% specs.Tt4 = 2700; % R
% specs.Tt0 = 518.69; % R
% specs.Force = 17.8e3; %lbf
% specs.BPR = 0.36; % bypass ratio
% specs.FPR = 3.8; % Fan Pressure Ratio
% specs.OPR = 16.0;

% TL level 2
% TL.eta_b = 0.94;
% TL.e_c = 0.84;
% TL.e_f = 0.82;
% TL.e_t = 0.89; % uncooled
% TL.pi_b = 0.92;
% TL.pi_n = 0.97;
% TL.pi_d = 0.95; % type C
% TL.eta_m = 0.97;  % shaft only

f_ratio = specs.SFC/3600 * (1+specs.BPR) * specs.Force / specs.mdot;

% Get cp_turbine
cp_t = fuel.cp_tAIR;%(cp_tAIR + (g_ratio * cp_gas + k_ratio * cp_kerosene)*f_ratio)/(1 + f_ratio);
Rt = 54.09/ 778.16; % BTU/lbm
gamma_t = 1 / (1 - Rt/cp_t);
fprintf('Turbine cp: %0.5f\n', cp_t);

% Get Tt3
% hPR_fuel = (k_ratio * hPR_kerosene + g_ratio*hPR_gas);
Tt3 = (cp_t / air.cp_c) * specs.Tt4 - (f_ratio/air.cp_c)*(TL.eta_b*fuel.hPR - cp_t * specs.Tt4);
fprintf('Tt3 = %4.2f R\n', Tt3);
tau_c = Tt3 / specs.Tt0;
pi_c = tau_c ^ (air.gamma_c * TL.e_c / (air.gamma_c - 1));
fprintf('pi_c: %0.5f\n', pi_c/specs.FPR); 

% Get Tau_T from pi_c and the pi_f
tau_f = specs.FPR ^ ((air.gamma_c - 1)/(air.gamma_c * TL.e_c));
tau_t = 1 - (air.cp_c * specs.Tt0)/(cp_t * specs.Tt4) * 1./(TL.eta_m * (1 + f_ratio))*(tau_c -1 + specs.BPR * (tau_f - 1));
fprintf('tau_t: %0.5f\n', tau_t);
pi_t = tau_t ^ (gamma_t/((gamma_t - 1) * TL.e_t));
fprintf('pi_t: %0.5f\n', pi_t);

% Output
OUT.pi_t = pi_t;
OUT.tau_t = tau_t;
OUT.pi_c = pi_c;

end