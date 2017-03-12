clc;
clear;
close all;

%% INFO:

% want to iterate weights to get new weight
% set weight tolerance
ctrl.tol = 0.01; % 5 percent

%% First Run of SAD.m

Supersonic_Aircraft_Design;
ctrl.Wt_old(1) = Wt.WTO; % Save old weight

% Get corrected weight
Class_2_Weights;

% Calculate Residual
ctrl.res(1) = WeightDiff / 100; %

clearvars -except Wt atm req ctrl
iterate = 1;

% Recalculate fuel fraction for new engines
% Cruise Fuel Consumption
Wt.fuel.sfc_cr = 0.6/3600; %0.572/3600; % 1/hr -> 1/s cruise, Table 4.6 (Sadrey)
[atm.delta,atm.theta,atm.sig_rho,Wt.fuel.a_snd] = AltTable(atm.alt,'h'); % speed of sound ratio
Wt.fuel.V_max_cr = Wt.fuel.a_snd * req.cr_M0(1) * 1116; % ft/s
Wt.fuel.LD_ratio = 7; % based on past, real aircraft (e.g. concorde)
Wt.fuel.w4_3 = exp(-req.range * 6076.12 * Wt.fuel.sfc_cr/(0.866*Wt.fuel.V_max_cr * Wt.fuel.LD_ratio));

% Loiter Fuel Consumption
Wt.fuel.t_loiter = 1; % hr of loiter
Wt.fuel.sfc_loiter = 0.7; % 1/hr
Wt.fuel.w5_4 = exp(-Wt.fuel.t_loiter*Wt.fuel.sfc_loiter / Wt.fuel.LD_ratio);

Wt.fuel.w6_5 = 0.99; % descent
Wt.fuel.w7_6 = 0.997; % approach & landing
Wt.fuel.w7_1 = Wt.fuel.w2_1 * Wt.fuel.w3_2 * Wt.fuel.w4_3 * Wt.fuel.w5_4 * Wt.fuel.w6_5 * Wt.fuel.w7_6;

Wt.fuel.Wf_Wto = 1 - Wt.fuel.w7_1; % Fuel ratio
fprintf('Fuel Ratio (Pre-Reserve): %0.5f\n', Wt.fuel.Wf_Wto);

Wt.fuel.reserve_ratio = 1.25; % Reserve fuel at least 20% FAR from Sadraey pg. 102 
Wt.fuel.w_max = (1/req.f_eff)*(Wt.pld.n_pass * req.range)*Wt.fuel.reserve_ratio;
% Wt.fuel.w_tot = Wt.fuel.w_max;
fprintf('Maximum Fuel Weight Allowed: %0.2f lb\n', Wt.fuel.w_max);

Wt.enginetype.name = 'PW-TF33-3-7';
Wt.enginetype.w_tot = 4650; 
Wt.enginetype.thr = 21e3;

% WF_TO = .47948; % taken from main code loiter 0.75%
% Trainer Jet Raymer Table 3.1 pg 31 
A = 1.59;
C = -0.1;
% Solve for WTO with Eq 3.4 and empirical Table 3.1 Eq
y = linspace(80000,170000,1000);
x1 = A*y.^C; % Raymer table 3.1
x2 = 1- Wt.fuel.Wf_Wto - (Wt.pld.w_tot + Wt.oew.crew)./y; % Sadrey eqn. 4.5
% plot(y,x1,y,x2) 
dif = abs(x1 - x2);
index = find(dif == min(abs(x1 - x2)));
Wt.WE_WTO = x1(index);
Wt.WTO = y(index);

fprintf('WE_WTO: %0.3f \n',Wt.WE_WTO);
fprintf('W_TO: %.2f lbs \n',Wt.WTO);

% Clear out other variables except for Wt
% clearvars -except Wt atm req ctrl

% Empty weight (lbs)
% Wt.WE = Wt.WTO - Wt.fuel.w_tot - Wt.pld.w_tot - Wt.oew.crew;
Wt.fuel.w_tot = Wt.fuel.Wf_Wto * Wt.WTO;
Wt.WE = Wt.WE_WTO * Wt.WTO; 

fprintf('WE_WTO: %0.3f \n', Wt.WE_WTO);
fprintf('W_TO: %.2f lbs \n', Wt.WTO);
fprintf('WE: %0.2f lbs \n', Wt.WE); 

while ctrl.res(end) > ctrl.tol
   iterate = iterate + 1;
    %% Constraints Plots
    fprintf('\n CONSTRAINT PLOTS: \n');
    Constraint_Plots;

    %% Wing Calculations

    fprintf('\n WING CALCULATIONS: \n');
    WINGS;

    %% Preliminary Fuselage Design
    fprintf('\n FUSELAGE CALCULATIONS: \n');
    Preliminary_Fuselage_Design;

    %% Tail Calculations

    fprintf('\n TAIL CALCULATIONS: \n');
    % Sample Tail Params
    Kc = 1.2; % tail calculation correction factor?
    Vh = 0.6; % horizontal tail volumen coefficient
    Vv = 0.05; % vertical tail volumen coefficient
    sweepWing = 28; % wing sweep degrees
    taperh = 0.6; % horizontal tail taper ratio
    cglocAC = -5.5; % ft cg location in front or behind AC Wing
    TAIL = TailCalc(0, Vh, Vv, Wt.WTO, atm.sig_rho * atm.rho_sl, Wt.fuel.V_max_cr, D_C, Kc, WING.geom.S_area, WING.geom.AR, WING.Cmwf, sweepWing, taperh, cglocAC, '', req.cr_M0(1));

    %% V-n diagram

    fprintf('\n V-N DIAGRAM CALCULATIONS: \n');
    V_n_diagram;
    ctrl.Wt_old(iterate) = Wt.WTO; % Save old weight
    
    % Get corrected weight
    Class_2_Weights;
    
    % Calculate Residual
    ctrl.res(iterate) = WeightDiff / 100; %
    
    if ctrl.res(iterate) > ctrl.tol
        clearvars -except Wt atm req ctrl iterate
        close all;
    end
    
end

fprintf('No. of Iterations: %i\n', iterate);
fprintf('Final Residual: %0.5f\n', ctrl.res(end));

if Wt.enginetype.thr*3 < constraints.req_Thr
    fprintf('Thrust requirement not met!!!\n');
else
    fprintf('Thrust requirement met\n');
end

% save meta information about current script run
meta.date = datetime('today');

% save variables to .mat file
clear ctrl;
save('aircraft_vars.mat'); 