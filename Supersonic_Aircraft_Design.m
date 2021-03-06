%% SETUP
% clear;
% if ~exist('ctrl', 'var')
%    clc;
%    clearvars -except Wt atm req ctrl
%    close all;
% else
%     clear;
%     clc;
%     close all;
% end

clc;
clearvars -except Wt atm req ctrl pass
close all;

%% File META - Preliminary Design Elements

% Requirements: 
    % Sonic Boom: 70-75 PLdB
    % Payload: 6-20 Passengers
    % Airport Noise: ICAO Ch.14 w/ margin
    % Equivalent emissions to current subsonic aircraft
    % Cruise Speed: Mach 1.6-1.8 
    % Range 4000 nmi
    % Fuel_Efficiency = 1; % passenger-miles/lb of fuel

% Atmosphere/Flight Conditions
    % Altitude: 40 kft
    
% Weight
    % number of passengers: 8
    % number of crew: 2

% Wings
    % Initial Wing Area - 930 ft^2

%% Supersonic Business Jet

atm.alt = 42; %kft
atm.rho_sl = 0.0023769; % slugs/ft^3

% N+1 REquirements
req.sb_vol = 75; % PLdB
req.pld_pass = [6 20]; % [min max] number of passengers
req.cr_M0 = [1.6 1.8]; % [min max] vel
req.range = 4000; % nmi, min
req.f_eff = 1.0; % passenger-miles/lb of fuel (min)

% FAR 25 Requirements
req.takeoffRun = 7000; % NASA takeoff field length of less than 7000 ft; based
% on FAR 25, aircraft must clear imaginary 35 ft obstacle

%% Payload Weight (pld)

% Total Payload Weight = Weight of passengers + weight of luggage + weight of fuel + weight crew*
% Total Weight = Total Payload Weight + Total Aircraft weight

fprintf('\t\t PRELIMINARY WEIGHT ESTIMATES: \n');

if ~exist('pass', 'var')
    Wt.pld.n_pass = 10; % Number of Passengers
    fprintf('Passenger study not on!\n');
else
    Wt.pld.n_pass = pass.num;
    warning('Pasenger study on!\n');
end
Wt.pld.apw = 170; % lbf, average passenger weight (apw)
Wt.pld.lug = 45; % lbf, luggage
Wt.pld.w_tot = Wt.pld.n_pass * (Wt.pld.apw + Wt.pld.lug); % lbf

%% Fuel Weight (fuel)

Wt.fuel.w2_1 = 0.98; % taxi and takeoff weight ratio - Table 4.3 (Sadrey)
Wt.fuel.w3_2 = 0.97; % climb weight ratio

% Cruise Fuel Consumption
Wt.fuel.sfc_cr = 0.8/3600; % 1/hr -> 1/s cruise, Table 4.6 (Sadrey)
[atm.delta,atm.theta,atm.sig_rho,Wt.fuel.a_snd] = AltTable(atm.alt,'h'); % speed of sound ratio
Wt.fuel.V_max_cr = Wt.fuel.a_snd * req.cr_M0(1) * 1116; % ft/s
Wt.fuel.LD_ratio = 7; % based on past, real aircraft (e.g. concorde)
Wt.fuel.w4_3 = exp(-req.range * 6076.12 * Wt.fuel.sfc_cr/(0.866*Wt.fuel.V_max_cr * Wt.fuel.LD_ratio));

% Loiter Fuel Consumption
Wt.fuel.t_loiter = 1; % hr of loiter
Wt.fuel.sfc_loiter = 0.8; % 1/hr
Wt.fuel.w5_4 = exp(-Wt.fuel.t_loiter*Wt.fuel.sfc_loiter / Wt.fuel.LD_ratio);

Wt.fuel.w6_5 = 0.99; % descent
Wt.fuel.w7_6 = 0.997; % approach & landing
Wt.fuel.w7_1 = Wt.fuel.w2_1 * Wt.fuel.w3_2 * Wt.fuel.w4_3 * Wt.fuel.w5_4 * Wt.fuel.w6_5 * Wt.fuel.w7_6;

Wt.fuel.Wf_Wto = 1 - Wt.fuel.w7_1; % Fuel ratio
fprintf('Fuel Ratio (Pre-Reserve): %0.5f\n', Wt.fuel.Wf_Wto);

Wt.fuel.reserve_ratio = 1.25; % Reserve fuel at least 20% FAR from Sadraey pg. 102 
Wt.fuel.w_max = (1/req.f_eff)*(Wt.pld.n_pass * req.range)*Wt.fuel.reserve_ratio;
% Wt.fuel.w_tot = Wt.fuel.w_max;
fprintf('Maximum Fuel Weight Allowed to be burned: %0.2f lb\n', Wt.fuel.w_max);

Wt.enginetype.name = 'Pegasus';
Wt.enginetype.w_tot = 3960; %lbm

%% Operating Empty Weight (oew)

% W_OE = W_E + W_TFO + W_CREW
% W_E = W_ME + W_FEQ

% Solve for WTO with Empirical Model 4.26 Saedray

% Crew Weight
Wt.oew.n_crew = 2; % Number of Crew Members
Wt.oew.baggage = 30; % crew has 30 lbs baggage
Wt.oew.crew = Wt.oew.n_crew*(Wt.pld.apw+Wt.oew.baggage); %lb

% Manufacturer's empty weight - empty airframe

% Fixed equipment weight
Wt.oew.seat = 33; % lbf  Aircraft Design - Sadre
Wt.oew.feq = Wt.oew.seat*(Wt.oew.n_crew + Wt.pld.n_pass);

% Trapped Fuel and Oil
Wt.oew.tfo = 0;

% Total Operating Empty Weight
Wt.oew.w_tot = Wt.oew.tfo + Wt.oew.crew + Wt.oew.feq;


%% Weight Display
% 
% wt_types = fieldnames(Wt);
% for i = 1:numel(wt_types)
%     fprintf('%s weight: %0.2f lb\n', wt_types{i}, Wt.(wt_types{i}).w_tot);
% end
% 
% clear wt_types;

%% WTO and WE/WTO Calculation

% WF_TO = .47948; % taken from main code loiter 0.75%
% Trainer Jet Raymer Table 3.1 pg 31 
A = 1.59;
C = -0.1;
% Solve for WTO with Eq 3.4 and empirical Table 3.1 Eq
y = linspace(80000,170000,1000);
x1 = A*y.^C;
x2 = 1- Wt.fuel.Wf_Wto - 1960./y;
plot(y,x1,y,x2) 
dif = abs(x1 - x2);
index = find(dif == min(abs(x1 - x2)));
Wt.WE_WTO = x1(index);
Wt.WTO = y(index);

fprintf('WE_WTO: %0.3f \n',Wt.WE_WTO);
fprintf('W_TO: %.2f lbs \n',Wt.WTO);

% Clear out other variables except for Wt
clearvars -except Wt atm req ctrl pass

% Empty weight (lbs)
% Wt.WE = Wt.WTO - Wt.fuel.w_tot - Wt.pld.w_tot - Wt.oew.crew;
Wt.fuel.w_tot = Wt.fuel.Wf_Wto * Wt.WTO;
Wt.WE = Wt.WE_WTO * Wt.WTO;

fprintf('WE_WTO: %0.3f \n', Wt.WE_WTO);
fprintf('W_TO: %.2f lbs \n', Wt.WTO);
fprintf('WE: %0.2f lbs \n', Wt.WE); 

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
cglocAC = -12; % ft cg location in front or behind AC Wing
TAIL = TailCalc(0, Vh, Vv, Wt.WTO, atm.sig_rho * atm.rho_sl, Wt.fuel.V_max_cr, D_C, Kc, WING.geom.S_area, WING.geom.AR, WING.Cmwf, sweepWing, taperh, cglocAC, '', req.cr_M0(1));

%% V-n diagram

fprintf('\n V-N DIAGRAM CALCULATIONS: \n');
V_n_diagram;

%% Get run-time meta info for future reference

% save meta information about current script run
% meta.date = datetime('today');

% save variables to .mat file
save('aircraft_vars.mat'); 