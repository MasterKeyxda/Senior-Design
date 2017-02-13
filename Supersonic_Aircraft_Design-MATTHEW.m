%% SETUP
clear;
clc;
close all;

%% File META 

% File that contains initial approximations for aircraft design. Includes
% approximations for weight of aircraft, and an overview of requirements.
% Below contains documented updates through each meeting.

% 2/8/2017:
    % Added Fuel Mass Ratio Calculations
    % Weight
    % Wing Loading (90 lb/ft^2)
    % Design Wing
        % pick AR
        % pick Cl for flight (between 0.1 - 0.3)
    % Design Tail
        % L-opt -> search in Sadrey
    % Verify if CG is between 30-35% MAC (aim for 35) -> prelim stability
    % done
    % Design fuselage shape
    % Check volume... make sure payload CG matches aircraft CG
    % Iterate to meet fuel requirments

% 2/6/2017:
    % Determine Empty Weight,Fuel Weight and MTOW (Max Takeoff Weight)
    % Changes made:
    % Changed number of crew members from 3 to 2.
    % Will assume that crew members have the same luggage weight
    % Number of passengers changed to 8
    % Current avg weight of passengers 200 lbf and 45 lbf luggage
    % Found out for airworthyness that excess fuel must be at least 20% for FAR
    % Searched for possible low bypass turbofan engines

% 2/2/2017:
    % File Created
    % Want to carry 20 people - burn 20 lbs/miles
    % Currently assume crew members not part of passengers, but part of vehicle
    % system
    % Keep in mind that geometry (specifically volume) drives aerodynamic performance
    % This is also true for fuel tank volume
    % Fuel consumption and fuel weight?
    % Redefine altitude?

%% Supersonic Business Jet

% Requirements: 
% Sonic Boom: 70-75 PLdB
% Payload: 6-20 Passengers
% Airport Noise: ICAO Ch.14 w/ margin
% Equivalent emissions to current subsonic aircraft
% Cruise Speed: Mach 1.6-1.8 
% Range 4000 nmi
% Fuel_Efficiency = 1; % passenger-miles/lb of fuel

req.sb_vol = 75; % PLdB
req.pld_pass = [6 20]; % [min max] number of passengers
req.cr_M0 = [1.6 1.8]; % [min max] vel
req.range = 4000; % nmi, min
req.f_eff = 1.0; % passenger-miles/lb of fuel (min)

%% Payload Weight (pld)

% Total Payload Weight = Weight of passengers + weight of luggage + weight of fuel + weight crew*
% Total Weight = Total Payload Weight + Total Aircraft weight

Wt.pld.n_pass = 8; % Number of Passengers
Wt.pld.apw = 200; % lbf, average passenger weight (apw)
Wt.pld.lug = 45; % lbf, luggage
Wt.pld.w_tot = Wt.pld.n_pass * (Wt.pld.apw + Wt.pld.lug); % lbf

%% Fuel Weight (fuel)

Wt.fuel.w2_1 = 0.98; % taxi and takeoff weight ratio - Table 4.3 (Sadrey)
Wt.fuel.w3_2 = 0.97; % climb weight ratio

% Cruise Fuel Consumption
Wt.fuel.sfc_cr = 0.7/3600; % 1/hr -> 1/s cruise, Table 4.6 (Sadrey)
[~,~,~,Wt.fuel.a_snd] = AltTable(50,'h'); % speed of sound ratio
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
Wt.fuel.w_tot = Wt.fuel.w_max;
fprintf('Maximum Fuel Weight Allowed: %0.2f lb\n', Wt.fuel.w_max);

%% Operating Empty Weight (oew)

% W_OE = W_E + W_TFO + W_CREW
% W_E = W_ME + W_FEQ

% Solve for WTO with Empirical Model 4.26 Saedray

% Solve for Total Weight using Fuel Weight Ratio
WTO = Wt.fuel.w_tot/Wt.fuel.Wf_Wto; % Total Takeoff Weight


% Crew Weight
Wt.oew.n_crew = 2; % Number of Crew Members
Wt.oew.crew = Wt.oew.n_crew*(Wt.pld.apw+Wt.pld.lug); %lbf

% Manufacturer's empty weight - empty airframe


% Fixed equipment weight
Wt.oew.seat = 33; % lbf  Aircraft Design - Sadre
Wt.oew.feq = Wt.oew.seat*(Wt.oew.n_crew + Wt.pld.n_pass);

% Trapped Fuel and Oil
Wt.oew.tfo = 0;

Wt.oew.w_tot = Wt.oew.tfo + Wt.oew.crew + Wt.oew.feq;


%% Weight Display

wt_types = fieldnames(Wt);
for i = 1:numel(wt_types)
    fprintf('%s weight: %0.2f lb\n', wt_types{i}, Wt.(wt_types{i}).w_tot);
end
fprintf('WTO: %0.2f lb\n',WTO)

%% Iterate to Get WTO
% Solve for Wt.oew.me with Empirical Model 4.26 Saedray
a = 1.13e-6; b = 0.48; % From table 4.8
WTO_1 = WTO;
WE_TO = a*WTO_1 + b;
WTO_2 = (Wt.pld.w_tot + Wt.oew.crew)/(1 - Wt.fuel.Wf_Wto - WE_TO);
res = (WTO_2 - WTO_1);
ii = 1;
fprintf('Iteration Started\n');
while abs(res) >= 1
    WTO_1 = WTO_2;
    WE_TO = a*WTO_1 + b;
    WTO_2 = (Wt.pld.w_tot + Wt.oew.crew)/(1 - Wt.fuel.Wf_Wto - WE_TO);
    res = (WTO_2 - WTO_1);

    ii = ii + 1;
    if ii > 1e6
       fprintf('Warning: Iteration Broken\n');
       break;
    end
end

fprintf('Iterated WTO: %0.2f lb \n', WTO_1);
fprintf('Number of Iterations: %i \n', ii);

