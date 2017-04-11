%% Cost Analysis from Raymer CH 18 pg 501-517
% Using Variable Notation from Raymer
% Needed or new parameters include: 
% Engine max mach number, Turbine Inlet Temp, Avionics cost, Max Thrust
% How many cycles or flights will be conducted per year?
%% Inital Parameters
clc
clear 
close all
% Will be replaced by master variable file later on
% (*) symbol indicates things that need to be changed
WTO = 80000; % lbf*
W_E = 40000; % lbf*
V = 1.6*573; % maximum velocity in knots*
Q = 1; % Production Quantity
FTA = 2; % Number of flight test aircraft (typically 2-6)
N_e = 3; % Number of Engines
C_e = 1000000; % engine cost
N_Eng = (Q + FTA)*N_e; % total production quantity times number of engines per aircraft
M_max = 2.3; % engine maximum mach number
T_max = 4384; % engine max thrust (lbf)
T_4 = 2000; % turbine inlet temperature (Rankine)
C_avionics = 1000000; % avionics cost

% Hourly Rates in 1986 constant dollars
R_E = 59.10; % engineering
R_T = 60.70; % tooling
R_Q = 55.40; % quality control
R_M = 50.10; % manufacturing
%% Modified DAPCA IV Cost Model pg 508-517
%% RDTE and Production Costs
% Cost in constant 1986 dollars not including inflation

% Fudge Factors for different materials for hour calculations
% (aluminum: 1.0, graphite-epoxy = 1.5-2.0, fiberglass 1.1-1.2, steel
% 1.5-2.0, titanium 1.7-2.2)

FF = 1; % currently for aluminum aircraft

% Engineering Hours 
H_E = FF*(4.86*W_E^0.777*V^0.894*Q^0.163);

% Tooling Hours Eq 18.2
H_T = FF*(5.99*W_E^0.777*V^0.696*Q^0.263);

% Mfg (manufacturing) Hours Eq 18.3
H_M = FF*(7.37*W_E^0.82*V^0.484*Q^0.641);

% QC (quality control) Hours Eq 18.4
H_Q = 0.133*H_M; % (mfg hours)for other than cargo plane

% Development support costs Eq 18.5
C_Dev = 45.42*W_E^0.630*V^1.3;

% Flight Test Cost Eq 18.6
C_F = 1243.03*W_E^0.325*V^0.822*FTA^1.21;

% Mfg Materials Cost Eq 18.7
C_Mat = 11*W_E^0.921*V^0.621*Q^0.799;

% Eng Production Cost Eq 18.8
C_Eng = 1548*(0.043*T_max+ 243.25*M_max + 0.969*T_4 - 2228);

% RDT&E (research development test evaluation) + flyaway COST Eq 18.9 
RDTE = H_E*R_E + H_T*R_T + H_Q*R_Q + H_M*R_M + C_Dev + C_F + C_Mat + ...
    C_Eng*N_Eng + C_avionics;
fprintf('The total aircraft cost is %0.1f \n',RDTE)
%% Operations and Maintanence Costs
% Two man crew cost Eq 18.10
C_crew = 35*(V*WTO/10^5)^0.3 + 84;

% Table 8.1 pg 511
MMH_FH = 4.5; % Maintenance man hours per flight hour business jet (avg)
FH_YR = 1250; % Flight hour per year per aircraft


% Calculating aircraft cost one less engine

N_e = 2; % Number of Engines
N_Eng = (Q + FTA)*N_e;
% Engineering Hours 
H_E = FF*(4.86*W_E^0.777*V^0.894*Q^0.163);
% Tooling Hours Eq 18.2
H_T = FF*(5.99*W_E^0.777*V^0.696*Q^0.263);
% Mfg (manufacturing) Hours Eq 18.3
H_M = FF*(7.37*W_E^0.82*V^0.484*Q^0.641);
% QC (quality control) Hours Eq 18.4
H_Q = 0.133*H_M; % (mfg hours)for other than cargo plane
% Development support costs Eq 18.5
C_Dev = 45.42*W_E^0.630*V^1.3;
% Flight Test Cost Eq 18.6
C_F = 1243.03*W_E^0.325*V^0.822*FTA^1.21;
% Mfg Materials Cost Eq 18.7
C_Mat = 11*W_E^0.921*V^0.621*Q^0.799;
% Eng Production Cost Eq 18.8
C_Eng = 1548*(0.043*T_max+ 243.25*M_max + 0.969*T_4 - 2228);
% RDT&E (research development test evaluation) + flyaway COST Eq 18.9 
C_a = H_E*R_E + H_T*R_T + H_Q*R_Q + H_M*R_M + C_Dev + C_F + C_Mat + ...
    C_Eng*N_Eng + C_avionics;
fprintf('The aircraft cost minus one engine is %0.1f \n',C_a)

% material cost per flight hour Eq 18.12
MC_FH = 3.3*(C_a/10^6) + 7.04 + (58*(C_e/10^6) - 13)*N_e;
fprintf('The material cost per flight hour is %0.1f \n',MC_FH)

% material cost per cycle Eq 18.13
MC_cycle = 4.0*(C_a/10^6) + 4.6 + (7.5*(C_e/10^6) + 2.8)*N_e;
fprintf('The material cost per cycle is %0.1f \n',MC_cycle)

% "cycles" are the number of flights
% The total materials cost is the cost per flight hour times the 
% flight hours per year, plus the cost per cycle times the cycles per year.
