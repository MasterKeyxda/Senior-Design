%% Aircraft Component Weight Calculation

% Date Created: 2/14/17
% Summary: This script calculates the weight of individual components of
% the aircraft using the methodology outlined in Chapter 10 of Sadraey. The
% estimates are based on four sources:
    % 1. direct relationship between weight of an object and its average density (Table 10.6);
    % 2. actual published data on weight of various components (e.g., Table 10.5);
    % 3. derived empirical factors by the author;
    % 4. published empirical equations
% The results of this script are written to an excel file for ease of table
% creation. 

% Units: English / Imperial
% Note: All equations and tables referenced are from Sadraey unless
% otherwise stated.

clc;
clear;

%% ------------------------------ Wing --------------------------------- %%

% Section: 10.4.1
% Eqn(s) 10.3, 10.4

% Equation Parameters
Sw = 890; % wing planform area (ft^2)
MACw = 6; % mean aerodynamic chord (ft)
thickToChordRatioWing = 0.12; % max thickness to chord ratio
densityWing = 0.1015; % density of wing material (lb/in^3)
AR_Wing = 3; % aspect ratio
sweepQuarterWing = 15; % 1/4 chord sweep angle (degrees)
taperRatioWing = 2; % taper ratio
g = 32.17; % gravitational constant (32.17 ft/s^2)
KpWing = 0.0025; % wing density factor (Table 10.8)
nMax = 3; % max positive load factor (Table 10.9)
nUlt = 1.5 * nMax; % ultimate load factor

% Wing Weight
WtWing = Sw * MACw * thickToChordRatioWing * densityWing * KpWing * (((AR_Wing * nUlt) / cosd(sweepQuarterWing))^0.6) * (taperRatioWing^0.04) * g; 

%% ------------------------   Horizontal Tail  ------------------------- %%

% Section: 10.4.2
% Eqn(s) 10.5

% Equation Parameters
Sht = 200; % exposed planform area
MACht = 2; % mean aerodynamic chord (ft)
thickToChordRatioHT = 0.07; % max thickness to chord ratio
densityHT = 0.1015; % density of horizontal tail material (lb/in^3)
AR_HT = 2; % aspect ratio
sweepQuarterHT = 15; % 1/4 chord sweep angle (degrees)
taperRatioHT = 2; % taper ratio
g = 32.17; % gravitational constant (32.17 ft/s^2)
chordRatio = 2; % elevator-to-tail chord ratio
volumeRatioHT = 30; % horizontal tail volume ratio
KpHT = 0.02; % horizontal tail density factor (Table 10.10)

% Horizontal Tail Weight
WtHorizontalTail = Sht * MACht * thickToChordRatioHT * densityHT * KpHT * ((AR_HT / (cosd(sweepQuarterHT)))^0.6) * (taperRatioHT^0.04) * (volumeRatioHT^0.3) * (chordRatio^0.4) * g;

%% --------------------------- Vertical Tail --------------------------- %%

% Section: 10.4.3
% Eqn(s) 10.6

% Equation Parameters
Svt = 200; % exposed planform area
MACvt = 2; % mean aerodynamic chord (ft)
thickToChordRatioVT = 0.07; % max thickness to chord ratio
densityVT = 0.1015; % density of vertical tail material (lb/in^3)
AR_VT = 2; % aspect ratio
sweepQuarterVT = 15; % 1/4 chord sweep angle (degrees)
taperRatioVT = 2; % taper ratio
volumeRatioVT = 30; % vertical tail volume ratio
chordRatio = 2; % rudder-to-vertical tail chord ratio
KpVT = 0.035; % vertical tail density factor (Table 10.10)

% Vertical Tail Weight
WtVerticalTail = Svt * MACvt * thickToChordRatioVT * densityVT * KpVT * ((AR_VT / (cosd(sweepQuarterVT)))^0.6) * (taperRatioVT^0.04) * (volumeRatioVT^0.2) * (chordRatio^0.4) * g;

%% ----------------------------- Fuselage ------------------------------ %%

% Section: 10.4.4
% Eqn(s) 10.7

% Equation Parameters
Lf = 32; % fuselage length (ft)
DfMax = 8; % max diameter of equivalent circular cross-section (ft)
densityFuse = 0.1015; % density of fuselage material (lb/in^3)
KpFuse = 0.0032; % fuselage density factor (Table 10.11)
KInlet = 1.25; % 1.25 for inlets on fuselage or 1 for inlets elsewhere

% Fuselage Weight
WtFuselage = Lf * (DfMax^2) * densityFuse * KpFuse * (nUlt^0.25) * KInlet * g; 

%% ---------------------------- Landing Gear --------------------------- %%

% Section: 10.4.5
% Eqn(s) 10.8

% Equation Parameters
KL = 1; % landing place factor (1.8 for Navy aircraft and 1 otherwise)
Kret = 1.07; % 1 for fixed landing gear and 1.07 for retractable landing gear
KLandingGear = 0.28; % landing gear weight factor (Table 10.12)
WLanding = 70000; % aircraft weight at landing (lbf)
LGHeight = 10; % landing gear height (ft)
b = 20; % wing span (ft)
nMaxLG = 3; % max load factor for landing gear (Table 10.9)
nUltLG = 1.5 * nMaxLG; % ultimate load factor

% Landing Gear Weight
WtLandingGear = KL * Kret * KLandingGear * WLanding * (LGHeight / b) * (nUltLG^0.2);

%% -------------------------- Installed Engine  ------------------------ %%

% Section: 10.4.6
% Eqn(s) 10.9

% Equation Parameters
KE = 2.6; % engine weight factor (2.6 for British units, 3 for metric units)
NE = 2; % number of engines
WtEngine = 1; % weight of single engine (lbf)

% Installed Engine Weight
WtInstalledEngine = KE * NE * (WtEngine^0.9);

%% ----------------------------- Fuel System --------------------------- %%

% Section: 10.4.7
% Eqn(s): 10.10 or 10.11 or 10.12

% Equation Parameters

% Fuel System Weight
WtFuelSystem = 40000; 

%% ----------------------- Equipment + Subsystems  --------------------- %%

% Section: 10.4.8

% Recommended to be assessed using a variety of sources (see below).
% Sadraey - Table 10.13 gives weights of seat, lavatories, air instruments
% parachutes, stick, yoke, and wheels
% Roskam - Airplane Design, Part V
% Torenbeek, Synthesis of Subsonic Airplane Design

% Equipment and Subsystem Weight
WtEquipment = 5100; 

%% --------------------------  Total Weight ---------------------------- %%

% Finds total weight
WtTotal = WtWing + WtHorizontalTail + WtVerticalTail + WtFuselage + ...
    WtLandingGear + WtInstalledEngine + WtFuelSystem + WtEquipment; 

% Weight Breakdown
fprintf('Wing Weight = %0.2f lbf (%0.2f percent of total weight) \n', WtWing, (WtWing / WtTotal) * 100)
fprintf('Horizontal Weight = %0.2f lbf (%0.2f percent of total weight) \n', WtHorizontalTail, (WtHorizontalTail / WtTotal) * 100)
fprintf('Vertical Weight = %0.2f lbf (%0.2f percent of total weight) \n', WtVerticalTail, (WtVerticalTail / WtTotal) * 100)
fprintf('Fuselage Weight = %0.2f lbf (%0.2f percent of total weight) \n', WtFuselage, (WtFuselage / WtTotal) * 100)
fprintf('Landing Gear Weight = %0.2f lbf (%0.2f percent of total weight) \n', WtLandingGear, (WtLandingGear / WtTotal) * 100)
fprintf('Installed Engine Weight = %0.2f lbf (%0.2f percent of total weight) \n', WtInstalledEngine, (WtInstalledEngine / WtTotal) * 100)
fprintf('Fuel System Weight = %0.2f lbf (%0.2f percent of total weight) \n', WtFuelSystem, (WtFuelSystem / WtTotal) * 100)
fprintf('Equipment + Subsystem Weight = %0.2f lbf (%0.2f percent of total weight) \n', WtEquipment, (WtEquipment / WtTotal) * 100)

%% ------------------------ MATLAB --> Excel --------------------------- %%

% Populate excel spreadsheet with main component weights
%filename = 'Aircraft_Component_Weight.xlsx';
%sheet = 1; 
%xlRange = 'D1';
%xlswrite(filename,WtWing,sheet,xlRange)
