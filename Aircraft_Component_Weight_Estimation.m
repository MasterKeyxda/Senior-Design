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
% dfkjdf
% Units: British / Imperial
% Note: All equations and tables referenced are from Sadraey unless
% otherwise stated.

clc;
clear;

%% ------------------------------ Wing --------------------------------- %%

% Section: 10.4.1
% Eqn(s) 10.3, 10.4

% Equation Parameters
Sw = 890; % wing planform area (ft^2)
MAC = 6; % mean aerodynamic chord (ft)
thickToChordRatio = 0.12; % max thickness to chord ratio
densityWing = 0.1015; % density of wing material (lb/in^3)
AR = 3; % aspect ratio
sweepQuarter = 15; % 1/4 chord sweep angle (degrees)
taperRatio = 2; % taper ratio
g = 32.17; % gravitational constant (32.17 ft/s^2)
KpWing = 0.0025; % wing density factor (Table 10.8)
nMax = 3; % max positive load factor (Table 10.9, Transport Aircraft) 
nUlt = 1.5 * nMax; % ultimate load factor

% Wing Weight
WtWing = Sw * MAC * thickToChordRatio * densityWing * KpWing * (((AR * nUlt) / cosd(sweepQuarter))^0.6) * (taperRatio^0.04) * g; 

%% ------------------------   Horizontal Tail  ------------------------- %%

%% --------------------------- Vertical Tail --------------------------- %%

%% ----------------------------- Fuselage ------------------------------ %%

%% ---------------------------- Landing Gear --------------------------- %%

%% ---------------------- Installed Engine Weight ---------------------- %%

%% ----------------------- Equipment + Subsystems  --------------------- %%

%% ------------------------  Total Weight ------------------------------ %%


%% ------------------------ MATLAB --> Excel --------------------------- %%