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

% General / Common Inputs
g = 32.17; % gravitational constant (32.17 ft/s^2)
nMax = 3; % max positive load factor (Table 10.9) (Tentative, need V-n Diagram)
nUlt = 1.5 * nMax; % ultimate load factor
dens.Al7075 = 0.102; % Aluminum 7075 density (lb/in^3); Source: MIL-HDBK-5J
dens.Al2090 = 0.094; % Aluminum 2090 density (lb/in^3); Source: MIL-HDBK-5J
dens.Ti_6Al_4V = 0.160; % Titanium 6Al-4V density (lb/in^3); Source: MIL-HDBK-5J
dens.AS4_Epoxy = 0.056; % CFRP AS4/Epoxy density (lb/in^3); % Source: Mechanics of 
% Fibrous Composites, 1997 by C. T. Herakovich	

% Convert material densities to lb/ft^3
dens.Al7075 = dens.Al7075 * 1728;
dens.Al2090 = dens.Al2090 * 1728;
dens.Ti_6Al_4V = dens.Ti_6Al_4V * 1728;
dens.AS4_Epoxy = dens.AS4_Epoxy * 1728;

%% ------------------------------ Wing --------------------------------- %%

% Section: 10.4.1
% Eqn(s) 10.3, 10.4

% General / Common Inputs

KpWing = 0.0025; % wing density factor (Table 10.8)
chordRatioWing = 0.05; % max thickness to chord ratio

% Wing Material Composition
% Assume wing is constructed from 50% composite and 50% aluminum
Wing.dens = (dens.AS4_Epoxy + dens.Al7075) / 2; 

% Subsonic Equation Parameters
Wing.Subsonic.Sw = 890; % wing planform area (ft^2)
Wing.Subsonic.MAC = 6; % mean aerodynamic chord (ft)
Wing.Subsonic.AR = 3; % aspect ratio
Wing.Subsonic.sweep = 15; % quarter chord sweep angle (degrees)
Wing.Subsonic.taperRatio = 2; % taper ratio

% Supersonic Equation Parameters
Wing.Supersonic.Sw = 890; % wing planform area (ft^2)
Wing.Supersonic.MAC = 6; % mean aerodynamic chord (ft)
Wing.Supersonic.AR = 3; % aspect ratio
Wing.Supersonic.sweep = 15; % quarter chord sweep angle (degrees)
Wing.Supersonic.taperRatio = 2; % taper ratio

% Split weight calculation into two sections: a subsonic and supersonic
% portion; Add individual weights together

% Subsonic Wing Weight
Wt.Wing.Subsonic = Wing.Subsonic.Sw * Wing.Subsonic.MAC * (chordRatioWing) * Wing.dens * KpWing * (((Wing.Subsonic.AR * nUlt) / (cosd(Wing.Subsonic.sweep)))^0.6) * (Wing.Subsonic.taperRatio^0.04) * g;

% Supersonic Wing Weight
Wt.Wing.Supersonic = Wing.Supersonic.Sw * Wing.Supersonic.MAC * (chordRatioWing) * Wing.dens * KpWing * (((Wing.Supersonic.AR * nUlt) / (cosd(Wing.Supersonic.sweep)))^0.6) * (Wing.Supersonic.taperRatio^0.04) * g;

% Total Wing Weight
Wt.Wing = Wt.Wing.Subsonic + Wt.Wing.Supersonic; 
%% ------------------------   Horizontal Tail  ------------------------- %%

% Section: 10.4.2
% Eqn(s) 10.5

% General / Common Inputs
KpHT = 0.022; % horizontal tail density factor (Table 10.10), Transport,T-Tail

% Horizontal Tail Material Composition
HorizontalTail.dens = dens.Al7075; % Aluminum 7075 density (lb/in^3)

% Horizontal Tail Equation Parameters
HorizontalTail.Area = 200; % exposed planform area
HorizontalTail.MAC = 2; % mean aerodynamic chord (ft)
HorizontalTail.chordRatio = 0.07; % max thickness to chord ratio
HorizontalTail.AR = 2; % aspect ratio
HorizontalTail.sweep = 15; % quarter chord sweep angle (degrees)
HorizontalTail.taperRatio = 2; % taper ratio
HorizontalTail.ElevToTail = 2; % elevator-to-tail chord ratio
HorizontalTail.volumeRatio = 30; % horizontal tail volume ratio (see Ch 6 Sadraey)

% Horizontal Tail Weight
Wt.HorizontalTail = HorizontalTail.Area * HorizontalTail.MAC * HorizontalTail.chordRatio * HorizontalTail.dens * KpHT * ((HorizontalTail.AR / (cosd(HorizontalTail.sweep)))^0.6) * (HorizontalTail.taperRatio^0.04) * (HorizontalTail.volumeRatio^0.3) * (HorizontalTail.ElevToTail^0.4) * g;

%% --------------------------- Vertical Tail --------------------------- %%

% Section: 10.4.3
% Eqn(s) 10.6

% General / Common Inputs
KpVT = 0.04; % vertical tail density factor (Table 10.10), Transport,T-Tail

% Vertical Tail Material Composition
VerticalTail.dens = dens.Al7075; % Aluminum 7075 density (lb/in^3)

% Vertical Tail Equation Parameters
VerticalTail.Area = 200; % exposed planform area
VerticalTail.MAC = 2; % mean aerodynamic chord (ft)
VerticalTail.chordRatio = 0.07; % max thickness to chord ratio
VerticalTail.AR = 2; % aspect ratio
VerticalTail.sweep = 15; % quarter chord sweep angle (degrees)
VerticalTail.taperRatio = 2; % taper ratio
VerticalTail.volumeRatio = 30; % vertical tail volume ratio
VerticalTail.RudderToTail = 2; % rudder-to-vertical tail chord ratio

% Vertical Tail Weight
Wt.VerticalTail = VerticalTail.Area * VerticalTail.MAC * VerticalTail.chordRatio * VerticalTail.dens * KpVT * ((VerticalTail.AR  / (cosd(VerticalTail.sweep)))^0.6) * (VerticalTail.taperRatio^0.04) * (VerticalTail.volumeRatio^0.2) * (VerticalTail.RudderToTail^0.4) * g;

%% ----------------------------- Fuselage ------------------------------ %%

% Section: 10.4.4
% Eqn(s) 10.7

% General / Common Inputs
KpFuse = 0.0025; % fuselage density factor (Table 10.11)

% Fuselage Material Composition
Fuse.dens = dens.Al7075; % Aluminum 7075 density (lb/in^3)

% Fuselage Equation Parameters
Fuse.Length = 120; % fuselage length (ft)
Fuse.DiameterMax = 8; % max diameter of equivalent circular cross-section (ft)
KInlet = 1.25; % 1.25 for inlets on fuselage or 1 for inlets elsewhere

% Fuselage Weight
Wt.Fuselage = Fuse.Length * (Fuse.DiameterMax^2) * Fuse.dens * KpFuse * (nUlt^0.25) * KInlet * g; 
%% ---------------------------- Landing Gear --------------------------- %%

% Section: 10.4.5
% Eqn(s) 10.8

% Equation Parameters
KL = 1; % landing place factor (1.8 for Navy aircraft and 1 otherwise)
Kret = 1.07; % 1 for fixed landing gear and 1.07 for retractable landing gear
KLandingGear = 0.28; % landing gear weight factor (Table 10.12)
Wt.Landing = 70000; % aircraft weight at landing (lbf)
LG.Height = 10; % landing gear height (ft)
Wing.span = 60; % wing span (ft)

% Landing Gear Weights
Load.MainGear = 0.90; % percentage of load absorbed by main landing gear
Load.NoseGear = 0.10; % percentage of load absorbed by nose landing gear

% Weight of Nose Gear + Main LG (Tricycle Configuration)
Wt.LandingGear = KL * Kret * KLandingGear * Wt.Landing * (LG.Height / Wing.span) * (nUlt^0.2);

% Weight of Main Landing Gear
Wt.MainGear = Load.MainGear*Wt.LandingGear;

% Weight of Nose Landing Gear
Wt.NoseGear = Load.NoseGear*Wt.LandingGear;

%% -------------------------- Installed Engine  ------------------------ %%

% Section: 10.4.6
% Eqn(s) 10.9

% For installed engine weight of other types of aircraft 
% (e.g., fighter, transport), the interested reader is referred to Ref. [5]. 
% For fighter and transport aircraft, the propulsion system weight includes weight
% of engine cooling, weight of starter, weight of inlet system, weight of firewall, weight of
% nacelle, weight of engine control, and weight of auxiliary power unit.
% Equation Parameters
%KE = 2.6; % engine weight factor (2.6 for British units, 3 for metric units)
%NE = 3; % number of engines
%WtEngine = 3960; % weight of single engine (lbf)

% Installed Engine Weight
%WtInstalledEngine = KE * NE * (WtEngine^0.9);

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
