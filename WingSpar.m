clear all
close all
clc

%% Wing Spar Preliminary Calculations

% Wing Properties
span = 50 * 12; % inches
W_to = 102000; %lbf
force = W_to/2; % point load for one wing
y = span / 2; % inches, location of point load on beam

% Beam Properties
% Assumption: solid rectagular beam
h = 36; % height of beam , inches
b = h / 5; % base of beam, inches
I_x = (b * h^3) /12; % rectangular cross section

%Load Factors
n_lim = 2.50; % from v-n diagram
n_ult = 1.5 * n_lim;

% Shear Stress due to Bending
moment = force * y; %lbf.in
sigma = -(moment * y) / I_x;
sigma = sigma * n_ult / 3; % divide by number of spars
%"Aluminum is the most common material from which
%to construct wing" - FAA handbook

sigma_yield = 220000; %psi, Aircraft steel (5 Cr-M-V)
safetyFactor = sigma_yield / sigma;

fprintf('Bending Moment = %flbf.in \n', moment)
fprintf('Shear stress due to bending (1 spar) = %f psi \n',-sigma*3)
fprintf('Shear stress due to bending (3 spars) = %f psi \n',-sigma)
fprintf('Safety Factor including ultimate load factor (3 spars) = %f \n', safetyFactor)



