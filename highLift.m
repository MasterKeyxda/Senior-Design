clear all 
close all 
clc

%This code calculates the change in C_L from flaps and slats.
%User must change k1 and k2 which represent the percentages of wetted area
%for the subsonic and supersonic airfoils.

%Followed procedure from 'AircraftDesign_8_HighLift'. 
%File located in Box folder named 'High Lift Devices'. 
%'AircraftDesign_8_HighLift' document uses plots from DATCOM (1978).


%% DATCOM 1978: HLD for supersonic airfoil 

% p. 8-13 
dcl_max_b = 1.05; %for t/c = 5.25 
k1 = 1.2; % for flap-chord of 30%
k2 = 1; % 40 degree flap deflection
k3 = 1; % assumes actual flap angle equals reference flap angle
dcl_max_f = k1*k2*k3*dcl_max_b; % trailing edge flaps c_l

cl_d_max = 1.7; % with nose flap-chord ratio, cf/c = 0.30
n_max = 0.45; % corresponding to 0 leading edge radius
n_delta = 0.58; % corresponding to 30 degree deflection angle
delta_f = 30; % degree deflection from airfoil chord plane
chord_RAT = .30; % ratio of the chord with and without deflection of slat
dcl_max_s = cl_d_max * n_max * n_delta * delta_f * (chord_RAT); % lead edge slats c_l

% p. 8-18
%reference--> http://adg.stanford.edu/aa241/highlift/clmaxest.html
K_lambda = 0.92; % corresponding to wing sweep of 10 degrees
S_RAT_f = 0.58; % flaps-wing wetted area ratio, referenced to DC-6 which had ratio of 0.58
dCL_max_f = dcl_max_f * S_RAT_f * K_lambda; %  change in coefficient of lift from flaps

S_RAT_s = 0.20; % slat-wing wetted area ratio 
phi_HL = 10; %sweep angle of the hinge line of flap
d_CL_max_s = dcl_max_s * S_RAT_s * cosd(phi_HL);

dC_L1 = dCL_max_f + d_CL_max_s;
fprintf('dCL = %f for supersonic wing \n',dC_L1)
% For Reference...
% dC_L = 1.55 for DC-9-30
%http://adg.stanford.edu/aa241/highlift/highliftintro.html

%percentFlap = dCL_max_f/dC_L1
%percentSlat = 1-percentFlap

%% DATCOM 1978: HLD for supersonic airfoil 

% p. 8-13 
dcl_max_b = 1.22; %for t/c = 10 
k1 = 1.2; % for flap-chord of 30%
k2 = 1; % 40 degree flap deflection
k3 = 1; % assumes actual flap angle equals reference flap angle
dcl_max_f = k1*k2*k3*dcl_max_b; % trailing edge flaps c_l

cl_d_max = 1.7; % with nose flap-chord ratio, cf/c = 0.30
n_max = 0.45; % corresponding to 0 leading edge radius
n_delta = 0.58; % corresponding to 30 degree deflection angle
delta_f = 30; % degree deflection from airfoil chord plane
chord_RAT = .30; % ratio of the chord with and without deflection of slat
dcl_max_s = cl_d_max * n_max * n_delta * delta_f * (chord_RAT); % lead edge slats c_l

% p. 8-18
%reference--> http://adg.stanford.edu/aa241/highlift/clmaxest.html
K_lambda = 0.57; % corresponding to wing sweep of 60 degrees
S_RAT_f = 0.58; % flaps-wing wetted area ratio, referenced to DC-6 which had ratio of 0.58
dCL_max_f = dcl_max_f * S_RAT_f * K_lambda; %  change in coefficient of lift from flaps

S_RAT_s = 0.20; % slat-wing wetted area ratio 
phi_HL = 60; %sweep angle of the hinge line of flap
d_CL_max_s = dcl_max_s * S_RAT_s * cosd(phi_HL);

dC_L2 = dCL_max_f + d_CL_max_s;
fprintf('dCL = %f for subsonic wing \n',dC_L2)
% For Reference...
% dC_L = 1.55 for DC-9-30
%http://adg.stanford.edu/aa241/highlift/highliftintro.html

%percentf = dCL_max_f/dC_L2
%percents = 1-percentf

%% dCL for supersonic and subsonic wing combined

k1 = 0.80; % percent of wetted area for supersonic airfoil region
k2 = 0.20; % percent of wetted area for subsonic airfoil region

dCL = k1*dC_L1 + k2*dC_L2;
fprintf('dCL = %f for total wing   \n',dCL)

