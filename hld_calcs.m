clear;
clc;
close all;

load('aircraft_vars.mat');

%% FLAP/SLAT CALCULATIONS FOR TAKEOFF

% Will be using fowler flap and slats to get a wing that will allow us to
% takeoff

%% LIFT REQUIREMENTS & DEFICIT
% Determine the lift coefficient needed

% WEIGHT: Wt.WTO
% REQUIRED THRUST: constraints.req_Thr

hld.gamma_TO = atan(35/100); % FAR 25 requirement for 35 ft obstacle, and with 100 flight path to clear it
hld.alpha_TO = 10 * pi / 180; % take-off AOA, deg -> radians
hld.V_TO = constraints.Vstall * 1.2; % fly at 20% safety margin 
hld.Cl_base = 0.61; % generated w/ OpenVSP
hld.sf = 1.2; % safety factor to obtain Cl_max
hld.dCl_req = (hld.sf*Wt.WTO * cos(hld.gamma_TO) - constraints.req_Thr * sin(hld.alpha_TO))/(0.5*atm.rho_sl*(hld.V_TO^2)*WING.geom.S_area) - hld.Cl_base;
fprintf('Cl deficit: %0.5f\n', hld.dCl_req);

%% DATCOM 1978: HLD for supersonic airfoil 

% p. 8-13 
dcl_max_b = 1.05; %for t/c = 5.25 
k1 = 1.05; % for flap-chord of 25%
k2 = 1; % 40 degree flap deflection
k3 = 1; % assumes actual flap angle equals reference flap angle
dcl_max_f = k1*k2*k3*dcl_max_b; % trailing edge flaps c_l
% p. 8-18
%reference--> http://adg.stanford.edu/aa241/highlift/clmaxest.html
K_lambda = 0.92; % corresponding to wing sweep of 10 degrees
S_RAT_f = 0.470; % flaps-wing wetted area ratio, referenced to DC-6 which had ratio of 0.58
dCL_max_f = dcl_max_f * S_RAT_f * K_lambda; %  change in coefficient of lift from flaps

% FLAPERON
% p. 8-13 
dcl_max_b = 1.05; %for t/c = 5.25 
k1 = 1.05; % for flap-chord of 25%
k2 = 0.8; % 20 degree flap deflection
k3 = 1; % assumes actual flap angle equals reference flap angle
dcl_max_f = k1*k2*k3*dcl_max_b; % trailing edge flaps c_l
% p. 8-18
%reference--> http://adg.stanford.edu/aa241/highlift/clmaxest.html
K_lambda = 0.92; % corresponding to wing sweep of 10 degrees
S_RAT_f = 0.40; % flaps-wing wetted area ratio, referenced to DC-6 which had ratio of 0.58
dCL_max_flaperon = dcl_max_f * S_RAT_f * K_lambda; %  change in coefficient of lift from flaps


dC_L1 = dCL_max_f ;%+ d_CL_max_s;
fprintf('dCL = %f for supersonic wing \n',dC_L1)
fprintf('dCL = %f for supersonic wing flaperon \n',dCL_max_flaperon)
fprintf('Total dCL = %f  \n',dCL_max_flaperon+dC_L1)

% cl_d_max = 1.7; % with nose flap-chord ratio, cf/c = 0.30
% n_max = 0.45; % corresponding to 0 leading edge radius
% n_delta = 0.58; % corresponding to 30 degree deflection angle
% delta_f = 20; % degree deflection from airfoil chord plane
% chord_RAT = .10; % ratio of the chord with and without deflection of slat
% dcl_max_s = cl_d_max * n_max * n_delta * delta_f * (chord_RAT); % lead edge slats c_l
% S_RAT_s = 0.80; % slat-wing wetted area ratio 
% phi_HL = 10; %sweep angle of the hinge line of flap
% d_CL_max_s = dcl_max_s * S_RAT_s * cosd(phi_HL);

% For Reference...
% dC_L = 1.55 for DC-9-30
%http://adg.stanford.edu/aa241/highlift/highliftintro.html

%percentFlap = dCL_max_f/dC_L1
%percentSlat = 1-percentFlap

%% DATCOM 1978: HLD for subsonic airfoil 

% p. 8-13 
% dcl_max_b = 1.22; %for t/c = 10 
% k1 = 1.0; % for flap-chord of 30%
% k2 = 1; % 40 degree flap deflection
% k3 = 1; % assumes actual flap angle equals reference flap angle
% dcl_max_f = k1*k2*k3*dcl_max_b; % trailing edge flaps c_l
% % p. 8-18
% %reference--> http://adg.stanford.edu/aa241/highlift/clmaxest.html
% K_lambda = 0.57; % corresponding to wing sweep of 60 degrees
% S_RAT_f = 1; % flaps-wing wetted area ratio, referenced to DC-6 which had ratio of 0.58
% dCL_max_f = dcl_max_f * S_RAT_f * K_lambda; %  change in coefficient of lift from flaps

% cl_d_max = 1.6; % with nose flap-chord ratio, cf/c = 0.30
% n_max = 0.45; % corresponding to 0 leading edge radius
% n_delta = 0.58; % corresponding to 30 degree deflection angle
% delta_f = 30; % degree deflection from airfoil chord plane
% chord_RAT = .25; % ratio of the chord with and without deflection of slat
% dcl_max_s = cl_d_max * n_max * n_delta * delta_f * (chord_RAT); % lead edge slats c_l
% S_RAT_s = 0.20; % slat-wing wetted area ratio 
% phi_HL = 60; %sweep angle of the hinge line of flap
% d_CL_max_s = dcl_max_s * S_RAT_s * cosd(phi_HL);

dC_L2 = dCL_max_f; %+ d_CL_max_s;
%fprintf('dCL = %f for subsonic wing \n',dC_L2)
% For Reference...
% dC_L = 1.55 for DC-9-30
%http://adg.stanford.edu/aa241/highlift/highliftintro.html

%percentf = dCL_max_f/dC_L2
%percents = 1-percentf

%% dCL for supersonic and subsonic wing combined

%k1 = 0.80; % percent of wetted area for supersonic airfoil region
%k2 = 0.20; % percent of wetted area for subsonic airfoil region

%dCL = k1*dC_L1 + k2*dC_L2;
%fprintf('dCL = %f for total wing   \n',dCL)
