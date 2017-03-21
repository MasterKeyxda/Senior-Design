clear all 
close all 
clc

<<<<<<< HEAD
load aircraft_vars.mat; 
close all
% Find: flap chord, flap span, flap deflection during TO and TD.


%V_stall = 211; % ft/s
%V_to = 1.2 * V_stall; % Sadraey p. 230,255

%C_Lto = 2 * WTO / (atm.rho_sl * V_to^2 * S_w );

% Lifting-line theory can be employed to calculate the lift increment for 
% each HLD deflection. You can then adjust the HLD chord (Cf) to achieve 
% the required lift increment. 

flapChordRAT = 0.30; % p. 237 from Gulfstream II

% recommended flap deflections, 20deg:TO, 50deg:LD
delta_f = 60; % flap deflection
alpha_dflap = -1.15 * flapChordRAT * delta_f;

%% Sadraey's MATLAB code, p. 256
N = 100; % (number of segments-1)
S = S_w*0.092903; % m^2
AR = e; % Aspect ratio
lambda = taperh; % Taper ratio
alpha_twist = 0.0001; % Twist angle (deg)
i_w = 0; % wing setting angle (deg)
a_2d = 2*pi/sqrt(1.6^2 - 1); % lift curve slope (1/rad)
a_0 = alpha_dflap; % flap up zero-lift angle of attack (deg)
a_0_fd = alpha_dflap*2; % flap down zero-lift angle of attack (deg)
b = sqrt(AR*S); % wing span
bf_b=1; %flap-to-wing span ratio
MAC = S/b; % Mean Aerodynamic Chord
Croot = (1.5*(1+lambda)*MAC)/(1+lambda+lambda^2); % root chord
theta = pi/(2*N):pi/(2*N):pi/2;
alpha=i_w+alpha_twist:-alpha_twist/(N-1):i_w;% segment?s angle of attack
for i=1:N
if (i/N)>(1-bf_b)
       alpha_0(i)=a_0_fd; %flap down zero lift AOA

       else
      alpha_0(i)=a_0; %flap up zero lift AOA
end
end
z = (b/2)*cos(theta);
c = Croot * (1 - (1-lambda)*cos(theta)); % MAC at each segment
mu = c * a_2d / (4 * b);
LHS = mu .* (alpha - alpha_0)/57.3; % Left Hand Side
% Solving N equations to find coefficients A(i):
for i=1:N
    for j=1:N
    B(i,j) =  sin((2*j-1) * theta(i)) * (1 + (mu(i) *(2*j-1)) / sin(theta(i)));
    end
end
A=B\transpose(LHS);
for i = 1:N
    sum1(i) = 0;
    sum2(i) = 0;
    for j = 1 : N
        sum1(i) = sum1(i) + (2*j-1) * A(j)*sin((2*j-1)*theta(i));
        sum2(i) = sum2(i) + A(j)*sin((2*j-1)*theta(i));
    end
end
   

CL_TO = pi * AR * A(1);

% Account for wing-sweep
omega = 35; % sweep angle
CL_TO = CL_TO / (cos(omega))^2
%http://adg.stanford.edu/aa200b/potential3d/sweeptheory.html

%% 

velocity = 230; % ft/s
Lift = 0.5 * CL_TO * atm.rho_sl * velocity^2 * S_w;
disp(Lift)

diff = WTO - Lift

%% AircraftDesign_8_HighLift, DATCOM

% p. 8-13 
% dcl_max_b = 1.05 for t/c = 5.25
% dcl_max_b = 1.22 for t/c = 10  
dcl_max_b = 1.05;
k1 = 1.2; 
k2 = 1; % 40 degreee flap deflection
k3 = 1; % needs to be changed
dcl_max_f = k1*k2*k3*dcl_max_b; % trailing edge flaps c_l

cl_d_max = 1.7; % with nose flap-chord ratio, cf/c = 0.30
n_max = 0.45; % corresponding to 0 leading edge radius
n_delta = 0.58; % corresponding to 30 degree deflection angle
delta_f = 30; % degree deflection from airfoil chord plane
chord_RAT = 0.35; % ratio of the chord with and without deflection of slat
dcl_max_s = cl_d_max * n_max * n_delta * delta_f * (chord_RAT); % lead edge slats c_l


% p. 8-18
K_lambda = 0.82; % corresponding to wing sweep of 35 degrees
S_RAT_f = 0.6; % area ratio for the flaps
dCL_max_f = dcl_max_f * S_RAT_f * K_lambda; %  change in coefficient of lift from flaps

S_RAT_s = 0.30; % area ratio of the slats
phi_HL = 30; %sweep angle of the hinge line of flap
d_CL_max_s = dcl_max_s * S_RAT_s * cos(phi_HL)

CL_maxc = 1.2;

dC_L = dCL_max_f + d_CL_max_s
%C_Lmax = CL_maxc + dCL_max_f + d_CL_max_s










=======
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

%cl_d_max = 1.6; % with nose flap-chord ratio, cf/c = 0.20
%n_max = 0.45; % corresponding to 0 leading edge radius
%n_delta = 0.58; % corresponding to 20 degree deflection angle
%delta_f = 30; % degree deflection from airfoil chord plane
%chord_RAT = .20; % ratio of the chord with and without deflection of slat
%dcl_max_s = cl_d_max * n_max * n_delta * delta_f * (chord_RAT); % lead edge slats c_l

% From DATCOM 1978, p.1885,1920
dcl_s = 0.020; % for cf/c = 0.50
delta_f = 30; % degrees
c_RAT = 1.30;
dcl_max_s = dcl_s * delta_f * c_RAT;

% p. 8-18
%reference--> http://adg.stanford.edu/aa241/highlift/clmaxest.html
K_lambda = 0.92; % corresponding to wing sweep of 10 degrees
S_RAT_f = 0.65; % flaps-wing wetted area ratio, referenced t?o DC-6 which had ratio of 0.58
dCL_max_f = dcl_max_f * S_RAT_f * K_lambda; %  change in coefficient of lift from flaps

S_RAT_s = 1; % slat-wing wetted area ratio 
phi_HL = 10; %sweep angle of the hinge line of flap
dCL_max_s = dcl_max_s * S_RAT_s * cosd(phi_HL);

dC_L1 = dCL_max_f + dCL_max_s;
fprintf('dCL = %f for supersonic wing \n',dC_L1)
% For Reference...
% dC_L = 1.55 for DC-9-30
%http://adg.stanford.edu/aa241/highlift/highliftintro.html

<<<<<<< Updated upstream
<<<<<<< Updated upstream
%percentFlap = dCL_max_f/dC_L1
%percentSlat = 1-percentFlap
=======

>>>>>>> Stashed changes
=======

>>>>>>> Stashed changes

%% DATCOM 1978: HLD for supersonic airfoil 

% p. 8-13 
dcl_max_b = 1.22; %for t/c = 10 
k1 = 1.2; % for flap-chord of 30%
k2 = 1; % 40 degree flap deflection
k3 = 1; % assumes actual flap angle equals reference flap angle
dcl_max_f = k1*k2*k3*dcl_max_b; % trailing edge flaps c_l

%cl_d_max = 1.7; % with nose flap-chord ratio, cf/c = 0.30
%n_max = 0.45; % corresponding to 0 leading edge radius
%n_delta = 0.58; % corresponding to 30 degree deflection angle
%delta_f = 30; % degree deflection from airfoil chord plane
%chord_RAT = .30; % ratio of the chord with and without deflection of slat
%dcl_max_s = cl_d_max * n_max * n_delta * delta_f * (chord_RAT); % lead edge slats c_l

% From DATCOM 1978, p.1885,1920
dcl_s = 0.020; % for cf/c = 0.50
delta_f = 30; % degrees
c_RAT = 1.30;
dcl_max_s = dcl_s * delta_f * c_RAT;


% p. 8-18
%reference--> http://adg.stanford.edu/aa241/highlift/clmaxest.html
K_lambda = 0.57; % corresponding to wing sweep of 60 degrees
S_RAT_f = 0.60; % flaps-wing wetted area ratio, referenced to DC-6 which had ratio of 0.58
dCL_max_f = dcl_max_f * S_RAT_f * K_lambda; %  change in coefficient of lift from flaps

S_RAT_s = 1; % slat-wing wetted area ratio 
phi_HL = 65; %sweep angle of the hinge line of flap
d_CL_max_s = dcl_max_s * S_RAT_s * cosd(phi_HL);

dC_L2 = dCL_max_f + d_CL_max_s;
fprintf('dCL = %f for subsonic wing \n',dC_L2)
% For Reference...
% dC_L = 1.55 for DC-9-30
%http://adg.stanford.edu/aa241/highlift/highliftintro.html

<<<<<<< Updated upstream
<<<<<<< Updated upstream
%percentf = dCL_max_f/dC_L2
%percents = 1-percentf
=======
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

%% dCL for supersonic and subsonic wing combined

k1 = 0.80; % percent of wetted area for supersonic airfoil region
k2 = 0.20; % percent of wetted area for subsonic airfoil region

dCL = k1*dC_L1 + k2*dC_L2;
fprintf('dCL = %f for total wing   \n',dCL)
>>>>>>> master

