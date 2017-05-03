clear all
close all
clc

%% Aileron Design
% Aircraft must comply with maneuverability and roll requirements set by
% MIL-STD. time to achieve a specified bank angle change

% Given
m_TO = 87000; %lbs
S = 825; % ft^2
S_w = S;
b = 49.8;  %span
AR = 3;
lambda = 0.6; % taper ratio
S_h = 164.6;
S_vt = 74;
V_s = 211 ;% Vstall ft/s
C_L_alphaWing = 2*pi; % per rad
I_xx = 7040260.3;
tau = 0.48; %25 percent chord

y_out = b/2;
y_in = 0.60 * b/2;

c_bar = b / AR;
degrees = 20;
deltaA = degrees * (180/pi); %delta in RADIANS
rho = 23.77*10^(-4); %slugs/ft^3

C_r = 21.75;%3/2 * c_bar * (1+lambda)/(1 + lambda + lambda^2); % wing root chord

C_l_deltaA = (2*C_L_alphaWing*tau*C_r) / (S*b) * ...
    ((y_out^2/2 + 2/3*(lambda-1)/b*y_out^3)-(y_in^2/2 + 2/3*...
    (lambda-1)/b*y_in^3)); % 1/rad



C_l = C_l_deltaA * deltaA; % aircraft roll moment coefficient

V_app = 1.3*V_s;
L_A = 0.5*rho*V_app^2*S*C_l*b; % aircraft rolling moment

C_D_R = 0.9; % average value of wing, horizontal tail, vertical tail...
             % rolling drag coefficient 
y_D = 0.4* b/2; % drag momement arm assumed to be 40% of wing span
P_ss = sqrt((2*L_A)/(rho*(S_w + S_h + S_vt)*C_D_R*y_D^3)); % Steady state roll rate [rad/s]

phi_1 = I_xx / (rho*y_D^3*(S_w + S_h + S_vt)*C_D_R) * log(P_ss^2);

P_dot = P_ss^2 / (2*phi_1); % aircraft rate of roll rate

phi_des = 30; %degrees
t = sqrt((2*phi_des)/(P_dot)) % time in seconds

% Geometry of each aileron
C_w = 16.6 % MAC
b_A = y_out - y_in % span of aileron
c_A = 0.3*C_w % chord of aileron
A_A = 2*b_A*c_A % area of both ailerons 


%% Elevator Sizing

C_L_cruise = 2*85000 / (0.00067595*(0.8*968)^2*825);
C_L_TO = 1.75;

C_D_TO = 0.0458;
V_s = 211; %ft/s

Drag_TO = 0.5*0.0002378*V_s^2*S*C_D_TO;
Lift_TO = 0.5*0.0002378*V_s^2*S*C_L_TO;

M_ac_wf =  0.5*0.002378*211^2*825*-.0085*16.6;

mu = 0.04; % assume concrete runway
W = 87000;
F_runway = mu*(W - Lift_TO);

mass = W; % lbf
T = 21900*3; % lbf
a = (T - Drag_TO - F_runway) / mass

M_W = W*(x_mg - x_cg);
M_D = D*(z_D - z_mg);
M_T = T*(z_T - z_mg);
M_L_wf = L_wf*(x_mg - x_ac_wf);
M_a = m*a*(z_cg - z_mg);

theta_ddot = 7/(180/pi); % take off pitch angular accleration for small...
                % transport aircraft,  rad/s^2, Table 12.9
I_yymg = 
                
L_h = (L*(x_mg-x_ac_wf) + M_ac_wf + m*a*(z_cg-z_mg) + W*(x_mg-x_cg) + ...
    D*(z_D-z_mg) + T*(z_T-z_mg) - I_yymg*theta_ddot) / (x_ac_h - x_mg)
C_L_h = (2*L_h) / (rho*V_R^2*S_h);

delta_E_max = -25 / (180/pi); % recommended maximum deflection angle, Table 12.3
alpha_h = TAIL.hAngle; %from iter_weights, TAIL.AlphaVt
tau_e = (alpha_h + (C_L_h/C_L_alpha_h)) / (delta_E_max);

% Step 11



