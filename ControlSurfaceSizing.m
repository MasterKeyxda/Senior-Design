    %clear all
%close all
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
I_xx = 6730260.3;
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
t = sqrt((2*phi_des)/(P_dot)); % time in seconds

% Geometry of each aileron
C_w = 16.6; % MAC
b_A = y_out - y_in; % span of aileron
c_A = 0.3*C_w; % chord of aileron
A_A = 2*b_A*c_A; % area of both ailerons 

fprintf('Time to rotate 30 degrees (s) = %0.5f \n',t)
fprintf('Aileron Span (ft) = %0.5f \n',b_A)
fprintf('Aileron Chord (ft) = %0.5f \n',c_A)
fprintf('Aileron Area  (ft^2)  = %0.5f \n\n',A_A)

% Plot: Variations of bank angle versus time
n = 50;
phi2 = linspace(0,40,n);
t2 = sqrt((2*phi2)/(P_dot));
xx1 = t*ones(n);
yy1 = linspace(0,30,n);
xx2 = linspace(0,t,n);
yy2 = 30*ones(n);

figure
plot(t2,phi2,'o', xx1,yy1,'r',xx2,yy2,'r')
xlabel('time (sec)','fontsize',20)
ylabel('\phi (deg)','fontsize',20)

%% Elevator Sizing
% Satisifes take-off angular acceleration requirement of 7 deg/s^2 for small
% transport 

% Follows Sadraey's elevator-chord and span ratio recommendations
% Longitudinal Trim requirements for given flight envelope


C_L_cruise = 2*85000 / (0.00067595*(0.8*968)^2*825);
C_L_TO = 1.75;
C_D_TO = 0.0458;
V_s = 211; %ft/s

%Longitudinal Aerodynamic Forces
Drag_TO = 0.5*0.0002378*V_s^2*S*C_D_TO;
Lift_TO = 0.5*0.0002378*V_s^2*S*C_L_TO;
M_ac_wf =  0.5*0.002378*211^2*825*-.0085*16.6;

mu = 0.04; % friction coefficient of runway. assume concrete runway
W = 87000;
F_runway = mu*(W - Lift_TO);
mass = W; % lbf
T = 21900*3; % lbf

% Aircraft Linear Acceleration at the time of take-off rotation
a = (T - Drag_TO - F_runway) / mass;

% Step 7: Calculation of the contributing pitching moments in the take-off
% rotation
x_mg = 94.75; % ft
x_cg = 88.5; % ft
z_Dmg = 9.54; 
z_Tmg = 10.83;
x_mgacwf = 9.31;
z_cgmg = 10.08;
x_ac_h = 104; % UPDATE
z_mg = 10; %UPDATE
z_T = z_Tmg + z_mg;

rho = 23.77*10^(-4); %sea level
v = 231; % ft/s
A = WING.geom.S_area; %Reference area
% must run iter_weights
D = 0.5 * Cd_TO * rho * v^2 * A;
L_wf = C_L_TO * 0.5 * rho * v^2 * A;

M_W = W*(x_mg - x_cg);
M_D = D*(z_Dmg);
M_T = T*(z_Tmg);
M_L_wf = L_wf*(x_mgacwf);
M_a = m*a*(z_cgmg);

% take off pitch angular accleration for small transport aircraft,  
% rad/s^2, Table 12.9
theta_ddot = 7/(180/pi); 

I_yymg = 21000000.0;
% I_yymg = 2.1e7;
V_R = v; % aircraft linear speed at instant of rotation at TO             
L_h = (L_wf*(x_mgacwf) + M_ac_wf + m*a*(z_cgmg) + W*(x_mg-x_cg) + ...
    D*(z_Dmg) + T*(z_Tmg) - I_yymg*theta_ddot) / (x_ac_h - x_mg);
C_L_h = (2*L_h) / (rho*V_R^2*S_h);

C_L_alpha_h = TAIL.CLalphah;
delta_E_max = -25 / (180/pi); % recommended maximum deflection angle, Table 12.3
alpha_h = TAIL.hAngle; %from iter_weights
tau_e = (alpha_h + (C_L_h/C_L_alpha_h)) / (delta_E_max);

% Step 11
CE_Ch = 0.3; % from tau from Figure 12.12
Dalpha_oe = -1.14*CE_Ch*(-25);

% Check if elevator deflection stalls the horizontal tail during take-off
% rotation
alpha_S_TO = 10; % degrees
alpha_TO = alpha_S_TO - 2; % assumed fueselage is lifted up 2 deg. below the wing stall angle
%Horizontal Tail Take-off Angle with elevator
alpha_h_TO = alpha_TO * (1 - TAIL.deps) + TAIL.ih - TAIL.eps0;
%compare with Tail stall angle of attack during takeoff rotation
alpha_hs_deltaE0 = 13; % horizontal tail stall angle with no deflection
delta_alpha_hE = 5.3; % Table 12.19
% horizontal tail stall angle during takeoff rotation
alpha_h_s = (alpha_hs_deltaE0 - delta_alpha_hE); 
e_diff = alpha_h_s-alpha_h_TO;
fprintf('Horizontal tail take-off angle is less than\n the tail stall angle, thus the elevator is\n acceptable = %0.5f \n',e_diff)


l_h = x_mg;
% Step 16
V_H = (l_h * S_h) / (S* WING.geom.MAC);

bE_bh = 1;
eta_h = 0.9 ; %horizontal tail efficiency ,horizontal tail dynamic pressure ratio
C_mdeltaE = -TAIL.CLalphah * eta_h * V_H * bE_bh * tau_e; % 1/rad
C_LdeltaE = TAIL.CLalphah * eta_h * S_h/S*bE_bh * tau_e;
C_LhdeltaE = TAIL.CLalphah * tau_e; % 1/rad

C_L0 = 0.127;
C_L1 = C_L_cruise;
C_mo = 0.05; %UPDATE
rho = 23.77*10^(-4); %sea level
v = 231; % ft/s
q_bar = 0.5 * rho*v^2;
C_malpha = TAIL.Cm_alpha;
delta_E = - (((T*z_T/(q_bar * S * WING.geom.MAC) + C_mo)*C_L_alphaWing...
    + (C_L1-C_L0)*C_malpha) / (C_L_alphaWing*C_mdeltaE -C_malpha*C_LdeltaE));

b_E = bE_bh * TAIL.bh;
C_h = TAIL.ch;

% SOLUTION
C_E = CE_Ch * C_h;
S_E = b_E * C_E;
fprintf('Elevator chord (ft) = %0.5f \n',C_E)
fprintf('Elevator Span (ft^2) = %0.5f \n',b_E)
fprintf('Elevator Area (ft^2) = %0.5f \n',S_E)



