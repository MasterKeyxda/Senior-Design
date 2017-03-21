%% Fuselage Preliminary Design % Sadraey %
% Meta
% FUSELAGE DESIGN CHANGED TO SEARS-HAACK BODY
% Total Length Extended to 150ft
%%
% Test Change in GITHUB
% Design Requirements 
% Main: Accomodate Passengers
% Secondary: Must accomodate crew, landing gear, fuel, systems, empennage
% Optimization: Low Drag, Positive Lift generation, low weight, comfort,
% structural flight, external symmetry, low wetted area
%% Load Aircraft Variables
if ~exist('Wt', 'var')
   clc;
   clear;
   load('aircraft_vars.mat'); 
   close all;
end

%% Determine Number of Passengers and Crew
n_p = 8; % number of passengers
n_c = 2; % number of crew members
%% Establish Human Size and Target Passenger
% Passengers will be from U.S.
% Figures 7.5 and 7.6
% Design for passengers that are 5' 2" - 6' 3"
% Arm length 2.6 ft
% Chair Height 1.6 ft
% Shoulder to Shoulder 2ft
%% Cockpit % Sadraey %
% Distance from spine to taxi instrument 3.5 ft
% 1ft minimum from butt to floor
% 4 - 10 cm
T_w = 0.23; % ft with 7 cm wall thickness
L_CP = 5.5; %ft preliminary choice
W_pilot = 15/12; % distance between pilot seats (aisle) ft 
% Table 7.4
W_CPs = 20/12; % seat width cockpit ft
W_CP = 2*W_CPs + W_pilot; % Cockpit Width interior
% Equation 7.9 
D_CP = W_CP + 2*T_w; % diameter of cockpit with thickness included
% Must look into instrumentation sizing
fprintf('Cockpit Diameter is %0.2f ft \n',D_CP)
fprintf('Cockpit Length is %0.2f ft \n',L_CP)
%% Passenger Cargo
% pg 369 Eq 7.4
v_b = 5.156; % average bag/luggage size ft^3
v_btotal = n_p * v_b; %ft^3


%% L/D length to diameter ratio
% Table 7.7
% L/D for Concorde is 23.
% shooting for this range
%% Passenger Cabin % Sadraey %
% Calculate Cabin Diameter
W_s = 25/12 ; % seat width ft
% Table 7.3
W_a = 23/12; % aisle distance ft
W_C = 2*W_s + W_a; % Interior Cabin Width ft
% Equation 7.9
D_C = W_C + 2*T_w; % Cabin Diameter ft

% Calculate Cabin Length
P_s = 40/12; % seat pitch or distance between seats (ft)

n_r = n_p/2; % number of rows
% Equation 7.1
L_C = n_r*P_s+8; % Cabin Length
fprintf('Passenger Cabin Diameter is %0.2f ft \n',D_C)
fprintf('Passenger Cabin Length is %0.2f ft \n',L_C)
% Headroom from Table 7.4
% Will we have carry-on over or same level?
Headroom = 66.912; %ft for first class
% Must also consider lavatory

%% Cargo
% Depends on volume and will most likely be beneath passenger cabin
%% Nose Design
% Look at supersonic design and weight and balance
L_n = 25; % arbitrary length of nose ft
% Looking at the SEARS-HACK Body Equation from WIKIPEDIA
% x = linspace(0,0.5,100); % range
% R_max = 2; % arbitrary maximum nose radius ft
% V_n = (3*pi^2/16 * R_max^2*L_n)/2; % Volume of nose
% S_x = pi*R_max^2*(4*x.*(1-x)).^(3/2); % nose cross sectional area
% r_x = R_max*(4*x.*(1-x)).^(3/4);% nose radius
% FigHandle = figure('Position', [250, 250, L_n*40, R_max*40]);
% figure(1)
% plot(2*x*L_n,r_x,'r',2*x*L_n,-r_x,'r')
% title('Sears-Hack Body Nose Design')
% Based on Wave Drag Equation, wave drag is reduced by increasing length
% and decreasing volume
% calculated from the derivative of the radius
% Theta = 16.4223 degrees
% Beta = 58.471 degrees
% Oblique Shock Ma2 = 1.06
%% Rear Section
% Will Include Tail and Wing as well for supersonic aircraft.

L_R = 99; % arbitrary length of rear section ft


%% Length to Diameter Ratio Fineness Ratio
L_N = 25; % arbitrary length of nose ft
%% Length of Fuselage
% Equation 7.10 
L_F = L_CP + L_C + L_N + L_R; % Total Fuselage Length
fprintf('The fuselage length is %0.2f ft \n',L_F)
%% Length to Diameter Ratio
Dmax = 100/12;
LF_DF = L_F/Dmax; 
fprintf('The L/D ratio for the fuselage is %0.2f  \n',LF_DF)
%% Fuel Tanks
m_f = 182377.09/9.81; % kg of fuel
% Table 7.8 for fuel types and density
rho_f = 840; % kg/m^3
V_f = m_f/rho_f; % m^3
V_f = V_f * 35.31; % ft^3
fprintf('The fuel tank volume is %0.2f ft^3 \n',1.2*V_f)
%% Volume of Bottom Half
% Equation 7.4a
V_bot = 0.5* (pi*W_C^2/4*(L_C+L_R));
fprintf('The bottom volume for the fuselage is %0.2f ft^3 \n',V_bot)
%% SEARS-HACK FUSELAGE
% Look at supersonic design and weight and balance
L_F2 = L_F; % arbitrary length of nose ft
% Looking at the SEARS-HACK Body Equation from WIKIPEDIA
x = linspace(0,1,100); % range
x_2 = linspace(0,140*12,100);
R_max = D_C/2; % arbitrary maximum nose radius ft
V_n = (3*pi^2/16 * R_max^2*L_F2)/2;
S_x = pi*R_max^2*(4*x.*(1-x)).^(3/2); % nose cross sectional area
r_x = R_max*(4*x.*(1-x)).^(3/4);
r_x2 = ((1/(12.*L_F2)).*(4.*x_2.*(12.*L_F2-x_2)).^.75)/12; % Solidworks
figure('Position', [250, 250, L_F2*40, R_max*40]);
%figure(1)
x_CP = L_n+L_CP;
x_C = x_CP + L_C;
plot(x*L_F2,r_x2,'r',x*L_F2,-r_x2,'r')
line([L_n L_n], [-3.25 3.25]);
line([x_CP x_CP], [-3.25 3.25]);
line([x_C x_C], [-3.25 3.25]);

%title('Haack Body Fuselage Design')
%%
% z = L_F*x;
% Z = meshgrid(z,z);
% Theta = linspace(0,2*pi,100);
% [R ,PHI] = meshgrid(r_x,Theta);
% FigHandle = figure('Position', [50, 50, 600, 600]);
% surf(R.*cos(PHI), R.*sin(PHI), Z);
% view(3)
% set(gca,'zdir','reverse')
% colormap hsv
% axis equal

