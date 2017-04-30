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
n_p = 10; % number of passengers
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
W_s = 22/12 ; % seat width ft
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

L_R = 98.33-3.5; % arbitrary length of rear section ft


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
m_f = 213514/9.81; % kg of fuel
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
L_F2 = L_F; % Length of Fuselage
% Looking at the SEARS-HACK Body Equation from WIKIPEDIA
x = linspace(0,1,100); % range
x_2 = linspace(0,150*12,100);
R_max = D_C/2; % maximum radius of Sears Haack
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

ylabel('Y (ft)')
xlabel('X (ft)')
title('Top View Fuselage Layout')
text(32,0,'Passenger Cabin')
text(26,0,'CP')
text(15,0,'Nose')
text(75,0,'Rear Section')
%title('Haack Body Fuselage Design')

% close all;
Xloc = [0, 13.33155, 100/3, 41.66667, 66.667, 100, 120.83333, 129.1667, 131.25, 141.6667, L_F2];
Z_UP = [0, 1.16652, 10/3, 3.62475, 10/3, 3.46122, 3.45386, 3.39437, 3.36786, 3.52812, 0];
ZX1_u = (Z_UP(2) - Z_UP(1))/(Xloc(2) - Xloc(1));
ZXn_u = (Z_UP(end) - Z_UP(end-1))/(Xloc(end) - Xloc(end-1));
Zdat_u = spline(Xloc, [ZX1_u, Z_UP, ZXn_u]);

Z_LOW = [0.0, 0.91758, 10/3, 4.65953, 5.00052, 5.03238, 4.98319, 4.98680, 4.99473, 3.77775, 0]; 
ZX1_l = (Z_LOW(2) - Z_LOW(1))/(Xloc(2) - Xloc(1));
ZXn_l = (Z_LOW(end) - Z_LOW(end-1))/(Xloc(end) - Xloc(end-1));
Zdat_l = spline(Xloc, [ZX1_l Z_LOW ZXn_l]);

figure();plot(1:L_F2, spline(Xloc, [ZX1_u, Z_UP, ZXn_u], 1:L_F2));
hold on;
plot(1:L_F2, -1.*spline(Xloc, [ZX1_l Z_LOW ZXn_l], 1:L_F2));
axis equal

edit_vsp = 0; % change this one if you want to modify the VSP file!

if edit_vsp == 1
    for i = 1:length(Xloc)
        [R_val, R_nx, R_nxx] = fuselage_geom(Xloc(i), R_max, L_F2);
       fprintf('X = %0.5f, Station %i\n', Xloc(i), i-1);
       fprintf('Height: %0.5f\n', Z_UP(i) + Z_LOW(i));
       fprintf('Z_Disp: %0.5f\n', Z_UP(i)-0.5*(Z_UP(i)+Z_LOW(i)));
       fprintf('MaxWLoc: %0.5f\n', 1 - Z_UP(i)/(0.5*(Z_UP(i)+Z_LOW(i))));
       fprintf('\n');
       fprintf('Radius (Right) = %0.5f\n', R_val*2);
       fprintf('Angle (Right) = %0.5f\n', R_nx);
       fprintf('Curvature (Right) = %0.5f\n', R_nxx);
       fprintf('\n');
       if i == length(Xloc)
    %        fprintf('Radius (Top) = %0.5f\n', Z_UP(i-1));
           fprintf('Angle (Top) = %0.5f\n', atand(Zdat_u.coefs(i-1,3)));
           fprintf('Curvature (Top) = %0.5f\n', 2*Zdat_u.coefs(i-1,2));
           fprintf('\n');
    %        fprintf('Radius (Bottom) = %0.5f\n', Z_UP(i-1));
           fprintf('Angle (Bottom) = %0.5f\n', atand(Zdat_l.coefs(i-1,3)));
           fprintf('Curvature (Bottom) = %0.5f\n', 2*Zdat_l.coefs(i-1,2));
       else
           fprintf('Angle (Top) = %0.5f\n', atand(Zdat_u.coefs(i,3)));
           fprintf('Curvature (Top) = %0.5f\n', 2*Zdat_u.coefs(i,2));
           fprintf('\n');
    %        fprintf('Radius (Bottom) = %0.5f\n', Z_UP(i));
           fprintf('Angle (Bottom) = %0.5f\n', atand(Zdat_l.coefs(i,3)));
           fprintf('Curvature (Bottom) = %0.5f\n', 2*Zdat_l.coefs(i,2));
       end

       fprintf('\n\n');
    end
end
% [R_cp, R_cpx, R_cpxx] = fuselage_geom(x_CP, R_max, L_F2);
% [R_c, R_cx, R_cxx] = fuselage_geom(x_C, R_max, L_F2);
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
%% Egress
% By FAR 25 Rules, aircraft with passengers between 10 and 19
% passengers must have emergency egress of rectangular shape no less than
% 20" wide and 36" tall with stepup no greater that 20" (27" if above wing)
