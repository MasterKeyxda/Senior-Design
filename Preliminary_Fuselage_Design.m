%% Fuselage Preliminary Design % Sadraey %
% Design Requirements 
% Main: Accomodate Passengers
% Secondary: Must accomodate crew, landing gear, fuel, systems, empennage
% Optimization: Low Drag, Positive Lift generation, low weight, comfort,
% structural flight, external symmetry, low wetted area
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
%% L/D length to diameter ratio
% Table 7.7
% L/D for Concorde is 23.
% shooting for this range
%% Passenger Cabin % Sadraey %
% Calculate Cabin Diameter
W_s = 20/12 ; % seat width ft
% Table 7.3
W_a = 20/12; % aisle distance ft
% Equation 7.2
W_C = 2*W_s + W_a; % Interior Cabin Width ft
% Equation 7.9
D_C = W_C + 2*T_w; % Cabin Diameter ft

% Calculate Cabin Length
P_s = 36.22/12; % seat pitch or distance between seats (ft)
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
L_N = 20; % arbitrary length of nose ft
%% Rear Section
% Will Include Tail and Wing as well for supersonic aircraft.
L_R = 70; % arbitrary length of rear section ft
%% Length of Fuselage
% Equation 7.10 
L_F = L_CP + L_C + L_N + L_R; % Total Fuselage Length
fprintf('The fuselage length is %0.2f ft \n',L_F)
%% Length to Diameter Ratio
LF_DF = L_F/D_C; 
fprintf('The L/D ratio for the fuselage is %0.2f  \n',LF_DF)
%% Fuel Tanks
m_f = 177928.86/9.81; % kg of fuel
% Table 7.8 for fuel types and density
rho_f = 840; % kg/m^3
V_f = m_f/rho_f; % m^3
V_f = V_f * 35.31; % ft^3
fprintf('The fuel tank volume is %0.2f ft^3 \n',V_f)
%% Volume of Bottom Half
% Equation 7.4a
V_bot = 0.5* (pi*W_C^2/4*L_C);
fprintf('The bottom volume is ratio for the fuselage is %0.2f ft^3 \n',V_bot)