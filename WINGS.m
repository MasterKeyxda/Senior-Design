if ~exist('Wt', 'var')
   clc;
   clear;
   load('aircraft_vars.mat'); 
   close all;
end


%% Wing Calculations - Preliminary
fprintf('\nWing Prelim Design\n');
addpath([pwd '/aero_tools/']);
addpath([pwd '/ANSYS_results/']);
addpath([pwd '/VSP_results/']);

% Define Initial Wing Characteristics
WING.S_area = S_w; % ft^2
WING.AR = AR;

% Changes as of 2/25/2017
    % Airfoil is now 5.25% of chord
    % Aspect Ratio is now 3
    % Inner sweep modified to 35 degrees

%% STEP 1-3: Number of Wings & Vertical Location & Configuration
WING.no = 1;
WING.vert_loc = 'low';
WING.config = 'swept_tapered';

%% STEP 4-6: Calculate average cruise weight, required cruise CL, required CL_TO

% Average Cruise Weight
WING.wt_avg = 0.5*WTO*(1 + 1- Wt.fuel.Wf_Wto);
[~,~,sig_rho,~] = AltTable(atm.alt, 'h');

% Cruise CL
WING.CL_cr = 2*WING.wt_avg/(sig_rho * 0.002378 * Wt.fuel.V_max_cr^2 * WING.S_area);

% Takeoff CL
WING.CL_max = 2.4; 
WING.V_TO = 1.2*sqrt(2*WTO / (0.002378 * WING.CL_max * WING.S_area)); % ft/s
WING.CL_TO = 0.85 * 2 * WTO / (0.002378 * WING.V_TO^2 * WING.S_area);

%% STEP 7-8: Select High-Lift Device and Geometry

% Heavy Lift Devices
WING.dCL_flaps = [1 1.3]; % fowler flaps -> Sadraey

%% STEP 11: Determine Subsonic Taper Angles and Dihedral Angle
nose_angle = 13.24 * pi / 180; % degrees -> radians
WING.M0_angle = pi/180*((60.54-54.89)/2 * (nose_angle*180/pi - 12) + 54.89); % radians
WING.M_LE = 0.8;
WING.M_cone = (1./(sin(WING.M0_angle - nose_angle))).*sqrt((2/(1.4-1)+(req.cr_M0.*sin(WING.M0_angle)).^2)./((2*1.4/(1.4-1)).*(req.cr_M0.*sin(WING.M0_angle)).^2 - 1.0));
WING.sweep_LE = acosd(WING.M_LE./WING.M_cone);
fprintf('Mach Angles: %0.3f\n', WING.M0_angle);
fprintf('Sweep Angles: %0.3f\n', WING.sweep_LE(1));

% Select the dihedral angle
WING.dihedral = 6; % degrees

%% STEP 12: Select Aspect Ratio, Taper Ratio, and Wing Twist Angle

% Aspect Ratio - based on modified Oswald span efficiency for swept wing
% (Sadraey P. 237) [defunct]
% WING.swept.eff = 0.9;
% WING.swept.AR = ((((3.1+WING.swept.eff)/(4.61*cosd(WING.sweep_angle(2))^0.15))-1)/-0.045)^(1/0.68);
% disp(WING.swept.AR);

% 2/15/2017 - efficiency based on Raymer's approximations with pre-defined aspect
% ratio

% Determine Span & MAC
WING.span = sqrt(S_w * WING.AR);
WING.MAC = sqrt(S_w/WING.AR);

% Wing Taper Ratios - Raymer's Approximations
WING.subsonic.taper = 0.3; 
WING.supersonic.taper = 0.45;

WING.CD0 = 0.02;
WING.df_b = 6.5 / WING.span;
WING.eff_planform = 1;
% Subsonic Wing
% WING.supersonic.eff = oswaldfactor(WING.AR, WING.sweep_angle(2), 'Raymer', WING.CD0, WING.df_b, WING.eff_planform);

fprintf('\n\tWing Geometry: \n');

fprintf('Aspect Ratio: %0.2f\n', WING.AR);
fprintf('Span: %0.3f\n', WING.span);
fprintf('MAC: %0.3f\n', WING.MAC);
fprintf('Wing Area: %0.5f\n', WING.S_area);


%% STEP 9-10: Select Airfoil and determine wing incidence angle

% Subsonic HSNLF Characteristics
WING.swept.name = 'BACNLF';
if ~exist('boeing_polars.mat', 'file')
    [swept.polar_inv, swept.foil_inv] = xfoil([WING.swept.name '.dat'], [-2:19], 10e6, 0.0, 'oper/iter 10000');
    WING.swept.polar_inv = swept.polar_inv;
    WING.swept.foil_inv = swept.foil_inv;
else
    load('boeing_polars.mat');
    WING.swept.polar_inv = swept.polar_inv;
    WING.swept.foil_inv = swept.foil_inv;
    clear swept;
end

figure(); plot(WING.swept.polar_inv.alpha, WING.swept.polar_inv.CL);
ylabel('C_l');
xlabel('\alpha');
title('Boeing HSNLF Lift-Curve Polar');
% WING.swept.i_w = spline(WING.swept.polar_inv.CL, WING.swept.polar_inv.alpha, WING.CL_cr); % Wing incidence or setting angle (i_w)

% Supersonic Airfoil Characteristics
WING.supersonic.name = 'biconvex';
WING.supersonic.t_c = 0.12; % thickness ratio
cd_data = dlmread('cd-history', ' ', 2, 0);
cl_data = dlmread('cl-history', ' ', 2, 0);
cm_data = dlmread('cm-history', ' ', 2, 0);
iter = find(cd_data(:,1) == 1);
alpha = 0:1:length(iter)-1;
i = 1;
ii = 1;
while i <= length(iter)
    if i < length(iter)
        if cd_data(iter(i+1)-1,1) ~= 5000
            WING.supersonic.Cd(ii) = cd_data(iter(i+1)-1,end);
            WING.supersonic.Cl(ii) = cl_data(iter(i+1)-1,end);
            WING.supersonic.Cm(ii) = cm_data(iter(i+1)-1,end);
            WING.supersonic.alpha(ii) = alpha(i);
            ii = ii + 1;
        end
    else
        if cd_data(end,1) ~= 5000
            WING.supersonic.Cd(ii) = cd_data(end,end);
            WING.supersonic.Cl(ii) = cl_data(end,end);
            WING.supersonic.Cm(ii) = cm_data(end,end);
            WING.supersonic.alpha(ii) = alpha(i);
        end
    end
    i = i + 1;
end

figure(); plot(WING.supersonic.alpha(WING.supersonic.Cl>=0), WING.supersonic.Cl(WING.supersonic.Cl>=0),'o-');

ylabel('C_l');
xlabel('\alpha');
title('Bi-convex Lift-Curve Polar');

figure(); plot(WING.supersonic.Cl(WING.supersonic.Cl>=0), WING.supersonic.Cd(WING.supersonic.Cl>=0),'o-');
ylabel('C_d');
xlabel('C_l');
title('Bi-convex Drag Polar');

%% STEP 13: Calculate Lift Distribution at Cruise (Lifting Line Theory?)

WING.OpenVSP_init = [1.0, 0.03444, -0.03805; 2.0, 0.06732, -0.06808]; % used OpenVSP to calculate CL and CM numbers
WING.x_AC = -((WING.OpenVSP_init(2,3) - WING.OpenVSP_init(1,3))/(WING.OpenVSP_init(2,2) - WING.OpenVSP_init(1,2)))/WING.MAC;
fprintf('Aerodynamic Center: %0.5f\n', WING.x_AC);

% Get lift distribution and determine elliptic behavior - get CL, CD, and
% CM
WING.m1_6.alpha = -5:19;
cells = linspace(0,(length(WING.m1_6.alpha)-1)*78, length(WING.m1_6.alpha));
[~, ~, raw1_6] = xlsread('mach_1.6.csv', 'mach_1.6', sprintf('F%i:F%i',1,length(WING.m1_6.alpha)*78));

WING.m1_8.alpha = WING.m1_6.alpha;

% Mach 1.6
WING.m1_6.CL = [raw1_6{cells+14}];
WING.m1_6.CD_tot = [raw1_6{cells+10}];
WING.m1_6.CM_y = [raw1_6{cells+16}];
WING.m1_6.E = [raw1_6{cells+19}];
WING.m1_6.i_w = spline(WING.m1_6.CL, WING.m1_6.alpha, WING.CL_cr); % Wing incidence or setting angle (i_w)

WING.Cmwf = -0.0085; % given from vsp  spline(WING.m1_6.alpha, WING.m1_6.CM_y, WING.m1_6.i_w);

fprintf('Incidence Angles: %0.5f (M 1.6)\n', WING.m1_6.i_w);

% Plot CL's
figure();
plot(WING.m1_6.alpha, WING.m1_6.CL);
hold on;
% plot(WING.m1_8.alpha, WING.m1_8.CL);
xlabel('\alpha');
ylabel('C_L');
title('Wing Lift Coefficient');
legend('Mach 1.6');

% Plot CD's
figure();
plot(WING.m1_6.alpha, WING.m1_6.CD_tot);
hold on;
% plot(WING.m1_8.alpha, WING.m1_8.CD_tot);
xlabel('\alpha');
ylabel('C_{D_{tot}}');
title('Drag Coefficient');
legend('Mach 1.6');

% % Plot CM_y
% figure();
% plot(WING.m1_6.alpha, WING.m1_6.CM_y);
% hold on;
% plot(WING.m1_8.alpha, WING.m1_8.CM_y);
% xlabel('\alpha');
% ylabel('C_{M_{y}}');
% title('Pitch Moment Coefficient');
% legend('Mach 1.6', 'Mach 1.8');
% 
% % Plot E
% figure();
% plot(WING.m1_6.alpha, WING.m1_6.E);
% hold on;
% plot(WING.m1_8.alpha, WING.m1_8.E);
% xlabel('\alpha');
% ylabel('E');
% title('Oswald Efficiency');
% legend('Mach 1.6', 'Mach 1.8');

%% WING OPTIMIZATION
% calculate actual winglift at cruise and iterate with necessary cruise
% coefficient

% check winglift coefficient at takeoff


% Calculate wing drag and optimize


% calculate pitching moment


% optimize
