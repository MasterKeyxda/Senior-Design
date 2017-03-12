if ~exist('Wt', 'var')
   clc;
   clear;
   load('aircraft_vars.mat'); 
   close all;
end


%% Wing Calculations - Preliminary
% fprintf('\nWing Prelim Design\n');
addpath([pwd '/aero_tools/']);
addpath([pwd '/ANSYS_results/']);
addpath([pwd '/VSP_results/']);

%% PRELIMINARY RESULTS OF WING DESIGN AND OPTIMIZATION

% Define Initial Wing Characteristics
WING.geom.S_area = S_w; % ft^2
WING.geom.AR = AR;

% Determine Span & MAC
WING.geom.span = sqrt(WING.geom.S_area * WING.geom.AR);
WING.geom.MAC = sqrt(WING.geom.S_area/WING.geom.AR);
WING.geom.sub_span = 0.2*WING.geom.span;
WING.geom.sup_span = 0.8 * WING.geom.span;


% Wing Taper Ratios - OPENVSP
WING.geom.taper_sub = 0.82489; 
WING.geom.taper_super = 0.74667;

WING.CD0 = 0.02;
WING.df_b = 6.5 / WING.geom.span;
WING.eff_planform = 1;
% Subsonic Wing
% WING.supersonic.eff = oswaldfactor(WING.AR, WING.sweep_angle(2), 'Raymer', WING.CD0, WING.df_b, WING.eff_planform);

WING.Cmwf = -0.0085; % given from vsp  spline(WING.m1_6.alpha, WING.m1_6.CM_y, WING.m1_6.i_w);
fprintf('\n\tWing Geometry: \n');

wt_types = fieldnames(WING.geom);
for i = 1:numel(wt_types)
    fprintf('%s: %0.4f\n', wt_types{i}, WING.geom.(wt_types{i}));
end
clear wt_types

%% STEP 1-3: Number of Wings & Vertical Location & Configuration
WING.no = 1;
WING.vert_loc = 'low';
WING.config = 'swept_tapered';

%% STEP 4-6: Calculate average cruise weight, required cruise CL, required CL_TO

% Average Cruise Weight
WING.wt_avg = 0.5*Wt.WTO*(1 + 1- Wt.fuel.Wf_Wto);
[~,~,sig_rho,~] = AltTable(atm.alt, 'h');

% Cruise CL
WING.req.CL_cr = 2*WING.wt_avg/(sig_rho * 0.002378 * Wt.fuel.V_max_cr^2 * WING.geom.S_area);

% Takeoff CL
WING.req.CL_max = 2.4; 
WING.req.V_TO = 1.2*sqrt(2*Wt.WTO / (0.002378 * WING.req.CL_max * WING.geom.S_area)); % ft/s
WING.req.CL_TO = 0.85 * 2 * Wt.WTO / (0.002378 * WING.req.V_TO^2 * WING.geom.S_area);

fprintf('\n\tWing Required Performance: \n');
wt_types = fieldnames(WING.req);
for i = 1:numel(wt_types)
    fprintf('%s: %0.4f\n', wt_types{i}, WING.req.(wt_types{i}));
end
clear wt_types

%% STEP 7-8: Select High-Lift Device and Geometry

% Heavy Lift Devices
% WING.dCL_flaps = [1 1.3]; % fowler flaps -> Sadraey

%% STEP 11: Determine Subsonic Taper Angles and Dihedral Angle

WING.dihedral = 6; % degrees

%% STEP 12: Select Aspect Ratio, Taper Ratio, and Wing Twist Angle

%% STEP 9-10: Select Airfoil and determine wing incidence angle

% Subsonic HSNLF Characteristics
WING.hsnlf.name = 'BACNLF';
if ~exist('boeing_polars.mat', 'file')
    [swept.polar_inv, swept.foil_inv] = xfoil([WING.swept.name '.dat'], [-2:19], 10e6, 0.0, 'oper/iter 10000');
    WING.hsnlf.polar_inv = swept.polar_inv;
    WING.hsnlf.foil_inv = swept.foil_inv;
else
    load('boeing_polars.mat');
    WING.hsnlf.polar_inv = swept.polar_inv;
    WING.hsnlf.foil_inv = swept.foil_inv;
    clear swept;
end

figure(); plot(WING.hsnlf.polar_inv.alpha, WING.hsnlf.polar_inv.CL);
ylabel('C_l');
xlabel('\alpha');
title('Boeing HSNLF Lift-Curve Polar');
% WING.swept.i_w = spline(WING.swept.polar_inv.CL, WING.swept.polar_inv.alpha, WING.CL_cr); % Wing incidence or setting angle (i_w)

% Supersonic Airfoil Characteristics
% WING.supersonic.name = 'biconvex';
% WING.supersonic.t_c = 0.12; % thickness ratio
% cd_data = dlmread('cd-history', ' ', 2, 0);
% cl_data = dlmread('cl-history', ' ', 2, 0);
% cm_data = dlmread('cm-history', ' ', 2, 0);
% iter = find(cd_data(:,1) == 1);
% alpha = 0:1:length(iter)-1;
% i = 1;
% ii = 1;
% while i <= length(iter)
%     if i < length(iter)
%         if cd_data(iter(i+1)-1,1) ~= 5000
%             WING.supersonic.Cd(ii) = cd_data(iter(i+1)-1,end);
%             WING.supersonic.Cl(ii) = cl_data(iter(i+1)-1,end);
%             WING.supersonic.Cm(ii) = cm_data(iter(i+1)-1,end);
%             WING.supersonic.alpha(ii) = alpha(i);
%             ii = ii + 1;
%         end
%     else
%         if cd_data(end,1) ~= 5000
%             WING.supersonic.Cd(ii) = cd_data(end,end);
%             WING.supersonic.Cl(ii) = cl_data(end,end);
%             WING.supersonic.Cm(ii) = cm_data(end,end);
%             WING.supersonic.alpha(ii) = alpha(i);
%         end
%     end
%     i = i + 1;
% end
% 
% figure(); plot(WING.supersonic.alpha(WING.supersonic.Cl>=0), WING.supersonic.Cl(WING.supersonic.Cl>=0),'o-');
% 
% ylabel('C_l');
% xlabel('\alpha');
% title('Bi-convex Lift-Curve Polar');
% 
% figure(); plot(WING.supersonic.Cl(WING.supersonic.Cl>=0), WING.supersonic.Cd(WING.supersonic.Cl>=0),'o-');
% ylabel('C_d');
% xlabel('C_l');
% title('Bi-convex Drag Polar');

% Flat Plate Formulas in Supersonic flow

% Symmetric Biconvex
WING.biconvex.alpha = (0:20).*pi ./ 180;
WING.biconvex.tc_u = 0.03*0.5; % 5.25 percent chord thickness
WING.biconvex.tc_l = WING.biconvex.tc_u;
WING.biconvex.chord = 1.0;
WING.biconvex.Cl = 4.* WING.biconvex.alpha ./ sqrt(req.cr_M0(1)^2 - 1); % section lift coefficient
WING.biconvex.Cd = WING.biconvex.Cl .* WING.biconvex.alpha + (2*WING.biconvex.tc_u^2 * pi^2 / (WING.biconvex.chord^2 * sqrt(req.cr_M0(1)^2 - 1))) + WING.CD0;
WING.biconvex.L_D = WING.biconvex.Cl ./WING.biconvex.Cd;

% Bottom Biased Modified Biconvex
WING.biconvex(2).alpha = WING.biconvex(1).alpha;
WING.biconvex(2).tc_u = 0.03*0.25; % 5.25 percent chord thickness
WING.biconvex(2).tc_l = 0.03*0.75;
WING.biconvex(2).chord = 1.0;
WING.biconvex(2).Cl = 4.* WING.biconvex(2).alpha ./ sqrt(req.cr_M0(1)^2 - 1); % section lift coefficient
WING.biconvex(2).Cd = WING.biconvex(2).Cl .* WING.biconvex(2).alpha + ((WING.biconvex(2).tc_u^2 + WING.biconvex(2).tc_l^2) * pi^2 / (WING.biconvex(2).chord^2 * sqrt(req.cr_M0(1)^2 - 1))) + WING.CD0;
WING.biconvex(2).L_D = WING.biconvex(2).Cl ./WING.biconvex(2).Cd;

figure(); % lift polars
plot(WING.biconvex(1).alpha, WING.biconvex(1).Cl);
hold on;
plot(WING.biconvex(2).alpha, WING.biconvex(2).Cl, 'o');
legend('Symmetric', 'Modified');
title('Lift Polars');

figure(); % drag polars
plot(WING.biconvex(1).Cl, WING.biconvex(1).Cd);
hold on;
plot(WING.biconvex(2).Cl, WING.biconvex(2).Cd, 'o');
legend('Symmetric', 'Modified');
title('Drag Polars');

figure(); % L/D polars
plot(WING.biconvex(1).alpha, WING.biconvex(1).L_D);
hold on;
plot(WING.biconvex(2).alpha, WING.biconvex(2).L_D, 'o');
legend('Symmetric', 'Modified');
title('Lift over Drag');
%% STEP 13: Calculate Lift Distribution at Cruise (Lifting Line Theory?)

% Plot Results from ANSYS analysis of the lift and drag polars

%% WING OPTIMIZATION
% calculate actual winglift at cruise and iterate with necessary cruise
% coefficient

% check winglift coefficient at takeoff


% Calculate wing drag and optimize


% calculate pitching moment


% optimize

