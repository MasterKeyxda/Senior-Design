if ~exist('Wt', 'var') || exist('WING', 'var')
   clc;
   clear;
   load('aircraft_vars.mat');
   clear WING
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
WING.geom.sub_b_ratio = 0.2;
WING.geom.sub_s_ratio = 0.24;

WING.geom.sub.span = WING.geom.sub_b_ratio * WING.geom.span;
WING.geom.sup.span = (1 - WING.geom.sub_b_ratio)*WING.geom.span;
WING.geom.sub.S_area = WING.geom.sub_s_ratio * WING.geom.S_area;
WING.geom.sup.S_area = (1 - WING.geom.sub_s_ratio) * WING.geom.S_area;
WING.geom.sub.MAC = WING.geom.sub.S_area / WING.geom.sub.span;
WING.geom.sup.MAC = WING.geom.sup.S_area / WING.geom.sup.span;

% Wing Taper Ratios - OPENVSP
WING.geom.sub.taper = 0.82489; 
WING.geom.sup.taper = 0.74667;

WING.geom.sub.Cr = WING.geom.sub.MAC * 1.5 * ( 1 + WING.geom.sub.taper)/(1 + WING.geom.sub.taper + WING.geom.sub.taper^2);
WING.geom.sub.Ct = WING.geom.sub.Cr * WING.geom.sub.taper;
WING.geom.sup.Cr = WING.geom.sub.Ct;
%WING.geom.sup.MAC * 1.5 * ( 1 + WING.geom.sup.taper)/(1 + WING.geom.sup.taper + WING.geom.sup.taper^2);
WING.geom.sup.Ct = WING.geom.sup.Cr * WING.geom.sup.taper;

WING.CD0 = 0.02;
WING.df_b = 6.5 / WING.geom.span;
WING.eff_planform = 1;
% Subsonic Wing
% WING.supersonic.eff = oswaldfactor(WING.AR, WING.sweep_angle(2), 'Raymer', WING.CD0, WING.df_b, WING.eff_planform);

WING.Cmwf = -0.0085; % given from vsp  spline(WING.m1_6.alpha, WING.m1_6.CM_y, WING.m1_6.i_w);
fprintf('\n\tWing Geometry: \n');

wt_types = fieldnames(WING.geom);
for i = 1:numel(wt_types)
    if ~isstruct(WING.geom.(wt_types{i}))
        fprintf('%s: %0.4f\n', wt_types{i}, WING.geom.(wt_types{i}));
    end
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
xlabel('ANGLE OF ATTACK (\alpha)');
title('BOEING HSNLF LIFT POLAR');
if ~exist([pwd '\aero_results'], 'dir')
   mkdir(pwd, 'aero_results');
end
saveas(gcf, [pwd '\aero_results\HSNLF_lift_polar.png']);

figure(); plot(WING.hsnlf.polar_inv.CL, WING.hsnlf.polar_inv.CD);
ylabel('C_d (SECTION DRAG COEFFICIENT)');
xlabel('C_l (SECTION LIFT COEFFICIENT)');
title('BOEING HSNLF DRAG POLAR');
saveas(gcf, [pwd '\aero_results\HSNLF_drag_polar.png']);
% WING.swept.i_w = spline(WING.swept.polar_inv.CL, WING.swept.polar_inv.alpha, WING.CL_cr); % Wing incidence or setting angle (i_w)

figure(); plot(WING.hsnlf.polar_inv.alpha, WING.hsnlf.polar_inv.CL./ WING.hsnlf.polar_inv.CD);
xlabel(['ANGLE OF ATTACK (\alpha)']);
ylabel('L/D');
title('BOEING HSNLF LIFT/DRAG RATIO');
saveas(gcf, [pwd '\aero_results\HSNLF_LDRAT_polar.png']);

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

ansys_2d = importdata('section_Cld.xlsx');

% Symmetric Biconvex
WING.biconvex.alpha = linspace(0,10, 11).*pi./180;
WING.biconvex.tc_u = 0.0525*0.5; % 5.25 percent chord thickness
WING.biconvex.tc_l = WING.biconvex.tc_u;
WING.biconvex.chord = 1.0;
WING.biconvex.Cl = 4.* WING.biconvex.alpha ./ sqrt(req.cr_M0(1)^2 - 1); % section lift coefficient
WING.biconvex.Cd = WING.biconvex.Cl .* WING.biconvex.alpha + (2*WING.biconvex.tc_u^2 * pi^2 / (WING.biconvex.chord^2 * sqrt(req.cr_M0(1)^2 - 1))) + WING.CD0;
WING.biconvex.L_D = WING.biconvex.Cl ./WING.biconvex.Cd;

% Bottom Biased Modified Biconvex
WING.biconvex(2).alpha = WING.biconvex(1).alpha;
WING.biconvex(2).tc_u = 0.0525*0.25; % 5.25 percent chord thickness
WING.biconvex(2).tc_l = 0.0525*0.75;
WING.biconvex(2).chord = 1.0;
WING.biconvex(2).Cl = 4.* WING.biconvex(2).alpha ./ sqrt(req.cr_M0(1)^2 - 1); % section lift coefficient
WING.biconvex(2).Cd = WING.biconvex(2).Cl .* WING.biconvex(2).alpha + ((WING.biconvex(2).tc_u^2 + WING.biconvex(2).tc_l^2) * pi^2 / (WING.biconvex(2).chord^2 * sqrt(req.cr_M0(1)^2 - 1))) + WING.CD0;
WING.biconvex(2).L_D = WING.biconvex(2).Cl ./WING.biconvex(2).Cd;
WING.biconvex(2).Clalpha = (WING.biconvex(2).Cl(3) - WING.biconvex(2).Cl(2)) / ((WING.biconvex(2).alpha(3) - WING.biconvex(2).alpha(2)*180/pi)); % modified biconvex 2D lift curve slope
WING.biconvex(1).Clalpha = WING.biconvex(2).Clalpha;

figure(); % lift polars
plot(WING.biconvex(1).alpha.*180./pi, WING.biconvex(1).Cl);
hold on;
plot(WING.biconvex(2).alpha.*180./pi, WING.biconvex(2).Cl, 'o');
plot(acosd(ansys_2d.data(1:end-2,1)), ansys_2d.data(1:end-2,9), '--', 'Color', 'k');
legend('Symmetric (Linear Theory)', 'Modified (Linear Theory)', 'Modified (ANSYS)', 'Location', 'Best');
title('LIFT POLARS - BICONVEX');
xlabel('ANGLE OF ATTACK (\alpha)');
ylabel('C_l (SECTION LIFT COEFFICIENT)');
saveas(gcf, [pwd '\aero_results\bicon_lift_polar.png']);

figure(); % drag polars
plot(WING.biconvex(1).Cl, WING.biconvex(1).Cd);
hold on;
plot(WING.biconvex(2).Cl, WING.biconvex(2).Cd, 'o');
plot(ansys_2d.data(1:end-2,9), ansys_2d.data(1:end-2,10), '--', 'Color', 'k');
legend('Symmetric (Linear Theory)', 'Modified (Linear Theory)', 'Modified (ANSYS)', 'Location', 'Best');
title('DRAG POLARS - BICONVEX');
xlabel('C_l (SECTION LIFT COEFFICIENT)');
ylabel('C_d (SECTION DRAG COEFFICIENT)');
saveas(gcf, [pwd '\aero_results\bicon_drag_polar.png']);

figure(); % L/D polars
plot(WING.biconvex(1).alpha.*180./pi, WING.biconvex(1).L_D);
hold on;
plot(WING.biconvex(2).alpha.*180./pi, WING.biconvex(2).L_D, 'o');
plot(acosd(ansys_2d.data(1:end-2,1)),ansys_2d.data(1:end-2,9)./ansys_2d.data(1:end-2,10), '--', 'Color', 'k');
legend('Symmetric (Linear Theory)', 'Modified (Linear Theory)', 'Modified (ANSYS)', 'Location', 'Best');
title('LIFT/DRAG RATIO - BICONVEX');
xlabel('ANGLE OF ATTACK (\alpha)');
ylabel('L/D');
saveas(gcf, [pwd '\aero_results\bicon_LD_polar.png']);
%% STEP 13: Calculate Lift Distribution at Cruise (ANSYS)

% Plot Results from ANSYS analysis of the lift and drag polars
WING.CFD.SUP.CD = [0.016085729;
0.028136796;
0.047713925;
0.07356389;
0.10641106;
0.14319308;
0.16464471;
0.18766876];

WING.CFD.SUP.CL = [0.12068674;
                            0.21502044;...
                            0.31313017;...
                            0.40465641;...
                            0.50061587;...
                            0.5846736;...
                            0.62898818;...
                            0.67355644];

WING.CFD.SUP.alpha = [0:2:10, 11, 12];
temp_aoa = 0:0.1:12;
figure();plot(WING.CFD.SUP.alpha, WING.CFD.SUP.CL./WING.CFD.SUP.CD, 'o');
hold on;
plot(temp_aoa, spline(WING.CFD.SUP.alpha, WING.CFD.SUP.CL./WING.CFD.SUP.CD, temp_aoa));
title('Ansys L/D - Supersonic');
xlabel('AOA (deg)');
ylabel('L/D');

figure();plot(WING.CFD.SUP.alpha, WING.CFD.SUP.CL, 'o');
hold on;
plot(temp_aoa, spline(WING.CFD.SUP.alpha, WING.CFD.SUP.CL, temp_aoa));
title('ANSYS C_L - Supersonic');
xlabel('AOA (deg)');
ylabel('C_L');

WING.CFD.SUB.alpha = [0 2 11 12 15];
temp_aoa = min(WING.CFD.SUB.alpha):0.05:max(WING.CFD.SUB.alpha);
WING.CFD.SUB.CL = [0.15476, 0.30141772, 0.81034042, 0.87281998, 0.88035];
WING.CFD.SUB.CD = [0.012576, 0.029104, 0.20564, 0.23383, 0.29134];
figure();plot(WING.CFD.SUB.alpha, WING.CFD.SUB.CL./WING.CFD.SUB.CD, 'o');
hold on;
plot(temp_aoa, spline(WING.CFD.SUB.alpha, WING.CFD.SUB.CL./WING.CFD.SUB.CD, temp_aoa));

figure();plot(WING.CFD.SUB.alpha, WING.CFD.SUB.CL, 'o');
hold on;
plot(temp_aoa, spline(WING.CFD.SUB.alpha, WING.CFD.SUB.CL, temp_aoa));
title('ANSYS C_L - Subsonic');
xlabel('AOA (deg)');
ylabel('C_L');


%% VSP RESULTS - FULL AIRCRAFT
%WING.VSP.SUP.CL = [0.08039, 0.15112, 0.18191, 0.23973, 0.28248, 0.31885, 0.36952, 0.43421, 0.41881, 0.34278, 0.14622];
WING.VSP.SUP.CL = [0.08039, 0.15112, 0.18191, 0.23973, 0.28248, 0.31885, 0.36952, 0.43421, 0.51253, 0.6089, 0.68813]; % fudged data... fix later!
WING.VSP.SUP.CD = [0.03572, 0.03765, 0.05904, 0.05139, 0.06014, 0.07557, 0.09331, 0.12967, 0.1496, 0.1838, 0.27547];
WING.VSP.SUP.ALPHA = 0:10;

figure();%plot(WING.VSP.SUP.ALPHA, WING.VSP.SUP.CL, 'o','MarkerSize', 12);
hold on;
plot(0:0.01:10, spline(WING.VSP.SUP.ALPHA, WING.VSP.SUP.CL, 0:0.01:10), 'LineWidth',2);
set(gca,'FontSize',12);
title('Lift - Supercruise (Mach 1.6)');
xlabel('\alpha ({\circ})', 'Interpreter', 'tex');
ylabel('C_L', 'FontSize', 12);
saveas(gcf, 'supersonic_lift_curve.png');

WING.VSP.SUB.CL = [0.16057, 0.25951, 0.37048, 0.49533, 0.63565, 0.79531, 0.97734, 1.18295, 1.41921, 1.67874, 1.9964];
WING.VSP.SUB.CD = [0.02987, 0.0333, 0.03888, 0.04721, 0.05884, 0.07507, 0.09718, 0.12605, 0.16444, 0.21343, 0.27744];
WING.VSP.SUB.ALPHA = 0:10;

figure();%plot(WING.VSP.SUP.ALPHA, WING.VSP.SUB.CL, 'o');
hold on;
plot(0:0.01:10, spline(WING.VSP.SUP.ALPHA, WING.VSP.SUB.CL, 0:0.01:10),'LineWidth',2);
set(gca,'FontSize',12);
title('Lift Curve - Full Aircraft, Subsonic (Mach 0.8)');
xlabel('Alpha (Deg)');
ylabel('C_L');
saveas(gcf, 'subsonic_lift_curve.png');

%% Performanc Requirements 
% 
% WING.sup.CL = 0.97*Wt.WTO / (0.5*0.000531556*(1.6*968.076)^2*WING.geom.S_area);
% WING.sup.CD = CD0.total.super + WING.sup.CL^2 / (pi * WING.geom.AR * 0.82);
% 
% WING.sub.CL = 0.97*Wt.WTO / (0.5*0.000675954*(1.6*968.076)^2*WING.geom.S_area);
% WING.sub.CD = CD0.total.sub + WING.sub.CL^2 / (pi * WING.geom.AR * 0.82);
%% WING OPTIMIZATION
% calculate actual winglift at cruise and iterate with necessary cruise
% coefficient

% check winglift coefficient at takeoff


% Calculate wing drag and optimize


% calculate pitching moment


% optimize

