%% Wing Calculations - Preliminary
fprintf('\nWing Prelim Design\n');
addpath([pwd '\aero_tools\']);
% Define Area 
WING.S_area = 1166; % ft^2

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

%% STEP 9-10: Select Airfoil and determine wing incidence angle

% Subsonic HSNLF Characteristics
WING.swept.name = 'BACNLF';
[WING.swept.polar_inv, WING.swept.foil_inv] = xfoil([WING.swept.name '.dat'], [-2:20], 10e6, 0.0, 'oper/iter 10000');
figure(); plot(WING.swept.polar_inv.alpha, WING.swept.polar_inv.CL);
ylabel('C_l');
xlabel('\alpha');
% WING.swept.i_w = spline(WING.swept.polar_inv.CL, WING.swept.polar_inv.alpha, WING.CL_cr); % Wing incidence or setting angle (i_w)

% Supersonic Airfoil Characteristics
WING.supersonic.name = 'biconvex';
WING.supersonic.t_c = 0.12; % thickness ratio
dx = 0.01;
WING.supersonic.coords = [[1:-dx:0, dx:dx:1]', 2.*WING.supersonic.t_c.*[(1-(1:-dx:0)).*(1:-dx:0),-(1-(dx:dx:1)).*(dx:dx:1)]'];
figure();
plot(WING.supersonic.coords(:,1), WING.supersonic.coords(:,2));
close(gcf);
if exist(sprintf('%s.dat', WING.supersonic.name)) ~= 2 % generate .dat file if not found
    fid = fopen(sprintf('%s.dat', WING.supersonic.name), 'w');
    fprintf(fid, '%s Supersonic Airfoil\n\n', WING.supersonic.name);
    fprintf(fid, '%12.7f  %12.7f\n', WING.supersonic.coords);
    fclose(fid);
end
[WING.supersonic.polar_inv, WING.supersonic.foil] = xfoil(WING.supersonic.coords, [-2:20], 10e6, req.cr_M0(2), 'oper/iter 1000');
figure(); plot(WING.supersonic.polar_inv.alpha, WING.supersonic.polar_inv.CL);
ylabel('C_l');
xlabel('\alpha');
WING.supersonic.i_w = spline(WING.supersonic.polar_inv.CL, WING.supersonic.polar_inv.alpha, WING.CL_cr); % Wing incidence or setting angle (i_w)


%% STEP 11: Determine Subsonic Taper Angles and Dihedral Angle
WING.M0_angle = asin(1./req.cr_M0)*180/pi;
nose_angle = 13.24; % degrees
WING.sweep_angle = 1.2*(90-WING.M0_angle);
fprintf('Mach Angles:\nMin: %0.3f\tMax: %0.3f\n', WING.M0_angle(1), WING.M0_angle(2));
fprintf('Sweep Angles:\nMin: %0.3f\tMax: %0.3f\n', WING.sweep_angle(1), WING.sweep_angle(2));

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
WING.AR = 4;

% Determine Span & MAC
WING.span = sqrt(S_w * WING.AR);
WING.MAC = sqrt(S_w/WING.AR);

% Wing Taper Ratios - Raymer's Approximations
WING.subsonic.taper = 0.3; 
WING.supersonic.taper = 0.45;

WING.CD0 = 0.02;
WING.df_b = 0;
WING.eff_planfom = 1;
% Subsonic Wing
WING.supersonic.eff = oswaldfactor(WING.AR, WING.sweep_angle(2), 'Raymer', WING.CD0, WING.df_b, WING.eff_planform);

% Supersonic wing
% WING.supersonic.eff = 1.78*(1-0.045*WING.supersonic.AR^0.68)-0.64;
WING.supersonic.eff = oswaldfactor(WING.AR, 0, 'Raymer', WING.CD0, WING.df_b, WING.eff_planform);
disp(WING.supersonic.eff);

%% STEP 13: Calculate Lift Distribution at Cruise (Lifting Line Theory?)
% Get lift distribution and determine elliptic behavior


% calculate actual winglift at cruise and iterate with necessary cruise
% coefficient


% check winglift coefficient at takeoff


% Calculate wing drag and optimize


% calculate pitching moment


% optimize
