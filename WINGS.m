%% Wing Calculations - Preliminary
fprintf('\nWing Prelim Design\n');

% Define Sweep Angle
WING.M0_angle = asin(1./req.cr_M0)*180/pi;
WING.sweep_angle = 1.2*(90-WING.M0_angle);
fprintf('Mach Angles:\nMin: %0.3f\tMax: %0.3f\n', WING.M0_angle(1), WING.M0_angle(2));
fprintf('Sweep Angles:\nMin: %0.3f\tMax: %0.3f\n', WING.sweep_angle(1), WING.sweep_angle(2));

% Define Area 
WING.S_area = 1166; % ft^2

% Number of Wings
WING.no = 1;
WING.vert_loc = 'low';
WING.config = 'swept_tapered';

WING.wt_avg = 0.5*WTO*(1 + 1- Wt.fuel.Wf_Wto);
[~,~,sig_rho,~] = AltTable(atm.alt, 'h');
WING.CL_cr = 2*WING.wt_avg/(sig_rho * 0.002378 * Wt.fuel.V_max_cr^2 * WING.S_area);

WING.CL_max = 2.4; 
WING.V_TO = 1.2*sqrt(2*WTO / (0.002378 * WING.CL_max * WING.S_area)); % ft/s
WING.CL_TO = 0.85 * 2 * WTO / (0.002378 * WING.V_TO^2 * WING.S_area);

% Heavy Lift Devices
WING.dCL_flaps = [1 1.3]; % fowler flaps -> Sadraey

% Determine airfoil design


% Select the dihedral angle


% Select AR, taper ratio, twist angle


% Get lift distribution and determine elliptic behavior


% calculate actual winglift at cruise and iterate with necessary cruise
% coefficient


% check winglift coefficient at takeoff


% Calculate wing drag and optimize


% calculate pitching moment


% optimize

% Taper Ratio
