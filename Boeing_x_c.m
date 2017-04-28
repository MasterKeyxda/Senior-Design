%% Boeing HSNLF Airfoil

fname = 'BACNLF';
x_c = 0:0.01:1; 
c_root = 21.78; % root chord (ft)
c_built = 19.87; % chord of built-in section (ft)

% Get BACNLF airfoil thickness distribution at chosen wing station
[tdist_u, tdist_l] = foil_t_get(fname, x_c);

% Determine upper and lower thicknesses at "built-in" section of BACNLF
t_built_u = tdist_u*c_built; % ft
t_built_l = tdist_l*c_built; % ft

% Get total thickness as function of x/c
t_built_dist = t_built_u + t_built_l;

% Spar Locations in terms of chord
spar_loc = [0.20, 0.36, 0.52, 0.68]; 

% Determine airfoil thickness at spar locations 
t_spar = t_built_dist((spar_loc*100)+1); % airfoil thickness in ft
t_spar = t_spar * 12; % convert thickness to inches


%% Biconvex Airfoil

fname = 'biconvex';
c_mid = 17.97; % chord length at start of biconvex airfoil

% Get biconvex airfoil thickness distribution at chosen wing station
[t_dist_u,t_dist_l] = foil_t_get(fname, x_c);

% Determine upper and lower thicknesses at chosen section of biconvex
t_bi_u = t_dist_u*c_mid; % ft
t_bi_l = t_dist_l*c_mid; % ft

% Get total thickness as function of x/c
t_bi_dist = t_bi_u + t_bi_l;

% Determine airfoil thickness at spar locations 
t_spar_bi = t_bi_dist((spar_loc*100)+1); % airfoil thickness in ft
t_spar_bi = t_spar_bi * 12; % convert thickness to inches

% Plot
figure(2)
plot(x_c, t_dist_u)
hold on;
plot(x_c, tdist_l)