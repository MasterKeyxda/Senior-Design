%% Boeing HSNLF Airfoil

fname = 'BACNLF';
x_c = 0:0.01:1; 
c_root = 21.78; % root chord (ft)
c_built = 19.87; % chord of built-in section (ft)

% Get BACNLF airfoil thickness distribution at chosen wing station
t_dist = foil_t_get(fname, x_c);

% Determine max thickness (W.S. = 0)
t_dist_max = max(t_dist);
t_max = t_dist_max*c_root; % max thickness (ft)

% Determine thicknesses at "built-in" section of BACNLF
t_dist_built = t_dist*c_built; % ft

% Spar Locations in terms of chord
spar_loc = [0.20, 0.36, 0.52, 0.68]; 

% Determine airfoil thickness at spar locations 
t_spar = t_dist_built((spar_loc*100)+1); % airfoil thickness in ft
t_spar = t_spar * 12 % convert thickness to inches

%% Biconvex Airfoil

fname = 'biconvex';
c_mid = 17.97; % chord length at start of biconvex airfoil

% Get biconvex airfoil thickness distribution at chosen wing station
t_dist_bi = foil_t_get(fname, x_c);

% Determine thicknesses at chosen section of biconvex
t_dist_biCon = t_dist_bi*c_mid; % ft

