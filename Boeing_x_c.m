%% Boeing HSNLF Airfoil
fname = 'BACNLF';
x_c = 0:0.01:1; 
c_root = 21.78; % root chord (ft)
c_built = 19.87; % chord of built-in section (ft)

% Get BACNLF airfoil thickness distribution at chosen wing station
t_dist = foil_t_get(fname, x_c);

% Determine max thickness (W.S. = 0)
tdist_max = max(t_dist);
t_max = tdist_max*c_root; % max thickness (ft)

% Determine thicknesses at "built-in" section of BACNLF
tdist_built = t_dist*c_built; % max thickness (ft)

% Spar Locations in terms of chord
spar_loc = [0.20, 0.36, 0.52, 0.68]; 

% Determine airfoil thickness at spar locations 
t_spar = tdist_built((spar_loc*100)+1); % airfoil thickness in ft
t_spar = t_spar * 12 % convert thickness to inches

%% Biconvex Airfoil
