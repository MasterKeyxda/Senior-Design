function [V, rho, q ,Cf_Turb,Cf_lam ,Re_Ratio] = dragCalc(h, M, visc, L_F, MAC, ch, cv, nacLength)
%UNTITLED2 Summary of this function goes here
%   % Inputs: 
    % altitude --> h (kft)
    % Mach Number --> M 
    % viscosity --> visc (slug-ft/s)
    % L_F --> fuselage length (ft)
    % MAC --> wing MAC (ft)
    % ch --> ht MAC (ft)
    % cv --> vt MAC (ft)
    % Outputs:
    % Dynamic Pressure --> q (psf)
    % Cf_turb --> array of skin friction coefficients (Order: fuse, wing,
    % ht, vt, nacelle)
    % Cf_lam --> laminar skin friction coeff
    % Re_Ratio --> array of Re ratios for laminar/turb mixed flow
    % V --> cruise speed (ft/s)
    % rho --> density at cruise alt (slugs/ft^3)
[~,~,rho,a] = AltTable(h,'h');

% Cruise Velocity (ft/s)
a = a * 1116.5;
rho = rho * 0.002378;
V = M*a;

% Reynolds Numbers
Re_fuse = (rho * V * L_F) / visc;
Re_wing = (rho * V * MAC) / visc;
Re_ht = (rho * V * ch) / visc;
Re_vt = (rho * V * cv) / visc;
Re_nacelle = (rho * V *  nacLength) / visc;
Re_crit = 1E6; % critical reynolds number

% Reynolds Number Ratios
Re_fuse_Ratio = Re_crit / Re_fuse;
Re_wing_Ratio = Re_crit / Re_wing;
Re_ht_Ratio = Re_crit / Re_ht;
Re_vt_Ratio = Re_crit / Re_vt;
Re_nac_Ratio = Re_crit / Re_nacelle;
Re_Ratio = [Re_fuse_Ratio, Re_wing_Ratio, Re_ht_Ratio, Re_vt_Ratio, Re_nac_Ratio];

% Skin Friction Coefficients 
Cf_lam = 1.327/sqrt(Re_crit); % laminar skin friction coeff

% Assume fully Turbulent Flow
Cf_fuseTurb = 0.455 / ((log10(Re_fuse))^2.58); % Fuselage
Cf_wingTurb = 0.455 / ((log10(Re_wing))^2.58); % Wing
Cf_htTurb = 0.455 / ((log10(Re_ht))^2.58); % Horizontal Tail
Cf_vtTurb = 0.455 / ((log10(Re_vt))^2.58); % Vertical Tail
Cf_nacelleTurb = 0.455 / ((log10(Re_nacelle))^2.58); % Nacelle
Cf_Turb = [Cf_fuseTurb, Cf_wingTurb, Cf_htTurb, Cf_vtTurb, Cf_nacelleTurb];

q = 0.5*rho*(V^2); % dynamic pressure, standard day, cruise altitude
end

