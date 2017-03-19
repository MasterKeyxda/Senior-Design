function [Tts4, f_ratio_T, SFC_T] = throttle_calcs(FC, ENG, TL, thrust)

% FC: flight conditions (P0, T0, M0, gamma, cp)
% ENG: engine parameters (pi_c, pi_f, Tt4_max, hPR, ENG.p0_9, ENG.p0_19, d (diameter, ft))
% THR: required thrust
% TL: tech level

% Get Speed of Sound
Rc = (FC.gamma - 1)/(FC.gamma) * FC.cp * 778.16; % BTU/lbm *R -> lbf. ft/ lbm * R
gc = 32.174; % Gravity Constant
a0 = sqrt(FC.gamma * Rc * gc * FC.T0);

area = pi * 0.25 * (ENG.d)^2;

% Ram
TR.r = 1 + 0.5*(FC.gamma - 1)*FC.M0^2;
PR.r = TR.r ^ (FC.gamma/ (FC.gamma - 1));

mdot = area * FC.M0 * a0 * FC.rho * gc;
if FC.M0 == 0.0
    mdot = ENG.mdot_sl; % lb/s
end
Fsp_req = thrust/mdot; % get required specific thrust
Tt4_vals = linspace(ENG.Tt4_min, ENG.Tt4_max, 1000);

% Diffuser
TR.d = 1;
if FC.M0<=1
    eta_r = 1;
else
    eta_r = 1 - 0.075*(FC.M0 - 1)^1.35;
end

PR.d = TL.d_max * eta_r;

% Compressor
TR.cR = ENG.pi_c ^ ((FC.gamma - 1)/(FC.gamma * TL.e_c)); % reference tau_c
TR.c = 1 + (TR.cR - 1) .* (Tt4_vals./(FC.T0 * TR.r))./(ENG.Tt4_sl / ENG.Tt0);
PR.c = TR.c.^ (FC.gamma * TL.e_c ./ (FC.gamma - 1));

Tt.s3 = FC.T0 * TR.r * TR.d * TR.c;
Pt.s3 = FC.P0 * PR.r * PR.d * PR.c;

% Fan Conditions
TR.f = ENG.pi_f ^ ((FC.gamma - 1)/(FC.gamma * TL.e_f));

% Turbine Conditions
%<<<<<<< Updated upstream
%=======
Tt4_vals = linspace(2200, ENG.Tt4_max, 1000);
%>>>>>>> Stashed changes
cp_t = [];
gamma_t = [];
for i = 1:length(Tt4_vals)
    [cp_t(i), gamma_t(i)] = cp_val(Tt4_vals(i), 'combustion');
end
TR.lambda = cp_t.*Tt4_vals ./ (FC.cp * FC.T0);

f_ratio = (cp_t.*Tt4_vals - (FC.cp * Tt.s3))./(TL.eta_b * ENG.hPR - cp_t.*Tt4_vals);

TR.t = 1 - (TR.r ./ TR.lambda) .* (1./(TL.eta_m .* (1 + f_ratio))) .* (TR.c - 1 + ENG.BPR * (TR.f - 1));
PR.t = TR.t .^(gamma_t ./ ((gamma_t - 1).* TL.e_t));

% Nozzle Conditions
PR.s9 = ENG.p0_9 * PR.r * PR.d * ENG.pi_c * TL.pi_b * PR.t * TL.pi_n; % Pt9/P9
TR.s9_0 = (Tt4_vals .* TR.t / FC.T0)./(PR.s9.^((gamma_t - 1)./gamma_t)); % T9/T0

M9 = sqrt(2./(gamma_t - 1) .* ((PR.s9.^((gamma_t - 1)./gamma_t)) - 1));
Rt = (gamma_t - 1)./gamma_t .* cp_t * 778.16; % BTU/lbm*R -> lbf ft/ lbm * R
V9_a0 = M9.*sqrt((gamma_t .* Rt) ./ (FC.gamma * Rc) .* TR.s9_0);

% Fan Nozzle Conditions
PR.s19 = ENG.p0_19 * PR.r * PR.d * ENG.pi_f * TL.pi_n; % Pt19/P19
TR.s19_0 = (TR.r * TR.f) ./ ( PR.s19 .^ ((FC.gamma- 1)./FC.gamma));
M19 = sqrt(2./(FC.gamma - 1) .* (PR.s19 .^ ((FC.gamma- 1)./FC.gamma) - 1));
V19_a0 = M19.*sqrt(TR.s19_0);

% Thrust and SFC Performance
Fsp = 1./(1 + ENG.BPR) .* (a0 / gc) .* ((1 + f_ratio) .* V9_a0 - FC.M0 + (1 + f_ratio).*(Rt.* TR.s9_0 ./ (Rc .* V9_a0) .* (1 - ENG.p0_9)./FC.gamma)) + ENG.BPR/(1 + ENG.BPR) .* (a0 / gc) .*(V19_a0 - FC.M0 +(TR.s19_0 ./ (V19_a0) .* (1 - ENG.p0_19)./FC.gamma));
SFC = f_ratio ./ ((1 + ENG.BPR) .* Fsp);

Tts4 = spline(Fsp, Tt4_vals, Fsp_req);
f_ratio_T = spline(Fsp, f_ratio, Fsp_req);
SFC_T = spline(Fsp, SFC, Fsp_req);


    function[cp, gamma] = cp_val(temp, type)
        if strcmp(type, 'air')
            if temp <= 1000
                cp = 0.24;
                gamma = 1.4;
            elseif temp>1000
                cp = 0.258;
                gamma = 1.36;
            end
        elseif strcmp(type, 'combustion')
            if temp < 1300
                cp = 0.258;
                gamma = 1.36;
            elseif (temp > 1300) && (temp<1800)
                cp = 0.276;
                gamma = 1.33;
            else
                cp = 0.295;
                gamma = 1.3;
            end
        end
    end


end