function [Tts4, f_ratio_T, SFC_T] = throttle_calcs(FC, ENG, TL, thrust)

% FC: flight conditions (P0, T0, M0, gamma, cp)
% ENG: engine parameters (pi_c, pi_f, Tt4_max, hPR, ENG.p0_9, ENG.p0_19, d (diameter, ft))
% THR: required thrust
% TL: tech level

% Calculate Through Mattingly's Equations
Rc = (FC.gamma - 1)/(FC.gamma) * FC.cp * 778.16; % BTU/lbm *R -> lbf. ft/ lbm * R
gc = 32.174; % Gravity Constant
a0 = sqrt(FC.gamma * Rc * gc * FC.T0);
V0 = a0*FC.M0;

% Init
Tt4_vals = linspace(ENG.Tt4_min, ENG.Tt4_max, 1000);

% Ram
TR.r = 1 + 0.5*(FC.gamma - 1)*FC.M0^2;
PR.r = TR.r ^ (FC.gamma/ (FC.gamma - 1));

if FC.M0<=1
    eta_r = 1;
else
    eta_r = 1 - 0.075*(FC.M0 - 1)^1.35;
end

% Diffuser
TR.d = 1;
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
cp_t = [];
gamma_t = [];
for i = 1:length(Tt4_vals)
    [cp_t(i), gamma_t(i)] = cp_val(Tt4_vals(i), 'combustion');
end
TR.lambda = cp_t.*Tt4_vals ./ (FC.cp * FC.T0);
Rt = (gamma_t - 1)./gamma_t .* cp_t * 778.16; % BTU/lbm*R -> lbf ft/ lbm * R

% Initial Values
TR.t = ENG.tau_t;
TR.fR = ENG.tau_f;
TR.f = TR.fR;
PR.t = TR.t .^(gamma_t ./ ((gamma_t - 1).* TL.e_t));
TR.cR = ENG.pi_c ^ ((FC.gamma - 1)/(FC.gamma * TL.e_c)); % reference tau_c

% Start ITeration Here

TR.c = 1 + (TR.cR - 1) .* (Tt4_vals./(FC.T0 * TR.r))./(ENG.Tt4_sl / ENG.Tt0).*(TR.fR ./ TR.f);
PR.c = TR.c .^ (FC.gamma * TL.e_c / (FC.gamma - 1));
PR.f = TR.f .^ (FC.gamma * TL.e_f / (FC.gamma - 1));

% Fan Conditions
PR.s19_0 = PR.r * PR.d * PR.f * TL.pi_n; % Pt19/P0
if PR.s19_0 < (0.5*(FC.gamma + 1))^(FC.gamma / (FC.gamma - 1))
    PR.s19 = PR.s19_0; %Pt19/P19
else
    PR.s19 = (0.5*(FC.gamma + 1))^(FC.gamma / (FC.gamma - 1));
end

M19 = sqrt(2./(FC.gamma - 1) .* (PR.s19 .^ ((FC.gamma- 1)./FC.gamma) - 1));

% Nozzle Conditions
PR.s9_0 = PR.r .* PR.d .* PR.c .* TL.pi_b .* PR.t .* TL.pi_n; % Pt9/P0
PR.s9 = zeros(size(gamma_t));
for i = 1:length(gamma_t)
    if PR.s9_0 < (0.5*(gamma_t(i) + 1))^(gamma_t(i)/(gamma_t(i) - 1))
        PR.s9(i) = PR.s9_0;
    else
        PR.s9(i) = (0.5*(gamma_t(i) + 1))^(gamma_t(i)/(gamma_t(i) - 1));
    end
end

M9 = sqrt(2./(gamma_t - 1) .* ((PR.s9.^((gamma_t - 1)./gamma_t)) - 1));

BPR = ENG.BPR .* (ENG.pi_c./PR.c).*sqrt((TR.lambda./(TR.r * TR.f))./(ENG.tau_lam ./ (ENG.tau_r * ENG.tau_f))) .* MFP(M19, FC.gamma, Rc) ./ MFP(ENG.M19, FC.gamma, Rc);
TR.f = 1 + (1 - TR.t)./(1 - ENG.tau_t) .* (TR.lambda ./ TR.r)./(ENG.tau_lam / ENG.tau_r) .* (1 + ENG.BPR)./(1 + BPR) .* (ENG.tau_f -1);
TR.t = 1 - (1 - PR.t.^((gamma_t - 1).*TL.e_t./gamma_t));
PR.t = PR.t .* sqrt(TR.t./ENG.tau_t) * MFP(ENG.M9, ENG.gam_t, ENG.Rt) ./ MFP(M9, gamma_t, Rt);
% end

if max(abs(TR.t - ENG.tau_t))/ENG.tau_t > 0.0001
   fprintf('Residual not met\n'); 
end

mdot = ENG.mdot_sl .* (1 + BPR)./(1 + ENG.BPR).*(FC.P0 .* PR.r .* PR.d .* PR.f .* PR.c)./(ENG.P0 * ENG.pi_f * ENG.pi_c) .* sqrt(ENG.Tt4_sl./Tt4_vals);
% Fsp_req = thrust./mdot;

Tt3 = FC.T0 * TR.r * TR.d * TR.c;
f_ratio = (cp_t.*Tt4_vals - (FC.cp * Tt3))./(TL.eta_b * ENG.hPR - cp_t.*Tt4_vals);
TR.s9_0 = (Tt4_vals .* TR.t / FC.T0)./(PR.s9.^((gamma_t - 1)./gamma_t)); % T9/T0

V9_a0 = M9.*sqrt((gamma_t .* Rt) ./ (FC.gamma * Rc) .* TR.s9_0);

% Fan Nozzle Conditions
TR.s19_0 = (TR.r * TR.f) ./ ( PR.s19 .^ ((FC.gamma- 1)./FC.gamma));
V19_a0 = M19.*sqrt(TR.s19_0);

% Thrust and SFC Performance
Fsp = 1./(1 + BPR) .* (a0 / gc) .* ((1 + f_ratio) .* V9_a0 - FC.M0 + (1 + f_ratio).*(Rt.* TR.s9_0 ./ (Rc .* V9_a0) .* (1 - ENG.p0_9)./FC.gamma)) + BPR/(1 + BPR) .* (a0 / gc) .*(V19_a0 - FC.M0 +(TR.s19_0 ./ (V19_a0) .* (1 - ENG.p0_19)./FC.gamma));
thr_avail = Fsp .* mdot;
SFC = f_ratio ./ ((1 + BPR) .* Fsp) * 3600;

Tts4 = spline(thr_avail, Tt4_vals, thrust);
f_ratio_T = spline(thr_avail, f_ratio, thrust);
SFC_T = spline(thr_avail, SFC, thrust);


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

    function MFPval = MFP(Ma, gam, Rgas)
       MFPval = sqrt(gam .* gc ./ Rgas) .* Ma .* (1 + 0.5.*(gam -1).*(Ma.^2)).^(-(gam + 1)./(2 * (gam - 1)));
    end

end