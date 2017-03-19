function [Tts4, f_ratio, SFC] = throttle_calcs(FC, ENG, TL, thrust)

% FC: flight conditions (P0, T0, M0)
% PR: pressure ratios (compressor, turbine, ram)
% TR: temperature ratios (compressor, turbine, ram)
% THR: required thrust
% TL: tech level

% Ram
TR.r = 1 + 0.5*(FC.gamma - 1)*FC.M0^2;
PR.r = TR.r ^ (FC.gamma/ (FC.gamma - 1));

% Diffuser
TR.d = 1;
if FC.M0<=1
    eta_r = 1;
else
    eta_r = 1 - 0.075*(FC.M0 - 1)^1.35;
end

PR.d = TL.d_max * eta_r;

% Compressor
TR.c = ENG.pi_c ^ ((FC.gamma - 1)/(FC.gamma * TL.e_c));

Tt.s3 = FC.T0 * TR.r * TR.d * TR.c;
Pt.s3 = FC.P0 * PR.r * PR.d * ENG.pi_c;


% Turbine Conditions
Tt4_vals = linspace(Tt.s3, ENG.Tt4_max, 1000);
cp_t = [];
gamma_t = [];
for i = 1:length(Tt4_vals)
    [cp_t(i), gamma_t(i)] = cp_val(Tt4_vals(i), 'combustion');
end
TR.lambda = cp_t.*Tt4_vals ./ (FC.gamma * FC.T0);

f_ratio = (cp_t.*Tt4_vals - (FC.cp * Tt.s3))./(TL.eta_b * ENG.hPR - cp_t(i).*(Tt4_vals));

TR.t = 1 - (TR.r ./ TR.lambda) .* (1./(TL.eta_m .* (1 + f_ratio))) .* (TR.c - 1 + ENG.BPR * (ENG.pi_f - 1));
PR.t = TR.t .^(gamma_t ./ ((gamma_t - 1).* TL.e_t));

% Nozzle Conditions
PR.s9 = ENG.p0_9 * PR.r * PR.d * PR.c * TL.pi_b * PR.t * TL.pi_n;
TR.s9_0 = 


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