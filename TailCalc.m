% clc
% clear all
% close all
%% FILE META
% File is now updated from AMAT calculation script into function calculations for senior design
function TAIL = TailCalc(alpha, Vh, Vv, Wt, rho, vel, Df, Kc, S, AR, Cmaf, sweepWing, taperh, cglocAC, vtail, M)
    %% Inputs:
%     alpha = 0; % deg AOA
%     Vh = 0.4; % horizontal tail volumen coefficient
%     Vv = 0.03; % vertical tail volumen coefficient
%     Wt = 4.5; % lbf Average Weight
%     rho = 0.002378; % slugs/ft^3 density
%     vel = 54; % ft/s Velocity
    cbar = sqrt(S/AR);%8/12; % ft MAC
    b = sqrt(AR*S);%48.125/12; % ft Wingspan
%     Df = 4/12; % ft diameter of fuselage
%     Kc = 1.2; % correction factor
%     S = b*cbar; % ft^2 planform area
%     AR = b^2/S; % aspect ratio wing
%     Cmaf = -0.181; % CM - 0.305
%     sweepWing = 6; % wing sweep angle (degrees)
%     cglocAC = -(3/8)/12; % ft cg location in front or behind AC Wing
    
    
    Lopt = Kc*((4*cbar*S*Vh)/(pi*Df))^0.5; % ft optimum tail moment arm
    Lfuselage = Lopt / 0.3; % fuselage length
    
    Sh = (cbar*S*Vh)/Lopt; % ft^2 horizontal tail area
    CL = (2*Wt)/(rho*vel^2*S); % coefficient of lift
    Cm0wf = Cmaf*((AR*cosd(sweepWing)^2) / (AR + 2*cosd(sweepWing))); % CM wing and fuselage
    
    Xcg = 0.25*cbar - cglocAC; % cg location from leading edge
    Xcgbar = Xcg / cbar; % nondimensional cg
    CLh = (Cm0wf + CL*(Xcgbar - 0.25))/Vh; % CL horizontal tail
    ARh = 2/3 * AR; % horizontal tail aspect ratio
    % ARh = 4;
    CL_alphah = 4/sqrt(M^2-1); % horizontal tail Cl_alpha (per radian) #Linearized Supersonic Flow
    CL_alpha = CL_alphah / ( 1 + (CL_alphah / (pi*ARh))); % tail lift curve slope
    Alphah = CLh/CL_alphah; % Horizontal Tail Alpha
    Alphah = Alphah*57.3;
    Alphah = Alphah - .1; % adjust angle of horizontal tail
    eps0 = 2 * CL / (pi*AR);
    % CL_alphaw = 3.763; % from xflr
    CL_alphaw = 1.72; % per radian
    deps = 2 * CL_alphaw / (pi*AR); % rad/rad
    eps = eps0 + deps*alpha/57.3;

    ih = Alphah + eps; % horizontal tail incidence angle
    ch = sqrt(Sh/ARh);
    bh = Sh/ch;
    cRooth = Sh / (bh * 0.5*(1 + taperh)); % root chord of horizontal tail
    cTiph = taperh * cRooth; % tip chord of horizontal tail
    effh = 0.9; % tail efficiency
    Cm_alpha = CL_alphaw*(Xcgbar - 0.25) - CL_alphah*effh*(Sh/S)*(Lopt/cbar - Xcgbar)*(1-deps);


    Sv = (b*S*Vv)/(Lopt);
    cv = sqrt(Sv/ARh);
    bv = sqrt(Sv * ARh);

    TotArea = Sv + Sh;
    vb = sqrt(ARh*TotArea);
    vbone = vb/2;
    % vbone = sqrt(ARh*(TotArea/2));

    cb = TotArea/vb;
    Angle = atand(sqrt(Sv/Sh));
    Lifth = 0.5*rho*vel^2*Sh*CLh; % total lift by total horizontal surface area
    Lifth = Lifth / 2; % split the lift into two parts for each arm
    LiftVtail = Lifth/cosd(Angle);
    CL_Vtail = (2*LiftVtail)/(rho*vel^2*(cb*vbone));
    AlphaVt = CL_Vtail/CL_alphah * 57.3;
    AlphaVt = AlphaVt + eps; % correction
    % CL_alphah = .11;
    % CL_alphaw = .1; %deg/deg
    % deps=0.4315;
    t = (CL_alphah/CL_alphaw)*effh*(Sh/S)*(1-deps);
    x_acwf = .25; % 25% MAC
    x_ach = Lopt/cbar;% RINA HOW DID YOU GET THIS NUMBER? CAN YOU CALCULATE THIS FROM OTHER PARAMETERS?
    x_ac = (x_acwf+t*x_ach)/(1+t);
    SM = x_ac-Xcgbar;

    % cma = - CL_alphaw*(x_ac-Xcgbar)
    
    %% Function Readouts
    fprintf('Tail Moment Arm = %0.3f ft\n', Lopt);
    fprintf('Horizontal Tail Area = %0.3f \n',Sh);
    fprintf('Horizontal Tail Chord = %0.3f ft \n', ch);
    fprintf('Horizontal Tail Root Chord = %0.3f \n', cRooth); 
    fprintf('Horizontal Tail Tip Chord = %0.3f \n', cTiph); 
    fprintf('Horizontal Tail Span = %0.3f ft \n', bh);
    fprintf('Horizontal Tail Angle = %0.3f deg \n',ih);
    fprintf('Cm_alpha = %0.3f  \n', Cm_alpha);
    fprintf('Vertical Tail Area = %0.3f ft^2 \n',Sv)
    fprintf('Vertical Tail Chord = %0.3f ft \n', cv);
    fprintf('Verical Tail Span = %0.3f ft \n', bv);
    fprintf('Total Tail Area = %0.3f ft^2\n',(Sh+Sv));
    if strcmp(vtail, 'yes')
        fprintf('One arm of V Tail')
        fprintf('V Tail Span = %.2f ft \n', vbone);
        fprintf('V Tail Chord = %.2f ft \n', cb);
        fprintf('V Tail Angle = %.2f deg \n',Angle)
        fprintf('Angle of V-tail with the x axis = %.2f deg\n',AlphaVt)
    end
    fprintf('X_AC_bar = %.2f\n', x_ac)
    fprintf('CG Location = %.2f ft behind LE.\n', Xcg)
    fprintf('Distance LE Wing to LE Tail = %0.2f ft\n', Xcg+Lopt-0.25*cv);
    fprintf('Static Margin (Percent) = %.2f \n', SM*100)
    
    % Save Important Performance Params
    TAIL.Lopt = Lopt;
    TAIL.Sh = Sh;
    TAIL.ch = ch;
    TAIL.Cr_h = cRooth;
    TAIL.bh = bh;
    TAIL.Cm_alpha = Cm_alpha;
    TAIL.Sv = Sv;
    TAIL.cv = cv;
    TAIL.Area_tot = (Sh+Sv);
    TAIL.vbone = vbone;
    TAIL.cb = cb;
    TAIL.Angle = Angle;
    TAIL.AlphaVt = AlphaVt;
    TAIL.x_ac = x_ac;
    TAIL.Xcg = Xcg;
    TAIL.SM = SM;
    TAIL.Xcgbar = Xcgbar;
    TAIL.hAngle = ih;
    TAIL.bv = bv;
end