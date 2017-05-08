rho = 0.002378; %slugs/ft^3
Lv = TAIL.Lopt+3.19; %ft
Vw = 40*1.68781; %Cross wind 40 knots -> ft/s
U1 = 1.1*211; %ft/s approach speed
Vt = sqrt(Vw^2+U1^2); % Total velocity
Ss = 1.02*(150*8 + TAIL.Sv); %Side area projected

% Find center of the fuselage and vertical tail
Xca = ((150*8)*(150/2)+ TAIL.Sv*(150 - TAIL.cv + (TAIL.cv/2))) / (150*8 + TAIL.Sv);

% Sadraey's method of Xcg
Xcg = 150 - Lv - 0.75 * TAIL.cv;

% OVERWRITE
Xcg = 91.6;

% Distance between aircraft cg and side view center
dc =  Xca - Xcg;

% Aircraft side force produced by cross wind
Cdy = 0.6; % Chosen from sadraey.
Fw = 0.5*rho*Vw^2*Ss*Cdy; % Force

% Sideslip angle
beta = atan(Vw/U1);

% Calculate Cnbeta
Kf1 = 0.75; % pg 742 Sadraey
Clalphav = 2*pi; % Subsonic CLalpha curve approximation
etav = 0.96; % Vertical Tail Efficiency
Cnbeta = Kf1*Clalphav*(1-0)*etav*(Lv*TAIL.Sv)/(WING.geom.span*WING.geom.S_area);
Kf2 = 1.35; % pg 742 Sadraey
Cybeta = -Kf2*Clalphav*(1-0)*etav*(TAIL.Sv)/(WING.geom.S_area);

%% Start from step 13
taur = 0.51; % Rudder effectiveness from Figure 12.12 --> 0.4 Cr/C
br = 1; % Replace with actual span of rudder control
bv = 1; % Replace with actual span of rudder
Cydeltar =  Clalphav*etav*taur*(br/bv)*(TAIL.Sv/WING.geom.S_area);
Cndeltar = -Clalphav*Vv*etav*taur*(br/bv);
% Coeffcients in static condition before start of control input
Cno = 0;
Cyo = 0;

%% Solve iteratively
sigma = 0; % Initial guess
sigmaold = 1;
deltarold = 1;
deltar = 0;
res = 1e-10;
i = 1;
count = 1e3;
while (sigmaold - sigma) >= res && (deltarold - deltar) >= res && i < count
    deltarold = deltar;
    sigmaold = sigma;
    f = @(deltar) (0.5*rho*Vt^2*WING.geom.S_area*WING.geom.span*(Cno+Cnbeta*(beta-sigma)+Cndeltar*deltar)+Fw*dc*cos(sigma) );
    deltar = fsolve(f,sigma);
% deltar = ((-Fw*dc*cos(sigma)/(0.5*rho*Vt^2*WING.geom.S_area*WING.geom.span)) - ...
%     (Cno+Cnbeta*(beta-sigma)))/Cndeltar;
g = @(sigma) (0.5*rho*Vw^2*Ss*Cdy - 0.5*rho*Vt^2*WING.geom.S_area*(Cyo+Cybeta*(beta-sigma)+Cydeltar*deltar));
sigma = fsolve(g,deltar);
% sigma =  - (0.5*rho*Vw^2*Ss*Cdy/0.5*rho*Vt^2*WING.geom.S_area - (Cyo + Cydeltar*deltar))/Cybeta - beta;
i = i +1;
end

deltarcrosswind = deltar*57.3;
% 0.5*rho*Vt^2*WING.geom.S_area*WING.geom.span*(Cno+Cnbeta*(beta-sigma)+Cndeltar*deltar)+Fw*dc*cos(sigma) = 0
% 0.5*rho*Vw^2*Ss*Cdy = 0.5*rho*Vt^2*WING.geom.S_area*(Cyo+Cybeta*(beta-sigma)+Cydeltar*deltar)

%% Check asymmetric thrust condition one engine out - OEO
Vmc = 0.8;
Vs = constraints.Vstall*Vmc;
T = 21900;
D_f = 6; % ft
D_e = 54/12; % ft
deltar = (T*(D_f+D_e/2))/(-0.5*rho*Vs^2*WING.geom.S_area*WING.geom.span*Cndeltar);
deltar = deltar*57.3;
if deltar < 30
    fprintf('Rudder Meets One Engine Out Condition \n')
end

if deltar < deltarcrosswind
    fprintf('Crosswind Control Max Deflection = %0.2f degrees\n', deltarcrosswind);
elseif deltar > deltarcrosswind 
    fprintf('Engine Out Condition Max Control Deflection = %0.2f degrees\n', deltar);
end

%% Rudder Geometry
cr = 0.3*TAIL.cv;
br = TAIL.bv*(br/bv);
Sr = br*cr;
fprintf('The rudder span is: %0.2f ft\n', br)
fprintf('The rudder chord is: %0.2f ft\n', cr)
fprintf('The rudder area is: %0.2f ft^2\n', Sr)