%% V-n Diagram 

% Date Created: 2/26/17
% Constructs a V-n diagram for a FAR 25 Certified Airplane
% Equivalent Airspeed (KEAS) is plotted on the horizontal axis while load
% factor, n, is plotted on the vertical axis. 
% Source: Airplane Aerodynamics and Performance by Lan and Roskam
% Equations come from Section 12.4

% Load in variables from other scripts
if ~exist('Wt', 'var')
   clc;
   clear;
   load('aircraft_vars.mat'); 
   close all;
end
%% Maneuever V-n Diagram

% Info pulled from other scripts
% MTOW (lbf)
wingLoading = designPoint(1); % from Constraint_Plots.m script
% rhoSL = 0.002378; % density slugs/ft^3

% Positive +1g stall speed 
CLMaxPos = 1.75; % max positive coefficient of lift
% CDatCL = 0.314; % drag coefficient at CLmax
% CNMax = sqrt((CLMax^2) + (CDatCL^2)); % Roskam part 5, eqn 4.5 (flaps up)
CNMaxPos = 1.1*CLMaxPos; % eqn 12.17, preliminary value; update once CD is known at CLmax
Vs1 = sqrt(2*(wingLoading) / (atm.rho_sl * CNMaxPos));  % ft/s, eqn 12.15
Vs1 = Vs1 * 0.592484; % convert to knots
fprintf('The positive +1g stall speed is %0.2f KEAS \n', Vs1)

% Negative 1g stall speed
CLMaxNeg = -1; % max negative coefficient of lift
CNMaxNeg = 1.1*CLMaxNeg; % eqn 12.24, preliminary value; update once CD is known at CLmax
nNeg = -1:0.05:0; % negative load factor range, eqn 12.41
Vsneg = sqrt((2*nNeg*wingLoading) / (atm.rho_sl * CNMaxNeg));
Vsneg = Vsneg * 0.592484; % convert to knots
fprintf('The negative +1g stall speed is %0.2f KEAS \n', max(Vsneg))

% Positive Design limit load factor, nlim
nlimPos = 2.1 + (24000 / (Wt.WTO + 10000));
% nlimPos must be greater than or equal to 2.5
if nlimPos < 2.5
    nlimPos = 2.5; 
end

% Ultimate Load Factor
nUlt = nlimPos * 1.5; 
fprintf('The max limit load factor is %0.2f \n', nlimPos)
fprintf('The max ultimate load factor is % 0.2f \n', nUlt)

% Positive stall curve for load factor ranging from n = 0 to nlimPos
nPos = 0:0.05:nlimPos; % positive load factor range
Vspos = sqrt((2*nPos*wingLoading) / (atm.rho_sl * CNMaxPos)); % ft/s, eqn 12.18
Vspos = Vspos * 0.592484; % convert to knots

% Design Maneuvering Speed, VA
VA = Vs1*sqrt(nlimPos); % knots, eqn 12.40 
fprintf('The design maneuvering speed is %0.2f KEAS \n', VA)

% Design Cruise Speed, VC
M = req.cr_M0(1); % cruise Mach Number
[~,~,sigma,aRatio] = AltTable(atm.alt,'h'); % speed of sound ratio
aSound = 1116.5 * aRatio; % speed of sound at cruise alt, ft/s
Vtrue = M * aSound; % cruise speed (TAS), ft/s
Vtrue = Vtrue * 0.592484; % convert true airspeed to knots (KTAS)
VC = sqrt(sigma) * Vtrue; % cruise speed (KEAS)
fprintf('The design cruise speed is %0.2f KEAS \n', VC )

% Design Diving Speed, VD
VD = 1.25*VC; % knots (KEAS)
fprintf('The design diving speed is %0.2f KEAS \n', VD)
% Linear Curve from VC to VD, eqn 12.42
m = 0 - min(nNeg) / (VD - VC); % slope of line
x = VC:1:VD; % vector from VC to VD
y = m*x - 4.9934; % equation of line

% Plot 
figure()
ylim([min(nNeg)-0.5, nlimPos+0.5])
xlabel('Equivalent Airspeed (KEAS)','FontSize', 14)
ylabel('Load factor, n', 'FontSize', 14)
title('V-n Diagram', 'FontSize', 14)
hold on; 
% Plot positive stall speed curve
plot(Vspos, nPos,'k', 'LineWidth',2.5);

% Plot negative stall speed curve
vline1 = plot(Vsneg, nNeg,'k','LineWidth',2.5);

% Plot horizontal line from VA to VD at max positive n
plot([VA, VD], [nlimPos, nlimPos],'k','LineWidth',2.5) 

% Plot horizontal line from max(Vs1neg) to VC at max negative n
plot([max(Vsneg), VC], [min(nNeg), min(nNeg)],'k', 'LineWidth',2.5) 

% Plot linear Curve from VC to VD
plot(x,y,'k','LineWidth',2.5) 

% Plot Vertical line from VD to nlimPos
plot([VD,VD],[0,nlimPos],'k', 'LineWidth',2.5)

% Plot horizontal line from 0 to VD
plot([0, VD],[0,0],'k')

% Plot vertical line of 1g stall speed, Vs1
stallLine = plot([Vs1,Vs1],[0,1],'r:', 'LineWidth',1.6);

% Plot vertical line of design maneuvering speed, VA to nlimPos
ManvLine = plot([VA,VA], [0,nlimPos],'b:', 'LineWidth',1.6);

% Plot vertical line of design cruising speed, VC from nNeg to nlimPos
CruiseLine = plot([VC, VC], [min(nNeg), nlimPos],'g:', 'LineWidth',1.6);

%% Gust Diagram

e = 0.8; 
%CLalpha = WING.biconvex(2).Clalpha / (1 + (WING.biconvex(2).Clalpha / pi * e * WING.geom.AR));
% CLalpha = 1/deg; 
CLalpha = (WING.CFD.SUP.CL(3) - WING.CFD.SUP.CL(2))/(WING.CFD.SUP.alpha(3) - WING.CFD.SUP.alpha(2));
% Convert to 1/rad
CLalpha = CLalpha * (180/pi); 
h = (atm.alt * 10^3); % cruise altitude in ft
g = 32.17; % acceleration of gravity ft/s^2
atm.rho_VN = atm.rho_sl * sig_rho;
% Mass Ratio (eqn 12.31)
massRatio = (2 * wingLoading) / (atm.rho_VN * WING.geom.MAC * g * CLalpha);

% Gust Alleviation Factor (eqn 12.30)
Kg = (massRatio^1.03) / (6.9 + (massRatio^1.03));

% Gust Lines
% VB gust line
%VB = VC - 43; % Maximum Gust Intensity Speed (KEAS)
%UdeB = 84.67 - (0.000933 * h); % derived gust velocity
%UdeB = UdeB * 0.592484; % convert to knots
%VGustB = 0:1:VB; % gust line varies from 0 to VB
%VGustB = VGustB * 1.68781; % convert to ft/s
%jB = (Kg * UdeB * VGustB * CLalpha) / (498 * wingLoading); % limit factor
%nlimBPos = 1 + jB; % positive slope gust factor line
%nlimBNeg = 1 - jB; % negative slope gust factor line
%plot(VGustB*0.592484, nlimBPos, 'b')
%plot(VGustB*0.592484, nlimBNeg, 'b')

% VC gust line
UdeC = 66.67 - (0.000833 * h); % derived gust velocity (ft/s)
% UdeC = UdeC * 0.592484; % convert to knots
VGustC = 0:1:VC; % gust line varies from 0 to VC
VGustC = VGustC * 1.68781; % ft/s
jC = (Kg * UdeC * VGustC * CLalpha) / (498 * wingLoading); % limit factor
nlimCPos = 1 + jC; % positive slope gust factor line
nlimCNeg = 1 - jC; % negative slope gust factor line
plot(VGustC*0.592484, nlimCPos, 'k--')
plot(VGustC*0.592484, nlimCNeg, 'k--')

% VD gust line
UdeD = 33.35 - (0.000417 * h); % derived gust velocity
%UdeD = UdeD * 0.592484; % convert to knots
VGustD = 0:1:VD; % gust line varies from 0 to VD
VGustD = VGustD * 1.68781;
jD = (Kg * UdeD * VGustD * CLalpha) / (498 * wingLoading); % limit factor
nlimDPos = 1 + jD; % positive slope gust factor line
nlimDNeg = 1 - jD; % negative slope gust factor line
vline2 = plot(VGustD*0.592484, nlimDPos, 'k--');
plot(VGustD*0.592484, nlimDNeg, 'k--')

% Extend line from VC to VD
m = (min(nlimDNeg) - min(nlimCNeg)) / (VD - VC);
x = VC:1:VD; % vector from VC to VD
y = -m*x+4.02; % equation of line
% Plot linear Curve from VC to VD for gust diagram
plot(x,y,'k--') 

y = m*x-2.01; % equation of line
plot(x,y,'k--') 

% Legend
% legend([vline1, vline2, stallLine, ManvLine, CruiseLine],'Manuever Lines',...
%     'Gust Lines', '1g-stall speed','Manuever Speed',...
%     'Cruise Speed','Location', 'Northwest')

legend([vline1, vline2],'Manuever Lines',...
    'Gust Lines','Location', 'Northwest')
