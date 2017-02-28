%% V-n Diagram 

% Date Created: 2/26/17
% Constructs a V-n diagram for a FAR 25 Certified Airplane
% Equivalent Airspeed (KEAS) is plotted on the horizontal axis while load
% factor, n, is plotted on the vertical axis. 
% Source: Airplane Aerodynamics and Performance by Lan and Roskam
% Equations come from Section 12.4

%% Maneuever V-n Diagram

WTO = 102703; % MTOW (lbf)
wingLoading = 108.5; % wing loading design point (lb/ft^2)
rhoSL = 0.002378; % density slugs/ft^3

% Positive +1g stall speed 
CLMaxPos = 2.2; % max positive coefficient of lift
% CDatCL = 0.314; % drag coefficient at CLmax
% CNMax = sqrt((CLMax^2) + (CDatCL^2)); % Roskam part 5, eqn 4.5 (flaps up)
CNMaxPos = 1.1*CLMaxPos; % eqn 12.17, preliminary value; update once CD is known at CLmax
Vs1 = sqrt(2*(wingLoading) / (rhoSL * CNMaxPos));  % ft/s, eqn 12.15
Vs1 = Vs1 * 0.592484; % convert to knots

% Negative 1g stall speed
CLMaxNeg = -1; % max negative coefficient of lift
CNMaxNeg = 1.1*CLMaxNeg; % eqn 12.24, preliminary value; update once CD is known at CLmax
nNeg = -1:0.05:0; % negative load factor range, eqn 12.41
Vsneg = sqrt((2*nNeg*wingLoading) / (rhoSL * CNMaxNeg));
Vsneg = Vsneg * 0.592484; % convert to knots

% Positive Design limit load factor, nlim
nlimPos = 2.1 + (24000 / (WTO + 10000));
% nlimPos must be greater than or equal to 2.5
if nlimPos < 2.5
    nlimPos = 2.5; 
end

% Positive stall curve for load factor ranging from n = 0 to nlimPos
nPos = 0:0.05:nlimPos; % positive load factor range
Vspos = sqrt((2*nPos*wingLoading) / (rhoSL * CNMaxPos)); % ft/s, eqn 12.18
Vspos = Vspos * 0.592484; % convert to knots

% Design Maneuvering Speed, VA
VA = Vs1*sqrt(nlimPos); % knots, eqn 12.40 

% Design Cruise Speed, VC
M = 1.6; % cruise Mach Number
h = 42; % cruise altitude (kft)
[~,~,sigma,aRatio] = AltTable(h,'h'); % speed of sound ratio
aSound = 1116.5 * aRatio; % speed of sound at cruise alt, ft/s
Vtrue = M * aSound; % cruise speed (TAS), ft/s
VC = sqrt(sigma) * Vtrue; % cruise speed (EAS), ft/s
VC = VC * 0.592484; % convert to knots (KEAS)

% Design Diving Speed, VD
VD = 1.25*VC; % knots (KEAS)

% Linear Curve from VC to VD, eqn 12.42
m = 0 - min(nNeg) / (VD - VC); % slope of line
x = VC:1:VD; % vector from VC to VD
y = m*x - 4.9934; % equation of line

% Plot 
figure()
ylim([min(nNeg)-0.5, nlimPos+0.5])
xlabel('Speed, V (KEAS)')
ylabel('Load factor, n')
title('V-n Diagram (Maneuver and Gust)')
hold on; 
% Plot positive stall speed curve
plot(Vspos, nPos,'b') 

% Plot negative stall speed curve
plot(Vsneg, nNeg,'b') 

% Plot horizontal line from VA to VD at max positive n
plot([VA, VD], [nlimPos, nlimPos],'b') 

% Plot horizontal line from max(Vs1neg) to VC at max negative n
plot([max(Vsneg), VC], [min(nNeg), min(nNeg)],'b') 

% Plot linear Curve from VC to VD
plot(x,y,'b') 

% Plot Vertical line from VD to nlimPos
plot([VD,VD],[0,nlimPos],'b')

% Plot horizontal line from 0 to VD
plot([0, VD],[0,0],'k')

% Plot vertical line of 1g stall speed, Vs1
plot([Vs1,Vs1],[0,1],'k--')

% Plot vertical line of design maneuvering speed, VA to nlimPos
plot([VA,VA], [0,nlimPos],'k--')

% Plot vertical line of design cruising speed, VC from nNeg to nlimPos
plot([VC, VC], [min(nNeg), nlimPos],'k--')
%% Gust Diagram

MAC = 16; % mean aerodynamic chord, ft
CLalpha = 4.87; % lift-curve slope, 1/rad
h = 42000; % cruise altitude ft
g = 32.17; % acceleration of gravity ft/s^2

% Mass Ratio (eqn 12.31)
massRatio = (2 * wingLoading) / (rhoSL * MAC * g * CLalpha);

% Gust Alleviation Factor (eqn 12.30)
Kg = (massRatio^1.03) / (6.9 + (massRatio^1.03));

% VB gust line
VB = VC - 43; % Maximum Gust Intensity Speed (KEAS)
UdeB = 84.67 - (0.000933 * h); % derived gust velocity
UdeB = UdeB * 0.592484; % convert to knots
VGustB = 0:1:VB; % gust line varies from 0 to VB
jB = (Kg * UdeB * VGustB * CLalpha) / (498 * wingLoading); % limit factor
nlimBPos = 1 + jB; % positive slope gust factor line
nlimBNeg = 1 - jB; % negative slope gust factor line
%plot(VGustB, nlimBPos)
%plot(VGustB, nlimBNeg)

% VC gust line
UdeC = 66.67 - (0.000833 * h); % derived gust velocity (ft/s)
UdeC = UdeC * 0.592484; % convert to knots
VGustC = 0:1:VC; % gust line varies from 0 to VC
jC = (Kg * UdeC * VGustC * CLalpha) / (498 * wingLoading); % limit factor
nlimCPos = 1 + jC; % positive slope gust factor line
nlimCNeg = 1 - jC; % negative slope gust factor line
%plot(VGustC, nlimCPos)
%plot(VGustC, nlimCNeg)

% VD gust line
UdeD = 33.35 - (0.000417 * h); % derived gust velocity
UdeD = UdeD * 0.592484; % convert to knots
VGustD = 0:1:VD; % gust line varies from 0 to VD
jD = (Kg * UdeD * VGustD * CLalpha) / (498 * wingLoading); % limit factor
nlimDPos = 1 + jD; % positive slope gust factor line
nlimDNeg = 1 - jD; % negative slope gust factor line
%plot(VGustD, nlimDPos)
%plot(VGustD, nlimDNeg)