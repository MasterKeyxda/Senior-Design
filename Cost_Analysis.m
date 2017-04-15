%% Cost Analysis from Raymer CH 18 pg 501-517
% Using Variable Notation from Raymer
% How many cycles or flights will be conducted per year?
% Work on Maintanence
% Flyaway cost is approximately 49 million 2025 dollars-comparable to Gulfstream
% asking price between 32 and 52 million in 2017. -Flyaway cost confirmed.
%% Inflation Rate 1986-2025 US Bureau of Statistics and Calculator
% 1986-2025
Inf = 2.65; % inflation rate calcuated from avg 2.2 rate from 2017-2025
%% Inital Parameters
% Must Run Iter Weights before 
% (*) symbol indicates things that need to be changed
WTO = Wt.WTO; % lbf
W_E = Wt.WE; % lbf
V = Wt.fuel.V_max_cr * 0.592484; % maximum velocity in knots
Q = 300; % Production Quantity based on paper "Case for Small SS AC"
FTA = 2; % Number of flight test aircraft (typically 2-6)
N_e = 3; % Number of Engines
C_e = 1000000; % engine cost* NEEDS UPDATE and move down later
N_Eng = Q*N_e; % total production quantity times number of engines per aircraft
M_max = 2.3; % engine maximum mach number
T_max = 4384; % engine max thrust (lbf)* Check for sea level
T_4 = 2000; % turbine inlet temperature (Rankine)
C_avionics = 2000*Wt.WTO; % avionics cost

% Hourly Rates in 1986 constant dollars
R_E = 59.10; % engineering
R_T = 60.70; % tooling
R_Q = 55.40; % quality control
R_M = 50.10; % manufacturing
%% Modified DAPCA IV Cost Model pg 508-517
%% RDTE and Production Costs
% Cost in constant 1986 dollars not including inflation

% Fudge Factors for different materials for hour calculations
% (aluminum: 1.0, graphite-epoxy = 1.5-2.0, fiberglass 1.1-1.2, steel
% 1.5-2.0, titanium 1.7-2.2)

FF = 1; % currently for aluminum aircraft

% Engineering Hours 
H_E = FF*(4.86*Wt.WE^0.777*V^0.894*Q^0.163);

% Tooling Hours Eq 18.2
H_T = FF*(5.99*Wt.WE^0.777*V^0.696*Q^0.263);

% Mfg (manufacturing) Hours Eq 18.3
H_M = FF*(7.37*Wt.WE^0.82*V^0.484*Q^0.641);

% QC (quality control) Hours Eq 18.4
H_Q = 0.133*H_M; % (mfg hours)for other than cargo plane

% Development support costs Eq 18.5
C_Dev = 45.42*Wt.WE^0.630*V^1.3;

% Flight Test Cost Eq 18.6
C_F = 1243.03*Wt.WE^0.325*V^0.822*FTA^1.21;

% Mfg Materials Cost Eq 18.7
C_Mat = 11*Wt.WE^0.921*V^0.621*Q^0.799;

% Eng Production Cost Eq 18.8
C_Eng = 1.15*(1548*(0.043*T_max+ 243.25*M_max + 0.969*T_4 - 2228));

% RDT&E (research development test evaluation) + flyaway COST Eq 18.9 
RDTE = Inf*((H_E*R_E + H_T*R_T + H_Q*R_Q + H_M*R_M + C_Dev + C_F + C_Mat + ...
    C_Eng*N_Eng + C_avionics)/Q);

fprintf('The total aircraft RDTE and Flyaway cost per aircraft is $%0.0f in 2025  \n',RDTE)

Flyaway = Inf*(H_M*R_M+H_T*R_T + C_F + C_Mat +  C_Eng*N_Eng + C_avionics)/Q;
% Flyaway cost per aircraft
fprintf('The Flyaway cost per aircraft is $%0.0f in 2025  \n',Flyaway)
%% Learning Curve
H_1 = H_E + H_T + H_M + H_Q; % total number of labor hours to produce aircraft
Q_range = linspace(1,1000,1000); % range of quantity production
x = 0.678; % 80% learning curve
H = H_1*(Q_range).^(x-1); % hour learning curve equation Fig 18.1
figure()
Qmin = 250; Qmax = 400; % Production Quantity range based on paper "Case for Small SS AC"
plot(Q_range,H,'LineWidth',1) % plot curve
hold on
% Plotting maximum, minimum, and estimated production quantities
plot(linspace(1,Qmin,100),H(Qmin)*ones(1,100),'r','LineWidth',0.5)
plot(linspace(1,Qmax,100),H(Qmax)*ones(1,100),'c','LineWidth',1)
plot(linspace(1,Q,100),H(Q)*ones(1,100),'g','LineWidth',1)
plot(Qmin*ones(1,100),linspace(0,H(Qmin),100),'r','LineWidth',0.5)
plot(Qmax*ones(1,100),linspace(0,H(Qmax),100),'c','LineWidth',0.5)
plot(Q*ones(1,100),linspace(0,H(Q),100),'g','LineWidth',0.5)
xlabel('Production Quantity (AC)')
ylabel('Production Labor Hours per Aircraft (Hrs)')
title('Aircraft Production Learning Curve')
legend('Production Learning Curve','Minimum Production Estimate', 'Maximum Production Estimate',...
    'Estimated Production Run')
hold off
%% Operations and Maintanence Costs
% Two man crew cost Eq 18.10
C_crew = 35*(V*WTO/10^5)^0.3 + 84;

% Table 8.1 pg 511
MMH_FH = 4.5; % Maintenance man hours per flight hour business jet (avg)
FH_YR = 1250; % Flight hour per year per aircraft


% Calculating aircraft cost one less engine


C_a =Inf*((H_M*R_M+H_T*R_T + C_F + C_Mat +  C_Eng*Q*(N_e-1) + C_avionics)/Q);
%fprintf('The aircraft cost minus one engine is %0.1f \n',C_a)

% material cost per flight hour Eq 18.12
MC_FH = Inf*(3.3*(C_a/10^6) + 7.04 + (58*(C_e/10^6) - 13)*N_e);
%fprintf('The material cost per flight hour is %0.1f \n',MC_FH)

% material cost per cycle Eq 18.13
MC_cycle = Inf*(4.0*(C_a/10^6) + 4.6 + (7.5*(C_e/10^6) + 2.8)*N_e);
%fprintf('The material cost per cycle is %0.1f \n',MC_cycle)

% "cycles" are the number of flights
% The total materials cost is the cost per flight hour times the 
% flight hours per year, plus the cost per cycle times the cycles per year.
%% Bar Graph of All Costs

figure()
Cost = Inf*[H_E*R_E, H_T*R_T, H_Q*R_Q, H_M*R_M,C_Dev,C_F,C_Eng*N_Eng,C_avionics,C_Mat]/Q;  % array of costs'
x = 1:length(Cost); % number of costs to show
barh(x,Cost) % bar graph
Labels = {'Engineering', 'Tooling', 'Quality Assurance',...
    'Manufacturing','Development', 'Flight Testing','Turbofan Engines',...
    'Avionics','Materials'};

set(gca, 'YTick', 1:length(Cost), 'YTickLabel', Labels);
xlabel('Cost per Airplane (2025 $)')
xlim([0,30e6]);
for i1=1:numel(Cost)
    text(Cost(i1),x(i1),num2str(Cost(i1),'%0.1d'),...
               'HorizontalAlignment','left',...
               'VerticalAlignment','middle')
end
title('RDTE and Flyaway Cost Breakdown')