clc;
clear;
close all;

%% Initialize Study Parameters

pass.test_var = 6:18; % limit to 18 to optimize crew weight

pass.ranges = zeros(size(pass.test_var));
pass.fuel_eff = pass.ranges;
pass.weight = pass.fuel_eff;
pass.fuel_weight = pass.weight;
pass.thrust_req = pass.weight;

% Constants
    % engine
    % crew count
    % 
    
% Variables
    % fuel efficiency (lbm fuel / (nmi * passenger))
    % weight/thrust
    % range
    
% Process
    % Set passenger count -> 
    % run iter_weights
    % run TO_TD_perf (maybe not)
    % do drag analysis (?)
    % run Payload_Range
    
    % save range data
    % save fuel efficiency data (get specific numbers for lbm/pass*hr)
    % save actual thrust needed
    
%% Run scripts

pass.ct = 1;
while pass.ct <= length(pass.test_var)

    pass.num = pass.test_var(pass.ct);
    iter_weights;
    
    % Save Variables
    pass.weight(pass.ct) = Wt.WTO;
    pass.fuel_weight(pass.ct) = Wt.fuel.w_tot;
    pass.thrust_req(pass.ct) = constraints.req_Thr;
    %Wt.fuel.w_max = (1/req.f_eff)*(Wt.pld.n_pass * req.range)*Wt.fuel.reserve_ratio;
    
    Drag_Analysis;
    Payload_Range;
    
    pass.ranges(pass.ct) = Range.cr_climb.array(1);
    
    clearvars -except pass
    disp(pass.num);
    pass.ct = pass.ct + 1;
end

%% Post-Process

close all;
pass.fuel_eff = (pass.num .* pass.ranges)./pass.fuel_weight;

figure();
plot(pass.test_var, pass.fuel_weight);

figure();
plot(pass.test_var, pass.fuel_weight);

figure();
plot(pass.test_var, pass.weight);

figure();
plot(pass.test_var, pass.ranges);
hold on;
plot(spline(pass.ranges, pass.test_var, 4000), 4000, 'o');

range_opt = floor(spline(pass.ranges, pass.test_var, 4000));

figure();
plot(pass.test_var, pass.fuel_eff);
hold on;
plot(spline(pass.fuel_eff, pass.test_var, 1.0), 1.0, 'o');

fuelEff_opt = floor(spline(pass.fuel_eff, pass.test_var, 1.0));

pass.opt_num = min(range_opt, fuelEff_opt);
fprintf('Optimum Passenger Ct: %i\n', pass.opt_num);