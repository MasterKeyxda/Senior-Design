function [ varargout ] = AltTable( in, type )
%AltTable Uses altitude data for 0 to 100 kft for the Mattingly altitude
%tables. It takes in a numeric value and a type string describing what the
%numeric value cooresponds to on the table. type can be 'h', 'delta', or
%'sigma' cooresponding to an altitude, pressure ratio, and density ratio
%respectively. An input of 'h' will give an output of the cooresponding
%ratios [delta, theta, sigma, a/a_ref]. An input of 'delta' or 'sigma' will
%output the altitude at which that value occurs by interpolation. Theta and
%a cannot be used to find altitude as they aren't strictly monotonic. 

%% Data From Tables
alt_dat = 0:2:100;
delta_dat = [10000 9298 8637 8014 7429 6878 6362 5877 5422 4997 4599 4227 ...
    3880 3557 3256 2975 2715 2474 2250 2044 1858 1688 1534 1394 1267 1151 ...
    1046 950.7 864.0 785.2 713.7 648.6 589.5 535.8 487.1 442.9 402.8 366.5 ...
    333.6 303.6 276.5 251.8 229.4 209.1 190.6 173.8 158.5 144.6 132.0 120.4 ...
    110.0]/10000;
theta_dat = [10000 9863 9725 9588 9450 9313 9175 9038 8901 8763 8626 8489 8352 ...
    8215 8077 7940 7803 7666 7529 7519*ones(1,14) 7520 7542 7563 7584 7605 ...
    7626 7647 7668 7689 7710 7731 7752 7772 7793 7814 7835 7856 7877]/10000;
sigma_dat = [10000 9428 8881 8359 7861 7386 6933 6502 6092 5702 5332 4980 4646 ...
    4330 4031 3747 3480 3227 2988 2719 2471 2245 2040 1854 1685 1531 1391 1265 ...
    1149 1044 949.2 862.7 784.1 712.5 645.9 585.7 531.2 482.0 437.4 397.1 360.6 ...
    327.6 297.6 270.5 245.9 223.6 203.4 185.1 168.4 153.3 139.6]/10000;

if type == 'h'
    delta = interp1(alt_dat, delta_dat, in);
    theta = interp1(alt_dat, theta_dat, in);
    rho = interp1(alt_dat, sigma_dat, in);
    a = sqrt(theta);
    varargout = {delta, theta, rho, a};
elseif strcmp('delta',type) == 1
    alt = interp1(delta_dat, alt_dat, in);
    varargout = {alt};
elseif strcmp('sigma',type) == 1
    alt = interp1(sigma_dat, alt_dat, in);
    varargout = {alt};
end

