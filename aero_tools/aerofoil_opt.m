%% AIRFOIL OPTIMIZATION

clear;
clc;

% Get Wing Data
% Supersonic_Aircraft_Design;
close all;

clearvars -except WING
% clc;

%% Airfoil Setup Controls

t_c = 0.0525; % global thickness ratio
taper = 0.45;
fname = 'biconvex';
ftype = 'dat';
dx = 1/50;

wt = 0;
sim_type = 'ansys';

% chord_len
chord = 21.84 * 12; % ft -> in

%% Modify the Biconvex
airfoil_design(fname, ftype, dx, 1.0, t_c, wt, sim_type);

%% Modify the HSNLF
fid = fopen('BACNLF.dat', 'r');
tline = fgets(fid);
coords = [];
i = 1;
while ischar(tline)
    if ~isempty(str2num(tline))
        coords(i,:) = str2num(tline);
%         fprintf([tline '\n']);
        i = i+1;
    end
    
    tline = fgets(fid);
end

fclose(fid);

% Get Upper Camber
camber_upper = 0.5*(spline(coords(1:67, 1),coords(1:67, 2), coords(1:67, 1)) + spline(coords(68:139, 1),coords(68:139, 2), coords(1:67, 1)));
diff_up = coords(1:67,2) - camber_upper;
figure();
plot(coords(1:67, 1), camber_upper);
hold on;
plot(coords(:,1), coords(:,2), '*');
axis equal


% Get Lower Camber
camber_lower = 0.5*(spline(coords(1:67, 1),coords(1:67, 2), coords(68:139, 1)) + spline(coords(68:139, 1),coords(68:139, 2), coords(68:139, 1)));
diff_low = coords(68:139,2) - camber_lower;
plot(coords(68:139, 1), camber_lower, 'o');
hold on;

t_up = max(coords(1:67, 2) - camber_upper)*2.0;
t_low = max(abs(coords(68:139, 2) - camber_lower))*2.0;

disp(t_up);
disp(t_low);

coords = coords .* chord;
% coords(68:139, 2) = coords(68:139, 2) .* chord;

plot(coords(:,1), coords(:,2));
legend('Camber Line (Upper)', 'Original Airfoil', 'Camber Line (Lower)', '7% Airfoil', 'eastoutside');

% % write file for openvsp
% fid = fopen('BACNLF_7.dat', 'w');
% 
% fprintf(fid, 'DEMO GEOM AIRFOIL FILE\n');
% fprintf(fid, 'BACNLF 7 Percent\n');
% fprintf(fid, '0\tSym Flag\n');
% fprintf(fid, '%i\t Num Pnts Upper\n', 67);
% fprintf(fid, '%i\t Num Pnts Lower\n', 72);
% 
% for i = 1:66
%     fprintf(fid, '%3.7f\t%4.7f\n', coords(67-i,1), coords(67-i,2));
% end
% fprintf(fid, '\n');
% 
% for i = 68:139
%     fprintf(fid, '%3.7f\t%4.7f\n', coords(i,1), coords(i,2));
% end


% fclose(fid);

% write file for selig format
fid = fopen('BACNLF_real_selig.dat', 'w');
fprintf(fid, 'BACNLF w/ Root Chord Percent\n\n');
for i = 1:length(coords)
    fprintf(fid, '%3.7f\t%4.7f\t%i\n', coords(i,1), coords(i,2), 0.0);
    if i == 0.5*length(coords(:,1))
        fprintf(fid, '\n');
    end
end

fclose(fid);
