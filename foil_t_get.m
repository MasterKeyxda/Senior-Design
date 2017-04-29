function [tu_val, tl_val] = foil_t_get(fname, x_c)
load('aircraft_vars.mat', 'WING');

if strcmp(fname, 'biconvex')
    tau_u = WING.biconvex(2).tc_u;
    tau_l = WING.biconvex(2).tc_l;
    
    tu_val = 2.*tau_u.*x_c.*(1 - x_c);
    tl_val = 2.*tau_l.*x_c.*(1-x_c);
    
else
    fid = fopen([fname '.dat'], 'r');
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

    nose = find(coords(:,1) == min(coords(:,1)));

    top = coords(1:nose, :);
    bot = coords((nose+1):end, :);

figure();plot(top(:,1), top(:,2));
hold on;
plot(bot(:,1), bot(:,2));
legend('Top', 'Bottom');
axis equal;

    tu_val = abs(spline(top(:,1), top(:,2), x_c));
    tl_val = abs(spline(bot(:,1), bot(:,2), x_c));
end

end