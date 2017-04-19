function t_val = foil_t_get(fname, x_c)

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

t_val = abs(spline(top(:,1), top(:,2), x_c)) + abs(spline(bot(:,1), bot(:,2), x_c));

end